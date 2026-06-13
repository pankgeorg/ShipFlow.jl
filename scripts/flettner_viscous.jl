#!/usr/bin/env julia
#
# Viscous (real-flow) rotating cylinder — the "Flettner rotor" sidebar to the
# inviscid panel/analytic solution in NavalArchitectToolbox (`flettner_panel`,
# `flettner_analytic`). A 2D WaterLily run of a smooth circular cylinder of
# radius R spinning at angular rate ω in a uniform stream V∞, sweeping the
# same ω∈{0,0.5,1,1.5,2,2.5} as the panel study. We extract the surface
# pressure Cp(θ) and the lift coefficient C_L(ω), and write them to CSV for
# the analytical/panel/viscous overlay (see RESULTS-flettner.md).
#
# Why this differs from potential flow: a real boundary layer separates and
# the spin convects the separation point, so the viscous C_L grows far more
# slowly with ω than the inviscid 4πωR²/V∞ — that *deviation* is the point of
# the comparison.
#
# --- THE FIX (2026-06-13) -------------------------------------------------
# The first cut spun the cylinder with a rotation *map* (the SDF of a circle
# is rotation-symmetric, so the map left the shape unchanged while −d/dt of
# the map gave the surface speed ωR). That polluted `pressure_force`: with a
# time-dependent map, the BDIM surface kernel `nds(body,x,t)` and the gauge
# pressure carry the rigid-body / unsteady-frame contribution, which the
# simple `F = −(fp+fv)` reduction cannot separate from the aerodynamic lift
# (it returned wrong-sign, ~10× lift).
#
# The fix is to keep the SDF *genuinely static* and impose the spin purely as
# the no-slip surface velocity via the body velocity field V. We define a
# custom `SpinningCylinder <: WaterLily.AbstractBody` whose `measure(body,x,t)`
# returns a TIME-INDEPENDENT circle (d, radial normal n) and a rigid-rotation
# velocity V = ω×(x−c). Because the geometry never moves, `nds` is constant in
# time and `pressure_force` integrates a fixed circle — the measured force is
# the aerodynamic one. This is the moving-body `measure!` convention
# (V from the body, geometry from the SDF) done the right way for a spin.
#
# Run (headless, single ω):
#   julia --project=. scripts/flettner_viscous.jl --omega 1.0 --R 32 --Re 200 --tend 60
# Full sweep:
#   julia --project=. scripts/flettner_viscous.jl --sweep --R 32 --Re 200 --tend 60

using WaterLily
using WaterLily: SA          # StaticArrays' SA, re-exported through WaterLily
using Printf
using LinearAlgebra: norm

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "flettner_viscous"))
mkpath(OUTDIR)

# ----------------------------------------------------------------------------
# A spinning cylinder as a STATIC-SDF body with a prescribed surface velocity
# ----------------------------------------------------------------------------
#
# measure(body,x,t) must return (d, n, V):
#   d = |x−c| − R                    static signed distance (circle)
#   n = (x−c)/|x−c|                  outward radial unit normal (static)
#   V = ω·(y−cy, −(x−cx))            rigid-rotation surface velocity (CW spin)
# V is time-independent, so the geometry/normal never move — the spin enters
# only as the no-slip velocity BDIM imposes inside the immersed band. This is
# what makes the recovered pressure force the aerodynamic lift.
#
# SIGN. The analytic/panel tool (`flettner_analytic`) defines a positive ω as
# producing positive (+y) lift: its Vt = 2V∞sinθ + Γ/(2πR) speeds the flow
# over the TOP of the cylinder, i.e. the surface there moves WITH the stream
# (+x at the top) → a CLOCKWISE surface rotation, bound circulation Γ<0 by the
# CCW-positive loop convention, and Magnus lift +y (L = −ρU·Γ). So to match
# that convention here the spin must be CLOCKWISE: V = ω·(dy, −dx) (a positive
# ω moves the top surface in +x). The circulation check then expects
# Γ ≈ −2πωR² (CW = negative), with the viscous value a fraction of that ideal.
struct SpinningCylinder{T} <: WaterLily.AbstractBody
    cx::T; cy::T; R::T; ω::T
end

@inline function WaterLily.measure(b::SpinningCylinder, x::AbstractVector{T}, t;
                                   fastd²=Inf) where T
    dx = x[1] - b.cx; dy = x[2] - b.cy
    r  = sqrt(dx*dx + dy*dy)
    d  = r - b.R
    d^2 > fastd² && return (d, zero(x), zero(x))
    invr = r > eps(T) ? one(T)/r : zero(T)
    n = SA[dx*invr, dy*invr]                 # outward radial normal
    # CLOCKWISE rigid rotation (top surface moves +x, with the stream) so a
    # positive ω gives +y Magnus lift — matching the analytic/panel sign.
    V = SA[b.ω*dy, -b.ω*dx]
    return (d, n, V)
end

# ----------------------------------------------------------------------------
# Simulation set-up
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Spin similarity parameter
# ----------------------------------------------------------------------------
# The physical, scale-free spin is the TIP-SPEED RATIO α = ωR/U (= surface
# speed / freestream). The inviscid Magnus lift in terms of α is, via
# Kutta–Joukowski with Γ = 2πωR² = 2π α R U and Cl = 2Γ/(U·D), D = 2R:
#       Cl_inviscid(α) = 2πα.
# The NAT panel/analytic tool runs ω∈{0,0.5,…,2.5} at R=0.5, U=1, so its tip-
# speed ratio is α = ωR/U = 0.5·ω ∈ {0,0.25,…,1.25} and its analytic
# Cl = 4πωR²/V∞ = πω = 2πα — the same curve. We therefore sweep α and report
# Cl(α) so the viscous overlay is directly comparable to the panel/analytic
# line at matching physical states (NOT at matching raw ω in mismatched units).
nat_omega_to_alpha(ω_nat; R_nat=0.5, U=1.0) = ω_nat * R_nat / U      # = 0.5 ω_nat
cl_inviscid(α) = 2π * α

"""
    flettner_simulation(; α, R_cells=32, Re=200, n_diam_x=16, n_diam_y=8)

A 2D rotating cylinder of radius `R_cells` (in cells) at tip-speed ratio
`α = ωR/U` (so the cell-units spin rate is `ω = α·U/R`). `Re = U·D/ν`,
`D = 2R`; the cylinder sits `4D` downstream of the inlet, mid-height. Built on
the static-SDF `SpinningCylinder` body.
"""
function flettner_simulation(; α::Real, R_cells::Int=32, Re::Real=200,
                              n_diam_x::Int=16, n_diam_y::Int=8)
    nx = n_diam_x * 2R_cells
    ny = n_diam_y * 2R_cells
    cx, cy = 4.0 * 2R_cells, ny / 2
    R  = Float64(R_cells)
    U  = 1.0
    D  = 2R
    ω_cell = α * U / R                      # cell-units angular rate
    body = SpinningCylinder{Float64}(cx, cy, R, ω_cell)
    Simulation((nx, ny), (U, 0.0), D; ν = U * D / Re, body = body, T = Float64)
end

# ----------------------------------------------------------------------------
# Force coefficients (lift positive in +y), normalised on the diameter
# ----------------------------------------------------------------------------

# pressure_force returns ∫ p·n dS = force ON THE FLUID; force on the body is
# its negative. C_L = 2·F_y/(ρ U² D) with ρ=U=1, D=sim.L.
function coefficients(sim)
    fp = WaterLily.pressure_force(sim)
    fv = WaterLily.viscous_force(sim)
    F  = .-(fp .+ fv)
    s  = 0.5 * sim.L * sim.U^2
    return F[1] / s, F[2] / s             # (Cd, Cl)
end

# ----------------------------------------------------------------------------
# Circulation sanity check — Γ = ∮ u·dl around a loop enclosing the cylinder.
# A positive ω here is a CLOCKWISE surface spin (top moves +x), which drives a
# clockwise bound circulation → Γ < 0 by the CCW-positive loop convention.
# The inviscid ideal is Γ = −2πωR² = −2π α R U (cell units). Checking Γ before
# trusting pressure_force isolates the flow solve from the force integration;
# the viscous Γ should be a fraction of the ideal (boundary-layer slip
# deficit). Sampled on a square loop of half-size `h` cells.
# ----------------------------------------------------------------------------
function circulation(sim; rfac::Real=1.6)
    u = sim.flow.u
    R = sim.L / 2
    ny = size(u, 2) - 2
    cx = 4.0 * sim.L + 1                       # +1 ghost offset (loc convention)
    cy = ny / 2 + 1
    # nearest-node staggered velocity components. u[I,1] is the x-velocity on
    # the x-face at the low side of cell I; u[I,2] the y-velocity. For a
    # circulation loop the half-cell stagger is a sub-grid correction we
    # neglect (the loop is many cells across).
    @inline ux(xi, yi) = (i = clamp(round(Int, xi), 2, size(u,1)-1);
                          j = clamp(round(Int, yi), 2, size(u,2)-1); u[i, j, 1])
    @inline uy(xi, yi) = (i = clamp(round(Int, xi), 2, size(u,1)-1);
                          j = clamp(round(Int, yi), 2, size(u,2)-1); u[i, j, 2])
    # Γ = ∮ u·dl on a circle of radius rfac·R, tangential (CCW positive):
    #   t̂ = (−sinθ, cosθ),  dl = r dθ
    rr = rfac * R; ns = 720; Γ = 0.0; dθ = 2π / ns
    for k in 0:ns-1
        θ = (k + 0.5) * dθ
        s, cθ = sincos(θ)
        xi = cx + rr*cθ; yi = cy + rr*s
        ut = -ux(xi, yi) * s + uy(xi, yi) * cθ     # tangential (CCW) component
        Γ += ut * rr * dθ
    end
    return Γ
end

# Surface pressure Cp(θ): sample the pressure field just outside the cylinder
# on a ring of `nθ` angles. Cp = (p − p∞)/(½ρU²); p∞≈0 in WaterLily's gauge so
# Cp ≈ 2p/(U²). Sampled at radius R+1.5 cells (first clear cell off the body).
function surface_cp(sim; nθ::Int=180)
    p = sim.flow.p
    R = sim.L / 2
    ny = size(p, 2) - 2
    cx = 4.0 * sim.L + 1; cy = ny / 2 + 1        # +1 ghost offset (loc convention)
    rr = R + 1.5
    θs = range(0, 2π; length=nθ+1)[1:nθ]
    Cp = Float64[]
    U = sim.U
    for θ in θs
        xi = cx + rr*cos(θ); yi = cy + rr*sin(θ)
        i = clamp(round(Int, xi), 2, size(p,1)-1)
        j = clamp(round(Int, yi), 2, size(p,2)-1)
        push!(Cp, 2 * p[i, j] / U^2)
    end
    return collect(θs), Cp
end

# ----------------------------------------------------------------------------
# Run one tip-speed ratio α
# ----------------------------------------------------------------------------

function run_one(α; R_cells::Int=32, Re::Real=200, t_end::Real=60.0, tag="")
    @info "Flettner viscous run" α R_cells Re
    sim = flettner_simulation(; α, R_cells, Re)
    sim_step!(sim)                          # settle BCs
    t_s = Float64[]; Cd_s = Float64[]; Cl_s = Float64[]
    step = 0
    while sim_time(sim) < t_end
        # static SDF + static V ⇒ no remeasure needed (geometry never moves).
        sim_step!(sim; remeasure = false)
        step += 1
        if step % 5 == 0
            Cd, Cl = coefficients(sim)
            push!(t_s, sim_time(sim)); push!(Cd_s, Cd); push!(Cl_s, Cl)
        end
    end
    # mean over the last half (statistically steady)
    half = max(1, length(Cl_s) ÷ 2)
    Cl_mean = sum(@view Cl_s[half:end]) / (length(Cl_s) - half + 1)
    Cd_mean = sum(@view Cd_s[half:end]) / (length(Cd_s) - half + 1)

    # circulation sanity check: ideal Γ = −2πωR² = −2π α R U (cell units, U=1)
    R = Float64(R_cells)
    Γ_meas = circulation(sim)
    Γ_ideal = -2π * α * R * sim.U
    Cl_inv = cl_inviscid(α)

    θ, Cp = surface_cp(sim)
    sfx = isempty(tag) ? @sprintf("a%.2f", α) : tag
    open(joinpath(OUTDIR, "cp_$sfx.csv"), "w") do io
        println(io, "theta,Cp")
        for (t, c) in zip(θ, Cp); @printf(io, "%.6f,%.6f\n", t, c); end
    end
    @printf("α=%.3f  Cl=%.4f (inviscid %.4f)  Cd=%.4f  Γ=%.1f (ideal %.1f)  (samples=%d)\n",
            α, Cl_mean, Cl_inv, Cd_mean, Γ_meas, Γ_ideal, length(Cl_s))
    return (; α, Cl = Cl_mean, Cl_inv, Cd = Cd_mean, Γ = Γ_meas, Γ_ideal, θ, Cp)
end

function run_sweep(; R_cells::Int=32, Re::Real=200, t_end::Real=60.0,
                   αs = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25))   # = 0.5·{0,0.5,1,1.5,2,2.5}
    rows = NamedTuple[]
    open(joinpath(OUTDIR, "cl_vs_alpha.csv"), "w") do io
        println(io, "alpha,Cl,Cl_inviscid,Cd,Gamma,Gamma_ideal")
        for α in αs
            r = run_one(α; R_cells, Re, t_end)
            @printf(io, "%.4f,%.6f,%.6f,%.6f,%.4f,%.4f\n",
                    r.α, r.Cl, r.Cl_inv, r.Cd, r.Γ, r.Γ_ideal); flush(io)
            push!(rows, (; r.α, r.Cl, r.Cl_inv, r.Cd, r.Γ, r.Γ_ideal))
        end
    end
    @info "Wrote $(joinpath(OUTDIR, "cl_vs_alpha.csv"))"
    return rows
end

# ----------------------------------------------------------------------------
# CLI — drive by tip-speed ratio α=ωR/U (the physical, scale-free spin).
# `--omega-nat` accepts the NAT panel/analytic ω (R=0.5) and converts to α.
# ----------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
    R  = ("--R"    in args) ? parse(Int,     args[findfirst(==("--R"),    args)+1]) : 32
    Re = ("--Re"   in args) ? parse(Float64, args[findfirst(==("--Re"),   args)+1]) : 200.0
    te = ("--tend" in args) ? parse(Float64, args[findfirst(==("--tend"), args)+1]) : 60.0
    if "--sweep" in args
        run_sweep(; R_cells=R, Re=Re, t_end=te)
    else
        α = if "--alpha" in args
            parse(Float64, args[findfirst(==("--alpha"), args)+1])
        elseif "--omega-nat" in args
            nat_omega_to_alpha(parse(Float64, args[findfirst(==("--omega-nat"), args)+1]))
        else
            0.5
        end
        run_one(α; R_cells=R, Re=Re, t_end=te)
    end
end
