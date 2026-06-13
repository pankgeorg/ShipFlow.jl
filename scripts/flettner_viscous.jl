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
# Spinning a cylinder in WaterLily: the SDF of a circle is rotationally
# symmetric, so a rotation `map` leaves the geometry unchanged but imposes the
# surface velocity via −d(map)/dt (BDIM reads the body velocity from the map's
# time-derivative — the same mechanism as the `rotate` maintest). With
# map(x,t) = Rot(−ωt)·(x−c) the material points spin at ω, giving a surface
# tangential speed ωR with no change to the immersed shape.
#
# Run (headless, single ω):
#   julia --project=. scripts/flettner_viscous.jl --omega 1.0
# Full sweep:
#   julia --project=. scripts/flettner_viscous.jl --sweep

using WaterLily
using WaterLily: SA          # StaticArrays' SA, re-exported through WaterLily
using Printf
using LinearAlgebra: norm

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "flettner_viscous"))
mkpath(OUTDIR)

# ----------------------------------------------------------------------------
# Simulation set-up
# ----------------------------------------------------------------------------

"""
    flettner_simulation(; ω, R_cells=32, Re=200, n_diam_x=16, n_diam_y=8)

A 2D rotating cylinder of radius `R_cells` (in cells), spinning at non-
dimensional rate `ω` (so the surface speed is `ω·R` in cell/time units, and
the tip-speed ratio relative to the freestream U=1 is `ω·R/U`). `Re = U·D/ν`
with `D = 2R`. The cylinder is placed `4D` downstream of the inlet, mid-height.
"""
function flettner_simulation(; ω::Real, R_cells::Int=32, Re::Real=200,
                              n_diam_x::Int=16, n_diam_y::Int=8)
    nx = n_diam_x * 2R_cells
    ny = n_diam_y * 2R_cells
    cx, cy = 4 * 2R_cells, ny / 2
    c  = SA[Float64(cx), Float64(cy)]
    R  = Float64(R_cells)
    U  = 1.0
    D  = 2R

    # circle SDF (body frame) and a spin map: material points rotate at ω.
    sdf(x, t) = √sum(abs2, x) - R
    function spin(x, t)
        s, cs = sincos(ω * t)
        Rot = SA[cs s; -s cs]      # Rot(−ωt) on (x−c)
        Rot * (x .- c)
    end
    body = AutoBody(sdf, spin)

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

# Surface pressure Cp(θ): sample the pressure field just outside the cylinder
# on a ring of `nθ` angles. Cp = (p − p∞)/(½ρU²); p∞≈0 in WaterLily's gauge so
# Cp ≈ 2p/(U²). Sampled at radius R+1.5 cells (first clear cell off the body).
function surface_cp(sim; nθ::Int=180)
    p = sim.flow.p
    R = sim.L / 2
    # recover the centre from the domain geometry used in flettner_simulation
    ny = size(p, 2) - 2
    cx = 4 * sim.L; cy = (ny) / 2 + 1        # +1 ghost offset (loc convention)
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
# Run one ω
# ----------------------------------------------------------------------------

function run_one(ω; R_cells::Int=32, Re::Real=200, t_end::Real=60.0, tag="")
    @info "Flettner viscous run" ω R_cells Re
    sim = flettner_simulation(; ω, R_cells, Re)
    sim_step!(sim)                          # settle BCs
    t_s = Float64[]; Cd_s = Float64[]; Cl_s = Float64[]
    step = 0
    while sim_time(sim) < t_end
        sim_step!(sim; remeasure = (ω != 0))   # rotating body must remeasure
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

    θ, Cp = surface_cp(sim)
    sfx = isempty(tag) ? @sprintf("w%.1f", ω) : tag
    open(joinpath(OUTDIR, "cp_$sfx.csv"), "w") do io
        println(io, "theta,Cp")
        for (t, c) in zip(θ, Cp); @printf(io, "%.6f,%.6f\n", t, c); end
    end
    @printf("ω=%.2f  Cl_mean=%.4f  Cd_mean=%.4f  (samples=%d)\n",
            ω, Cl_mean, Cd_mean, length(Cl_s))
    return (; ω, Cl = Cl_mean, Cd = Cd_mean, θ, Cp)
end

function run_sweep(; R_cells::Int=32, Re::Real=200, t_end::Real=60.0,
                   ωs = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5))
    rows = NamedTuple[]
    open(joinpath(OUTDIR, "cl_vs_omega.csv"), "w") do io
        println(io, "omega,Cl,Cd")
        for ω in ωs
            r = run_one(ω; R_cells, Re, t_end)
            @printf(io, "%.3f,%.6f,%.6f\n", r.ω, r.Cl, r.Cd); flush(io)
            push!(rows, (; r.ω, r.Cl, r.Cd))
        end
    end
    @info "Wrote $(joinpath(OUTDIR, "cl_vs_omega.csv"))"
    return rows
end

# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
    if "--sweep" in args
        # allow a coarser/cheaper config from the CLI
        R = ("--R" in args) ? parse(Int, args[findfirst(==("--R"), args)+1]) : 24
        Re = ("--Re" in args) ? parse(Float64, args[findfirst(==("--Re"), args)+1]) : 200.0
        te = ("--tend" in args) ? parse(Float64, args[findfirst(==("--tend"), args)+1]) : 50.0
        run_sweep(; R_cells=R, Re=Re, t_end=te)
    else
        ω = ("--omega" in args) ? parse(Float64, args[findfirst(==("--omega"), args)+1]) : 1.0
        R = ("--R" in args) ? parse(Int, args[findfirst(==("--R"), args)+1]) : 24
        Re = ("--Re" in args) ? parse(Float64, args[findfirst(==("--Re"), args)+1]) : 200.0
        te = ("--tend" in args) ? parse(Float64, args[findfirst(==("--tend"), args)+1]) : 50.0
        run_one(ω; R_cells=R, Re=Re, t_end=te)
    end
end
