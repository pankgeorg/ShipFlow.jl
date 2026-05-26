#!/usr/bin/env julia
#
# F2: two-way coupling amplification trial. Place the rudder inside
# the rotor race (within 0.5·R of the disk plane) and compare the
# rudder side-force / hull-side-force response between WL_TWOWAY=0
# (rudder sees freestream only) and WL_TWOWAY=1 (rudder sees the
# local flow via `trilinear_inflow`).
#
# In the original rudder_trial_sinusoidal.jl the rudder sits 2·R aft
# of the disk plane, well outside the race; the inflow correction was
# only a ~0.6 % effect. Inside the race the rotor jet roughly doubles
# the axial flow, so dynamic pressure quadruples and the inflow path
# should matter materially.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)

using WaterLily, VoF, ShipShapes, LiftingSurfaces
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf, Statistics

const NX = parse(Int, get(ENV, "WL_NX", "128"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "32"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "80"))
const δ_DEG  = parse(Float32, get(ENV, "WL_DELTA", "10"))
const RUDDER_OFFSET = parse(Float64, get(ENV, "WL_RUD_OFF", "0.5"))  # × R

const L_c = 36f0; const B_c = 8f0; const T_c = 5f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = 0.30f0
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ/2)
const hull_xc = Float32(NX/5); const hull_yc = Float32(NY/2); const hull_zc = H_w_c

const R_prop  = 3.0f0
const prop_xc = Float32(hull_xc + L_c/2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c/2)
const J_op = 0.32
const Ω = Float64(π) * U∞ / (J_op * R_prop)

# In-race rudder placement.
const rud_xc = Float32(prop_xc + RUDDER_OFFSET * R_prop)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c/2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

@printf "=== F2: rudder-in-race two-way amplification ===\n"
@printf "  Grid       = %d × %d × %d\n" NX NY NZ
@printf "  Steps      = %d   δ = %.1f°   rud_off = %.2f·R\n" NSTEPS δ_DEG RUDDER_OFFSET
@printf "  rud_xc - prop_xc = %.2f   (R_prop = %.2f)\n" (rud_xc - prop_xc) R_prop
flush(stdout)

rotor  = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)
rudder = Rudder(; chord=4.0, span=5.0, ns=12, nc=6)

r_rot = rotor_forces(rotor, U∞, Ω)
Sref_r = π * R_prop^2
thrust = abs(r_rot.CT * 0.5 * U∞^2 * Sref_r)
torque = r_rot.CQ * 0.5 * U∞^2 * Sref_r * R_prop
@printf "  Rotor at J=%.2f: |CT|=%.3f thrust=%.2f |CQ|=%.3f torque=%.2f\n" J_op abs(r_rot.CT) thrust abs(r_rot.CQ) torque

function vof_pois_ctor(vof)
    (flow) -> begin
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
end

function run_case(label::String; twoway::Bool)
    @info "Running $label"
    vof = VoFFlow((NX, NY, NZ);
        α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
    sim = WaterLily.Simulation((NX, NY, NZ),
        (U∞, 0f0, 0f0), L_c;
        T=Float32, ν=vof.ν,
        g=(i, x, t) -> i == 3 ? -G_c : 0f0,
        Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
        pois_ctor=vof_pois_ctor(vof), U=U∞,
    )

    CL_hist  = Float32[]
    Fhy_hist = Float32[]
    Vlocal_hist = Float32[]
    t0 = time()
    for step in 1:NSTEPS
        function combo_udf(flow, t; kwargs...)
            smear_force!(flow.f, SVector(Float32(thrust), 0f0, 0f0),
                         SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
            smear_torque!(flow.f, Float32(torque),
                          SVector(prop_xc, prop_yc, prop_zc),
                          SVector(1f0, 0f0, 0f0), Float32(R_prop * 0.7);
                          N=8, ε=2.0f0)
            inflow_fn = twoway ? trilinear_inflow(flow.u) : nothing
            # Sample once at rudder centre for diagnostics.
            if twoway
                v = inflow_fn(SVector(Float64(rud_xc), Float64(rud_yc), Float64(rud_zc)))
                push!(Vlocal_hist, Float32(sqrt(v[1]^2 + v[2]^2 + v[3]^2)))
            else
                push!(Vlocal_hist, U∞)
            end
            r_rud = rudder_forces(rudder, deg2rad(δ_DEG), U∞; inflow=inflow_fn)
            side  = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
            drag  = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
            smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
                         SVector(rud_xc, rud_yc, rud_zc); ε=2.0f0)
            push!(CL_hist, Float32(r_rud.CL))
            return nothing
        end

        WaterLily.mom_step!(sim.flow, sim.pois; udf=combo_udf, pois_tol=1f-6, pois_itmx=50)
        step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1], perdir=(2,))
        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        push!(Fhy_hist, Float32(Fp[2] + Fv[2]))

        if step % 20 == 0
            @info @sprintf("  [%s]  step=%3d  |u|=%.3f  F_hy=%+.3f  elapsed=%.1fs",
                label, step, maximum(abs, sim.flow.u), Fhy_hist[end], time()-t0)
            flush(stdout)
        end
    end
    n = max(1, length(CL_hist) ÷ NSTEPS)
    CL_step = CL_hist[n:n:end]
    Vloc_step = Vlocal_hist[n:n:end]
    # Average over last 25% of run
    tail = (NSTEPS - NSTEPS÷4 + 1):NSTEPS
    CL_m  = mean(CL_step[tail])
    Fhy_m = mean(Fhy_hist[tail])
    Vloc_m = mean(Vloc_step[tail])
    return (; label, CL_m, Fhy_m, Vloc_m,
              gain = abs(Fhy_m) / max(abs(CL_m), 1e-6))
end

r0 = run_case("TWOWAY=0"; twoway=false)
r1 = run_case("TWOWAY=1"; twoway=true)

println()
@printf "═══════════════════════════════════════════════════════════════════════════\n"
@printf "                          TWOWAY=0         TWOWAY=1\n"
@printf "  CL_rudder              %+.4f          %+.4f\n" r0.CL_m r1.CL_m
@printf "  F_hull_y               %+.3f           %+.3f\n" r0.Fhy_m r1.Fhy_m
@printf "  |V_local| at rudder     %.3f            %.3f\n" r0.Vloc_m r1.Vloc_m
@printf "  Manoeuvring gain        %.2f            %.2f\n" r0.gain r1.gain
@printf "  Δ gain (TWOWAY=1 − 0)   %+.2f%%\n" 100*(r1.gain/r0.gain - 1)
@printf "═══════════════════════════════════════════════════════════════════════════\n"

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "rudder_in_race"))
mkpath(OUTDIR)
open(joinpath(OUTDIR, "summary.csv"), "w") do io
    println(io, "case,CL_m,Fhy_m,Vloc_m,gain")
    println(io, "twoway0,$(r0.CL_m),$(r0.Fhy_m),$(r0.Vloc_m),$(r0.gain)")
    println(io, "twoway1,$(r1.CL_m),$(r1.Fhy_m),$(r1.Vloc_m),$(r1.gain)")
end
println("Summary written to $(joinpath(OUTDIR, "summary.csv"))")
