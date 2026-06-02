#!/usr/bin/env julia
#
# I1: Containership Cb sweep. Vary par_frac (and hence Cb) across
# {0.30, 0.50, 0.70, 0.85}, finding the self-propulsion J for each.
# Generalises F1 to a one-parameter family.
#
# Cb relation:  Cb ≈ (1 + par_frac) / 2  →  Cb ∈ {0.65, 0.75, 0.85, 0.93}.
# Hypothesis:   J_self decreases monotonically as Cb grows
#               (more drag → more thrust → lower J).

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
using Printf

const NX = 96; const NY = 48; const NZ = 48
const L_c = 48.0; const B_c = 10.0; const T_c = 6.0
const ρ_w = 10.0; const ρ_a = 1.0
const U∞ = 1.0
const Fr = 0.40
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 1000.0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = NZ/2
const hull_xc = NX/3; const hull_yc = NY/2; const hull_zc = H_w_c
const prop_xc = hull_xc + L_c/2 + T_c/2
const prop_yc = hull_yc
const prop_zc = H_w_c - T_c/2
const R_prop  = 0.8 * T_c / 2

const NSTEPS_PER = 100
const AVG_FRAC   = 0.25
# J grid per par_frac. We bracket the expected range; larger par_frac =>
# more drag => lower J_self, so leaner par_frac uses a higher-biased grid.
const J_GRIDS = Dict(
    0.30 => [0.20, 0.25, 0.30, 0.35, 0.45],
    0.50 => [0.15, 0.175, 0.20, 0.25, 0.30],
    0.70 => [0.12, 0.15, 0.175, 0.20, 0.25],
    0.85 => [0.10, 0.12, 0.15, 0.175, 0.20],
)
const PAR_LIST = [0.30, 0.50, 0.70, 0.85]

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "cb_sweep_Fr40"))
mkpath(OUTDIR)

α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0
hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)

rotor = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.2*R_prop,
    chord=(0.25 * R_prop, 0.18 * R_prop),
    twist=(deg2rad(35.0), deg2rad(15.0)),
    ns=12, nc=4)

function run_one(par_frac, J)
    Ω = π * U∞ / (J * R_prop)
    vof = VoFFlow((NX, NY, NZ);
        α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float64)
    hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c,
        par_frac = par_frac, map = hull_map)
    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ),
        (U∞, 0.0, 0.0), L_c;
        T = Float64, ν = VoF.viscosity(vof),
        g = (i, x, t) -> i == 3 ? -G_c : 0.0,
        Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
        pois_ctor = vof_pois_ctor, U = U∞,
    )
    r_rot = rotor_forces(rotor, U∞, Ω)
    Sref = π * R_prop^2
    thrust = abs(r_rot.CT) * 0.5 * U∞^2 * Sref
    torque = r_rot.CQ * 0.5 * U∞^2 * Sref * R_prop
    function rotor_udf(flow, t; kwargs...)
        smear_force!(flow.f, SVector(thrust, 0.0, 0.0),
                     SVector(prop_xc, prop_yc, prop_zc); ε=2.5)
        smear_torque!(flow.f, torque,
                      SVector(prop_xc, prop_yc, prop_zc),
                      SVector(1.0, 0.0, 0.0), 0.7 * R_prop;
                      N=8, ε=2.0)
        return nothing
    end
    drag_history = Float64[]
    force_every = max(1, NSTEPS_PER ÷ 25)
    for step in 1:NSTEPS_PER
        WaterLily.mom_step!(sim.flow, sim.pois;
            udf = rotor_udf, pois_tol = 1e-6, pois_itmx = 50)
        step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
        if mod(step, force_every) == 0 || step > NSTEPS_PER * (1 - AVG_FRAC)
            Fp = WaterLily.pressure_force(sim)
            Fv = WaterLily.viscous_force(sim)
            push!(drag_history, -(Fp[1] + Fv[1]))
        end
    end
    N_avg = max(1, Int(round(length(drag_history) * AVG_FRAC)))
    drag_avg = sum(drag_history[end-N_avg+1:end]) / N_avg
    return abs(r_rot.CT), thrust, drag_avg
end

function find_sign_change(scan)
    for i in 1:length(scan)-1
        da = scan[i].T - scan[i].D
        db = scan[i+1].T - scan[i+1].D
        if da * db < 0
            t = -da / (db - da)
            return scan[i].J + t * (scan[i+1].J - scan[i].J)
        end
    end
    return nothing
end

results = []
csv = open(joinpath(OUTDIR, "scans.csv"), "w")
println(csv, "par_frac,Cb,J,CT,thrust,drag,T-D")

for par_frac in PAR_LIST
    Cb = (1 + par_frac) / 2
    @printf "=== par_frac=%.2f  Cb=%.2f ===\n" par_frac Cb
    scan = []
    for J in J_GRIDS[par_frac]
        CT, T_val, D_val = run_one(par_frac, J)
        push!(scan, (J=J, CT=CT, T=T_val, D=D_val))
        @printf "  J=%.3f  CT=%.3f  T=%.2f  D=%.2f  T-D=%+.2f\n" J CT T_val D_val (T_val-D_val)
        println(csv, "$par_frac,$Cb,$J,$CT,$T_val,$D_val,$(T_val-D_val)")
        flush(csv)
    end
    J_sp = find_sign_change(scan)
    push!(results, (par_frac=par_frac, Cb=Cb, J_self=J_sp, scan=scan))
    @printf "  J_self ≈ %s\n\n" (J_sp === nothing ? "no bracket" : @sprintf("%.4f", J_sp))
end
close(csv)

println("=" ^ 70)
@printf "  %-9s  %-7s  %-9s\n" "par_frac" "Cb" "J_self"
println("=" ^ 70)
for r in results
    @printf "  %-9.2f  %-7.2f  %s\n" r.par_frac r.Cb (
        r.J_self === nothing ? "no bracket" : @sprintf("%.4f", r.J_self))
end
println("=" ^ 70)
println("Detail in $(joinpath(OUTDIR, "scans.csv"))")
