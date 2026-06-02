#!/usr/bin/env julia
#
# Self-propulsion via J-sweep with the VLM BladedRotor on the
# Containership hull (par_frac=0.5, Cb≈0.75). Mirrors
# wigley_self_propulsion_VLM.jl exactly except for the hull SDF.
#
# Procedure:
#   for J in J_LIST:
#     Ω = π·V∞ / (J·R)
#     run NSTEPS_PER steps with the VLM rotor's CT smeared in
#     record final hull drag (averaged over last AVG_FRAC of run)
#   write (J, |CT|, thrust, drag, T-D) and find self-prop J.
#
# Hypothesis: the higher block coefficient (0.75 vs Wigley's 0.44)
# produces more pressure drag at the same Fr, so the self-propulsion
# point should shift to higher J (more thrust required).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)

using WaterLily
using VoF
using ShipShapes
using ShipShapes: StaticArrays
using LiftingSurfaces
const SVector = StaticArrays.SVector
using Printf

const NX = parse(Int, get(ENV, "WL_NX", "96"))
const NY = parse(Int, get(ENV, "WL_NY", "48"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const L_c = parse(Float64, get(ENV, "WL_L", "48"))
const B_c = parse(Float64, get(ENV, "WL_B", "10"))
const T_c = parse(Float64, get(ENV, "WL_T", "6"))
const PAR_FRAC = parse(Float64, get(ENV, "WL_PAR_FRAC", "0.5"))
const ρ_w = 10.0
const ρ_a = 1.0
const U∞  = 1.0
const Fr  = parse(Float64, get(ENV, "WL_FR", "0.25"))
const G_c = U∞^2 / (Fr^2 * L_c)
const Re  = parse(Float64, get(ENV, "WL_RE", "1000"))
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const T_NUM = Float64

const H_w_c   = NZ/2
const hull_xc = NX/3; const hull_yc = NY/2; const hull_zc = H_w_c
const prop_xc = hull_xc + L_c/2 + T_c/2
const prop_yc = hull_yc
const prop_zc = H_w_c - T_c/2
const R_prop  = parse(Float64, get(ENV, "WL_RPROP_FAC", "0.8")) * T_c / 2

# Default J range biased lower than Wigley's default (0.4..1.0). At higher
# Cb the hull drags more, so the self-prop J should shift higher, but the
# scan needs to bracket the actual zero crossing.
const J_LIST = let s = get(ENV, "WL_J_LIST", "0.25,0.30,0.35,0.40,0.50,0.70,1.00")
    [parse(Float64, x) for x in split(s, ",")]
end
const NSTEPS_PER = parse(Int, get(ENV, "WL_NSTEPS_PER", "100"))
const AVG_FRAC   = 0.25

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "containership_selfprop_VLM"))
mkpath(OUTDIR)

α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0
hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)

# One blade definition reused per J.
rotor = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.2*R_prop,
    chord=(0.25 * R_prop, 0.18 * R_prop),
    twist=(deg2rad(35.0), deg2rad(15.0)),
    ns=12, nc=4)

function run_one(J)
    Ω = π * U∞ / (J * R_prop)

    vof = VoFFlow((NX, NY, NZ);
        α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM)
    hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c,
        par_frac = PAR_FRAC, map = hull_map)

    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ),
        (T_NUM(U∞), T_NUM(0), T_NUM(0)), L_c;
        T = T_NUM, ν = VoF.viscosity(vof),
        g = (i, x, t) -> i == 3 ? T_NUM(-G_c) : T_NUM(0),
        Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
        pois_ctor = vof_pois_ctor, U = T_NUM(U∞),
    )

    # Cache the rotor's CT and CQ once per J — VLM result doesn't depend
    # on the WaterLily state (no two-way coupling in this scan).
    r_rot = rotor_forces(rotor, U∞, Ω)
    Sref  = π * R_prop^2
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
            udf = rotor_udf, pois_tol = T_NUM(1e-6), pois_itmx = 50)
        step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
        if mod(step, force_every) == 0 || step > NSTEPS_PER * (1 - AVG_FRAC)
            Fp = WaterLily.pressure_force(sim)
            Fv = WaterLily.viscous_force(sim)
            push!(drag_history, -(Float64(Fp[1]) + Float64(Fv[1])))
        end
    end
    N_avg = max(1, Int(round(length(drag_history) * AVG_FRAC)))
    drag_avg = sum(drag_history[end-N_avg+1:end]) / N_avg
    return abs(r_rot.CT), thrust, drag_avg
end

csv = open(joinpath(OUTDIR, "scan.csv"), "w")
println(csv, "J,CT,thrust,drag_avg,T-D")
@printf "Containership (par_frac=%.2f, Cb≈%.2f)  Re=%.0f  NSTEPS_PER=%d  scanning %d J values…\n" PAR_FRAC ((1+PAR_FRAC)/2) Re NSTEPS_PER length(J_LIST)
println("J     |CT|      thrust    drag    T - D")
for J in J_LIST
    CT, T_val, D_val = run_one(J)
    diff = T_val - D_val
    println(csv, "$J,$(CT),$(T_val),$(D_val),$(diff)")
    flush(csv)
    @printf "%.2f  %8.4f  %8.3f  %8.3f  %+8.3f\n" J CT T_val D_val diff
end
close(csv)

# Linear-interpolate the self-propulsion J (T - D crosses zero)
data = open(joinpath(OUTDIR, "scan.csv"), "r") do io
    lines = readlines(io)[2:end]
    [parse.(Float64, split(l, ',')) for l in lines]
end
function find_sign_change(data)
    for i in 1:length(data)-1
        da, db = data[i][5], data[i+1][5]
        if da * db < 0
            t = -da / (db - da)
            return data[i][1] + t * (data[i+1][1] - data[i][1])
        end
    end
    return nothing
end
J_sp = find_sign_change(data)
println()
if J_sp === nothing
    println("No sign change observed — extend J_LIST.")
else
    @printf "Estimated self-propulsion J ≈ %.4f\n" J_sp
end
println("Scan written to $(joinpath(OUTDIR, "scan.csv"))")
