#!/usr/bin/env julia
#
# Self-propulsion via parameter scan.  The PI controller (sibling
# script wigley_self_propulsion.jl) hits transient balance but the
# system oscillates around it because the propeller wake modifies
# the apparent hull drag — see RESULTS-headline.md.
#
# This script instead runs the simulation to (approximate) steady
# state at several fixed C_T values, records the resulting drag and
# thrust, and emits a CSV so we can interpolate the self-propulsion
# point analytically.
#
# Procedure:
#   for C_T in {0.1, 0.2, …, 1.0}:
#     run NSTEPS_PER mom_steps with that fixed thrust
#     record final drag (smoothed over last 20% of run)
#   write (C_T, T, D, T-D) to runs/wigley_selfprop_scan/scan.csv

using WaterLily
using VoF
using ShipShapes
using ShipShapes: StaticArrays
using Turbulence
using Propellers
const SVector = StaticArrays.SVector
using Printf

# Same physics as wigley_self_propulsion.jl
const NX = parse(Int, get(ENV, "WL_NX", "96"))
const NY = parse(Int, get(ENV, "WL_NY", "48"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const L_c = parse(Float64, get(ENV, "WL_L", "48"))
const B_c = parse(Float64, get(ENV, "WL_B", "10"))
const T_c = parse(Float64, get(ENV, "WL_T", "6"))
const ρ_w = 10.0
const ρ_a = 1.0
const U∞  = 1.0
const Fr  = parse(Float64, get(ENV, "WL_FR", "0.25"))
const G_c = U∞^2 / (Fr^2 * L_c)
const Re  = parse(Float64, get(ENV, "WL_RE", "1000"))
const ν_w_c = U∞ * L_c / Re
const ν_a_c = ν_w_c * 18
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * ν_a_c
const T_NUM = Float64

const H_w_c   = NZ/2
const hull_xc = NX/3; const hull_yc = NY/2; const hull_zc = H_w_c
const prop_xc = hull_xc + L_c/2 + T_c/2
const R_prop  = parse(Float64, get(ENV, "WL_RPROP_FAC", "0.8")) * T_c / 2
const W_prop  = 1.5

const CT_LIST = let s = get(ENV, "WL_CT_LIST", "0.0,0.1,0.2,0.3,0.4,0.5,0.7,1.0")
    [parse(Float64, x) for x in split(s, ",")]
end
const NSTEPS_PER = parse(Int, get(ENV, "WL_NSTEPS_PER", "80"))
const AVG_FRAC   = 0.25   # average drag over last quarter of run

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_selfprop_scan"))
mkpath(OUTDIR)

α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0
hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)

function run_one(CT)
    A = π * R_prop^2
    thrust = 0.5 * 1.0 * A * U∞^2 * CT

    vof = VoFFlow((NX, NY, NZ);
        α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM)
    turb = WALE((NX, NY, NZ); Cw = T_NUM(0.5), ν₀ = T_NUM(0))
    disk = ActuatorDisk(
        center = SVector(prop_xc, NY/2, H_w_c - T_c/2),
        axis   = SVector(1.0, 0.0, 0.0),
        R = R_prop, w = W_prop, thrust = thrust,
    )
    hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)

    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ),
        (T_NUM(U∞), T_NUM(0), T_NUM(0)), L_c;
        T = T_NUM, ν = turb.ν,
        g = (i, x, t) -> i == 3 ? T_NUM(-G_c) : T_NUM(0),
        Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
        pois_ctor = vof_pois_ctor, U = T_NUM(U∞),
    )
    disk_udf = (flow, t; kwargs...) -> disk(flow, t; kwargs...)
    drag_history = Float64[]
    # Compute force only every FORCE_EVERY steps. pressure_force +
    # viscous_force are ~70% of step cost, but for an averaged scan
    # we only need a handful of samples over the last quarter.
    force_every = max(1, NSTEPS_PER ÷ 25)
    for step in 1:NSTEPS_PER
        WaterLily.mom_step!(sim.flow, sim.pois;
            udf = disk_udf, pois_tol = T_NUM(1e-6), pois_itmx = 50)
        step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
        update_νt!(turb, sim.flow.u, vof.ν)
        if mod(step, force_every) == 0 || step > NSTEPS_PER * (1 - AVG_FRAC)
            Fp = WaterLily.pressure_force(sim)
            Fv = WaterLily.viscous_force(sim)
            push!(drag_history, -(Float64(Fp[1]) + Float64(Fv[1])))
        end
    end
    # Average drag over the last AVG_FRAC of the (sampled) run
    N_avg = max(1, Int(round(length(drag_history) * AVG_FRAC)))
    drag_avg = sum(drag_history[end-N_avg+1:end]) / N_avg
    return thrust, drag_avg
end

csv = open(joinpath(OUTDIR, "scan.csv"), "w")
println(csv, "C_T,thrust,drag_avg,T-D")
@printf "Re=%.0f, NSTEPS_PER=%d, scanning %d C_T values…\n" Re NSTEPS_PER length(CT_LIST)
println("C_T     thrust    drag    T - D")
for CT in CT_LIST
    T_val, D_val = run_one(CT)
    diff = T_val - D_val
    println(csv, "$CT,$(T_val),$(D_val),$(diff)")
    flush(csv)
    @printf "%.2f  %8.3f  %8.3f  %+8.3f\n" CT T_val D_val diff
end
close(csv)

# Find the self-propulsion C_T (T = D) by linear interpolation
data = open(joinpath(OUTDIR, "scan.csv"), "r") do io
    lines = readlines(io)[2:end]
    [parse.(Float64, split(l, ',')) for l in lines]
end
function find_sign_change(data)
    for i in 1:length(data)-1
        a, b = data[i], data[i+1]
        da, db = a[4], b[4]
        # Strict sign change: use product test (handles zero edge case too)
        if da * db < 0
            t = -da / (db - da)
            return a[1] + t * (b[1] - a[1])
        end
    end
    return nothing
end
self_prop_CT = find_sign_change(data)
println()
if self_prop_CT === nothing
    println("No sign change observed — extend CT_LIST.")
else
    @printf "Estimated self-propulsion C_T ≈ %.4f\n" self_prop_CT
end
println("Scan written to $(joinpath(OUTDIR, "scan.csv"))")
