#!/usr/bin/env julia
#
# I2: rotor two-way coupling. So far the rotor's CT was pre-computed
# at the freestream velocity once per simulation; this script
# evaluates `rotor_forces(rotor, U∞, Ω; inflow=trilinear_inflow(flow.u))`
# every step, so the rotor sees the hull-wake-modified flow at its
# disk plane. The wake fraction w = 1 − CT_twoway / CT_freestream is
# the quantity of interest.

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

const NX = 128; const NY = 64; const NZ = 32
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "80"))
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

@printf "=== I2: rotor two-way coupling (wake-fraction trial) ===\n"
@printf "  Grid       = %d × %d × %d   NSTEPS=%d\n" NX NY NZ NSTEPS
@printf "  J=%.2f, Ω=%.3f\n\n" J_op Ω
flush(stdout)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

rotor = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.5,
    chord=(1.0, 0.5), twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)

# Freestream baseline (constant across the run).
r0 = rotor_forces(rotor, U∞, Ω)
CT_free = abs(r0.CT)
@printf "Freestream baseline: |CT|=%.4f thrust=%.3f\n" CT_free abs(r0.CT * 0.5 * U∞^2 * π * R_prop^2)

function vof_pois_ctor(vof)
    (flow) -> begin
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir=(2,))
    end
end

# Run 0: hull only (no rotor) — measure the true Taylor wake at the
# would-be disk plane. This is the freestream-corrected wake fraction
# the towing-tank measurement reports.
@info "Run 0: hull alone (no rotor)"
vof0 = VoFFlow((NX, NY, NZ);
    α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
sim0 = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
    T=Float32, ν=vof0.ν,
    g=(i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor(vof0), U=U∞,
)
for s in 1:NSTEPS
    WaterLily.mom_step!(sim0.flow, sim0.pois; pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vof0, sim0; dt=sim0.flow.Δt[end-1], perdir=(2,))
end
# Sample the wake at the disk location (and 1.5R upstream).
inflow_atdisk_hullonly = trilinear_inflow(sim0.flow.u;
    offset=SVector(prop_xc, prop_yc, prop_zc))
inflow_upstream_hullonly = trilinear_inflow(sim0.flow.u;
    offset=SVector(prop_xc - 1.5f0*R_prop, prop_yc, prop_zc))
v_disk_hullonly = inflow_atdisk_hullonly(SVector(0.0, 0.0, 0.0))
v_up_hullonly   = inflow_upstream_hullonly(SVector(0.0, 0.0, 0.0))
w_taylor = 1f0 - Float32(v_disk_hullonly[1] + U∞) / U∞
@printf "Hull-only wake: V_disk=%.3f  V_upstream(1.5R)=%.3f  → w_Taylor=%+.3f\n" (v_disk_hullonly[1]+U∞) (v_up_hullonly[1]+U∞) w_taylor

# Run 1: free-stream rotor (single precomputed CT, no two-way).
@info "Run A: freestream rotor"
vofA = VoFFlow((NX, NY, NZ);
    α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
simA = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
    T=Float32, ν=vofA.ν,
    g=(i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor(vofA), U=U∞,
)
thrustA = abs(r0.CT * 0.5 * U∞^2 * π * R_prop^2)
function udfA(flow, t; kwargs...)
    smear_force!(flow.f, SVector(Float32(thrustA), 0f0, 0f0),
                 SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
    return nothing
end
dragA = Float64[]
for s in 1:NSTEPS
    WaterLily.mom_step!(simA.flow, simA.pois; udf=udfA, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vofA, simA; dt=simA.flow.Δt[end-1], perdir=(2,))
    Fp = WaterLily.pressure_force(simA); Fv = WaterLily.viscous_force(simA)
    push!(dragA, -(Float64(Fp[1]) + Float64(Fv[1])))
end

# Run 2: two-way rotor (per-step CT from local inflow).
@info "Run B: two-way rotor"
vofB = VoFFlow((NX, NY, NZ);
    α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=Float32)
simB = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
    T=Float32, ν=vofB.ν,
    g=(i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt=0.25f0, body=hull, ϵ=1, perdir=(2,), exitBC=true,
    pois_ctor=vof_pois_ctor(vofB), U=U∞,
)
CT_hist = Float64[]
Vlocal_hist = Float64[]
dragB = Float64[]
# Re-evaluate the VLM only every VLM_EVERY steps — the inflow varies on a
# slow time scale relative to the dt, so updating less often is fine and
# makes the run ~`VLM_EVERY`× faster.
const VLM_EVERY = parse(Int, get(ENV, "WL_VLM_EVERY", "5"))
CT_cached = abs(r0.CT)
thrust_cached = abs(r0.CT * 0.5 * U∞^2 * π * R_prop^2)

# In rotor-frame, the disk plane is at axial = 0. If we sample the
# WaterLily field exactly there, the rotor sees its own induced
# velocity (≈ a·U∞ at the actuator-disk plane) and we double-count
# its own race — the resulting CT is *higher* than freestream, not
# lower. To match the towing-tank "Taylor wake" convention we sample
# 1.5·R upstream of the disk plane, where the rotor's induced
# velocity is small but the hull wake is essentially fully developed.
const PROBE_UPSTREAM = 1.5f0 * R_prop
function inflow_at_disk(flow)
    return trilinear_inflow(flow.u;
        offset=SVector(prop_xc - PROBE_UPSTREAM, prop_yc, prop_zc))
end

function udfB(flow, t; kwargs...)
    smear_force!(flow.f, SVector(Float32(thrust_cached), 0f0, 0f0),
                 SVector(prop_xc, prop_yc, prop_zc); ε=2.5f0)
    return nothing
end
for s in 1:NSTEPS
    if s == 1 || s % VLM_EVERY == 0
        inflow_fn = inflow_at_disk(simB.flow)
        r = rotor_forces(rotor, U∞, Ω; inflow=inflow_fn)
        global CT_cached = abs(Float64(r.CT))
        global thrust_cached = abs(r.CT * 0.5 * U∞^2 * π * R_prop^2)
        v = inflow_fn(SVector(0.0, 0.0, 0.0))
        push!(Vlocal_hist, Float64(v[1] + U∞))
    end
    push!(CT_hist, CT_cached)
    WaterLily.mom_step!(simB.flow, simB.pois; udf=udfB, pois_tol=1f-6, pois_itmx=50)
    step_vof_mules!(vofB, simB; dt=simB.flow.Δt[end-1], perdir=(2,))
    Fp = WaterLily.pressure_force(simB); Fv = WaterLily.viscous_force(simB)
    push!(dragB, -(Float64(Fp[1]) + Float64(Fv[1])))
    s % 10 == 0 && println("  step=$s CT=$CT_cached")
end

# Reduce: average over last 25 %.
CT_A = CT_free
CT_B = mean(CT_hist[end-NSTEPS÷4:end])
Vloc_B = mean(Vlocal_hist[max(1, end-length(Vlocal_hist)÷4):end])
drag_A = mean(dragA[end-NSTEPS÷4:end])
drag_B = mean(dragB[end-NSTEPS÷4:end])

wake = 1 - CT_B / CT_A

println()
@printf "═══════════════════════════════════════════════════════════════════════════\n"
@printf "                          freestream   two-way\n"
@printf "  CT_rotor               %8.4f    %8.4f\n" CT_A CT_B
@printf "  V_local at disk           %.3f       %.3f\n" U∞ Vloc_B
@printf "  hull drag               %8.3f    %8.3f\n" drag_A drag_B
@printf "  wake fraction w = 1 − CT_B/CT_A  =  %+.3f  (%+.1f%%)\n" wake 100*wake
@printf "═══════════════════════════════════════════════════════════════════════════\n"

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "rotor_twoway"))
mkpath(OUTDIR)
open(joinpath(OUTDIR, "summary.csv"), "w") do io
    println(io, "case,CT,Vlocal,drag")
    println(io, "freestream,$CT_A,$U∞,$drag_A")
    println(io, "twoway,$CT_B,$Vloc_B,$drag_B")
    println(io, "wake_fraction,$wake,,")
end
println("Summary written to $(joinpath(OUTDIR, "summary.csv"))")
