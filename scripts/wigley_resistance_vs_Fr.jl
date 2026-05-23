#!/usr/bin/env julia
#
# Bare-hull resistance vs Froude number — classical naval-arch result.
# Repeat the half-submerged Wigley simulation across several Fr; record
# steady-state drag; output C_W (wave-resistance coefficient) vs Fr.

using WaterLily
using VoF
using ShipShapes
using ShipShapes: StaticArrays
using Turbulence
const SVector = StaticArrays.SVector
using Printf

const NX = 96; const NY = 48; const NZ = 48
const L_c = 48; const B_c = 10; const T_c = 6
const ρ_w = 10.0; const ρ_a = 1.0
const U∞  = 1.0
const Re  = parse(Float64, get(ENV, "WL_RE", "5000"))
const ν_w_c = U∞ * L_c / Re
const ν_a_c = ν_w_c * 18
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * ν_a_c
const T_NUM = Float64

const H_w_c   = NZ/2
const hull_xc = NX/3; const hull_yc = NY/2; const hull_zc = H_w_c

const FR_LIST = let s = get(ENV, "WL_FR_LIST", "0.15,0.20,0.25,0.30,0.35,0.40")
    [parse(Float64, x) for x in split(s, ",")]
end
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "80"))
const A_wet = (8/9) * L_c * (B_c + 4*T_c) / 2

const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "wigley_resistance_Fr"))
mkpath(OUTDIR)

α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0
hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)

function run_one(Fr)
    G_c = U∞^2 / (Fr^2 * L_c)
    vof = VoFFlow((NX, NY, NZ);
        α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = T_NUM)
    turb = WALE((NX, NY, NZ); Cw = T_NUM(0.5), ν₀ = T_NUM(0))
    hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
    function vof_pois_ctor(flow)
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ),
        (T_NUM(U∞), T_NUM(0), T_NUM(0)), Float64(L_c);
        T = T_NUM, ν = turb.ν,
        g = (i, x, t) -> i == 3 ? T_NUM(-G_c) : T_NUM(0),
        Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
        pois_ctor = vof_pois_ctor, U = T_NUM(U∞),
    )
    drag_hist = Float64[]
    for step in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois;
            pois_tol = T_NUM(1e-8), pois_itmx = 100)
        step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
        update_νt!(turb, sim.flow.u, vof.ν)
        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        push!(drag_hist, -(Float64(Fp[1]) + Float64(Fv[1])))
    end
    N_avg = max(1, NSTEPS ÷ 4)
    drag = sum(drag_hist[end-N_avg+1:end]) / N_avg
    Fp_final = WaterLily.pressure_force(sim)
    Fv_final = WaterLily.viscous_force(sim)
    drag_press = -Float64(Fp_final[1])
    drag_visc  = -Float64(Fv_final[1])
    return drag, drag_press, drag_visc
end

@printf "Re=%.0f, NSTEPS=%d, scanning Fr=%s …\n" Re NSTEPS join(string.(FR_LIST), ",")
println("Fr     drag    drag_p   drag_v   C_T    C_P    C_V")

csv = open(joinpath(OUTDIR, "drag_vs_Fr.csv"), "w")
println(csv, "Fr,drag,drag_press,drag_visc,C_T,C_P,C_V")
for Fr in FR_LIST
    d, dp, dv = run_one(Fr)
    C_T = 2 * d  / (1.0 * U∞^2 * A_wet)
    C_P = 2 * dp / (1.0 * U∞^2 * A_wet)
    C_V = 2 * dv / (1.0 * U∞^2 * A_wet)
    println(csv, "$Fr,$d,$dp,$dv,$C_T,$C_P,$C_V")
    flush(csv)
    @printf "%.2f  %6.2f  %6.2f   %6.2f   %.4f  %.4f  %.4f\n" Fr d dp dv C_T C_P C_V
end
close(csv)
println("\nWrote $(joinpath(OUTDIR, "drag_vs_Fr.csv"))")
