#!/usr/bin/env julia
#
# Hydrostatic-balance smoke test for VoFFlow.
#
# Setup: water column filling left half of the box, air above. With one or
# more `mom_step!` calls under gravity, water should remain stationary and
# pressure should match ρ·g·H.  Velocity at any point in water should stay
# below ~1% of √(g·H).
#
# Sweeps ρ ratio ∈ {1, 10, 100, 1000} and prints velocity / pressure stats.

using WaterLily
using VoF
using Printf
using Statistics

const N        = parse(Int, get(ENV, "WL_N", "64"))
const L_BOX    = 1.0
const ΔX       = L_BOX / N
const G_p      = 9.81
const H_w_p    = 0.5            # water column height (half the box)
const COL_W_p  = 0.5            # water occupies the left half too — actually
                                # for the static test we'll fill the WHOLE
                                # bottom half with water (= flat free surface)
const U_ref    = sqrt(G_p * H_w_p)
const G_c      = G_p * ΔX / U_ref^2
const COL_W_c  = COL_W_p / ΔX
const H_w_c    = H_w_p / ΔX
const ν_w_p    = 1.0e-6
const ν_a_p    = 1.8e-5
const ν_w_c    = ν_w_p / (U_ref * ΔX)
const ν_a_c    = ν_a_p / (U_ref * ΔX)
const T_NUM    = Float64

# Flat free surface initial condition: water fills y ≤ H_w_c.
α₀(_i, x_cell) = (x_cell[2] ≤ H_w_c) ? 1.0 : 0.0

function run_one(rho_ratio; nsteps = 30, tol = 1e-10, itmx = 400)
    ρ_w = T_NUM(rho_ratio)
    ρ_a = T_NUM(1.0)
    μ_w_c = ρ_w * ν_w_c
    μ_a_c = ρ_a * ν_a_c
    vof = VoFFlow((N, N);
        α₀ = α₀,
        ρ_w = ρ_w, ρ_a = ρ_a,
        μ_w = μ_w_c, μ_a = μ_a_c,
        T = T_NUM,
    )
    flow = WaterLily.Flow((N, N), (T_NUM(0), T_NUM(0));
        T = T_NUM,
        ν = vof.ν,
        g = (i, x, t) -> i == 2 ? T_NUM(-G_c) : T_NUM(0),
        Δt = 0.5,
    )
    pois = WaterLily.MultiLevelPoisson(flow.p, copy(vof.L), flow.σ)
    sim = (flow = flow, pois = pois)

    for _ in 1:nsteps
        WaterLily.mom_step!(sim.flow, sim.pois;
            pois_tol = T_NUM(tol), pois_itmx = itmx)
        step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
    end

    # Bring p back to a meaningful absolute pressure by anchoring to air-top.
    # Air column max y is N (cells); take p there as reference.
    p = sim.flow.p
    p_ref = p[2, N+1]
    p_anchored = p .- p_ref

    # Diagnostics in interior cells only
    in_y = 2:N+1
    u_x = view(sim.flow.u, 2:N, in_y, 1)
    u_y = view(sim.flow.u, 2:N+1, 2:N, 2)
    α = view(vof.α, 2:N+1, 2:N+1)

    # Hydrostatic prediction in cell-units of pressure (ρ_w · g_cell):
    # p_hydro(y_cell) = ρ_w · g_cell · (H_w_c - y_cell)  for y_cell < H_w_c
    ymid_water = max(2, Int(round(H_w_c / 2)) + 1)
    y_water_top_cell = max(2, Int(round(H_w_c)))
    j_bot = 2                       # bottom interior row index in p (1-indexed)
    j_water_top = j_bot + y_water_top_cell - 1
    p_bottom = p_anchored[Int(round(N/2)), j_bot]
    p_water_top = p_anchored[Int(round(N/2)), j_water_top]
    # Hydrostatic: p_bottom relative to top of domain = weight of column
    # above per unit area = g·(ρ_a·H_air + ρ_w·H_water).
    H_air_c = N - H_w_c
    p_hydro_bottom = G_c * (ρ_a * H_air_c + ρ_w * H_w_c)
    err_pct = 100 * (p_bottom - p_hydro_bottom) / p_hydro_bottom

    @printf "%s ρ_w/ρ_a=%6.0f | u_x∈[%+.2e,%+.2e]  u_y∈[%+.2e,%+.2e]  |u|_max=%.2e  p_bot=%.3f  predict=%.3f  err=%+5.1f%%  niter=%d\n" (
        abs(err_pct)<5 ? "✓" : "⚠"
    ) rho_ratio minimum(u_x) maximum(u_x) minimum(u_y) maximum(u_y) maximum(abs,
        sim.flow.u) p_bottom p_hydro_bottom err_pct sim.pois.n[end]
end

println("=== Hydrostatic-balance check (N=$N, $(Int(H_w_c)) cells of water) ===")
println("  U_ref=$(round(U_ref,digits=4)) m/s   g_cell=$(round(G_c,digits=4))")
for rho in (1.0, 10.0, 100.0, 1000.0)
    run_one(rho)
end
