#!/usr/bin/env julia
#
# Archimedes test with a CYLINDER (smooth SDF) — no corner-gradient
# issues that an axis-aligned box has.  Same 2D box of water, cylinder
# fully submerged at depth.

using WaterLily
using VoF
using Printf

const N      = parse(Int, get(ENV, "WL_N", "128"))
const L_BOX  = 1.0
const ΔX     = L_BOX / N
const G_p    = 9.81
const H_w_p  = 0.5
const ρ_w_p  = 1000.0
const ρ_a_p  = 1.0
const μ_w_p  = 1e-3
const μ_a_p  = 1.8e-5
const U_ref  = sqrt(G_p * H_w_p)
const G_c    = G_p * ΔX / U_ref^2
const H_w_c  = H_w_p / ΔX
const ν_w_c  = (μ_w_p / ρ_w_p) / (U_ref * ΔX)
const ν_a_c  = (μ_a_p / ρ_a_p) / (U_ref * ΔX)
const ρ_w    = ρ_w_p / ρ_a_p; const ρ_a = 1.0
const μ_w_c  = ρ_w * ν_w_c;   const μ_a_c = ρ_a * ν_a_c
const T_NUM  = Float64

# Cylinder (circle in 2D) radius and centre.  Pick depth from WL_CONFIG.
const R_phys  = 0.10
const xc_phys = 0.50
const yc_phys = parse(Float64, get(ENV, "WL_YBOX",
        get(ENV, "WL_CONFIG", "submerged") == "half" ? "0.50" : "0.25"))
const R_c    = R_phys / ΔX
const xc_c, yc_c = xc_phys / ΔX, yc_phys / ΔX

# Smooth SDF: distance to circle.
cyl_sdf(x, t) = sqrt((x[1]-xc_c)^2 + (x[2]-yc_c)^2) - R_c
α₀(_i, x_cell) = (x_cell[2] ≤ H_w_c) ? 1.0 : 0.0

V_box    = π * R_phys^2
# Submerged area for a circle whose centre is at depth h_below = H_w_p - yc_phys.
# h is measured downward from waterline to cylinder centre, +ve if below water.
h_below = H_w_p - yc_phys
V_sub_an = if h_below ≥ R_phys
    V_box                                  # fully submerged
elseif h_below ≤ -R_phys
    0.0                                    # fully out of water
else
    # Standard circular segment formula.  Submerged area:
    #   A = R² · acos(-h/R) + h · √(R² - h²)
    R_phys^2 * acos(-h_below / R_phys) + h_below * sqrt(R_phys^2 - h_below^2)
end
F_buoy_phys = ρ_w_p * G_p * V_sub_an
F_buoy_on_body_cell = ρ_w * G_c * (V_sub_an / ΔX^2)

println("=== Archimedes cylinder test ===")
@printf("  Grid          = %d × %d  (Δx = %.3e m)\n", N, N, ΔX)
@printf("  Free surface  = y = %.3f m (= %.1f cells)\n", H_w_p, H_w_c)
@printf("  Cylinder      = centre (%.3f, %.3f), R = %.3f m (= %.1f cells)\n",
        xc_phys, yc_phys, R_phys, R_c)
@printf("  ρ_w/ρ_a       = %.0f\n", ρ_w/ρ_a)
@printf("  V = πR²       = %.4f m²\n", V_box)
@printf("  F_buoy_phys (on body) = %.2f N/m\n", F_buoy_phys)
@printf("  F_buoy_on_body_cell   = %.2f\n", F_buoy_on_body_cell)
@printf("  expected pressure_force[2] = %.2f (force on fluid)\n", -F_buoy_on_body_cell)

vof = VoFFlow((N, N); α₀=α₀, ρ_w=ρ_w, ρ_a=ρ_a, μ_w=μ_w_c, μ_a=μ_a_c, T=T_NUM)
body = WaterLily.AutoBody(cyl_sdf)
function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ)
end
sim = WaterLily.Simulation((N, N), (T_NUM(0), T_NUM(0)), R_c;
    T=T_NUM, ν=vof.ν,
    g=(i,x,t)-> i==2 ? T_NUM(-G_c) : T_NUM(0),
    Δt=0.5, body=body, ϵ=1, pois_ctor=vof_pois_ctor, U=T_NUM(1.0))

@info "Stepping…"
nsteps = parse(Int, get(ENV, "WL_NSTEPS", "60"))
for step in 1:nsteps
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol=T_NUM(1e-10), pois_itmx=400)
    step_vof!(vof, sim; dt=sim.flow.Δt[end-1])
    if mod(step, 20) == 0 || step ≤ 3
        u_max = maximum(abs, sim.flow.u)
        niter = isempty(sim.pois.n) ? -1 : sim.pois.n[end]
        @info @sprintf("step=%3d  Δt=%.3e  |u|=%.3e  pois=%d",
                       step, sim.flow.Δt[end-1], u_max, niter)
    end
end

F_on_fluid = WaterLily.pressure_force(sim)
F_on_body  = -F_on_fluid
@printf("\n  pressure_force = (%.4f, %.4f)\n", F_on_fluid[1], F_on_fluid[2])
@printf("  Force on body  = (%.4f, %.4f)\n", F_on_body[1], F_on_body[2])
@printf("  Predicted F_buoy_on_body_y = %.4f\n", F_buoy_on_body_cell)
err = (F_on_body[2] - F_buoy_on_body_cell) / F_buoy_on_body_cell
@printf("  Vertical error = %+5.2f %%\n", 100 * err)
exit(abs(err) < 0.10 ? 0 : 1)
