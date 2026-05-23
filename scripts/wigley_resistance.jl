#!/usr/bin/env julia
#
# 3D Wigley-hull resistance smoke-test.  Fully-submerged hull in a uniform
# stream — measures pressure + viscous drag.  No free surface (Fr = 0).
# Compares to published low-Re submerged-Wigley data (qualitative only —
# this is the first end-to-end ShipShapes + WaterLily integration test).
#
# Setup:
#   - 3D box NX × NY × NZ
#   - Uniform inflow U∞ = 1 along +x, periodic in y/z, exit BC in x
#   - Wigley hull centred at the midplane, oriented along x
#   - Re = U·L/ν (laminar at Re ~ 100; turbulent / transitional otherwise)
#   - Run to steady state, measure F = pressure_force + viscous_force
#   - C_D = F_x / (½ρU²A_wetted)  with A_wetted from ShipShapes

using WaterLily
using ShipShapes
using ShipShapes: wigley_volume, StaticArrays
const SVector = StaticArrays.SVector
using Printf

const NX = parse(Int, get(ENV, "WL_NX", "256"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "64"))
const Re = parse(Float64, get(ENV, "WL_RE", "1000"))
const U∞ = 1.0f0

# Hull dimensions in cell-units.  Pick L_c = NX/2 so hull occupies the
# central half of the domain in x.
const L_c = parse(Float64, get(ENV, "WL_L", string(NX/2)))    # length in cells
const B_c = parse(Float64, get(ENV, "WL_B", string(NY/4)))    # beam in cells
const T_c = parse(Float64, get(ENV, "WL_T", string(NZ/4)))    # draft in cells

# Hull centre (midship) and waterline (z=0 in hull coords)
const xc = NX/2
const yc = NY/2
const zc = NZ/2  # depth at hull waterline (centre of vertical extent)

@printf "=== Wigley hull resistance (fully submerged) ===\n"
@printf "  Grid     = %d × %d × %d\n" NX NY NZ
@printf "  L,B,T    = %.1f, %.1f, %.1f (cells)\n" L_c B_c T_c
@printf "  L/B,L/T  = %.2f, %.2f\n" L_c/B_c L_c/T_c
@printf "  Re_L     = %.0f → ν = %.3e\n" Re U∞ * L_c / Re
@printf "  V_hull   = %.2f cells³ (4·L·B·T/9)\n" wigley_volume(L_c, B_c, T_c)

# Wigley with translation/rotation map so midship is at (xc, yc, zc).
# ShipShapes coordinate convention: x along length, y across beam, z=0
# at waterline (positive up).  Place the hull with z=0 at the box mid-z.
hull_map = (x, t) -> SVector(x[1] - xc, x[2] - yc, x[3] - zc)
hull     = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)

# Wetted surface area (analytic for Wigley; integrating the surface metric).
# For a quick smoke test we use the simple "8/9 · L · (B + 4T)" approximation
# (Series 60 / Wigley empirical), exact enough for an order-of-magnitude C_D.
const A_wet = (8/9) * L_c * (B_c + 4*T_c)
@printf "  A_wetted ≈ %.2f cells² (Wigley empirical)\n" A_wet

ν_cell = U∞ * L_c / Re
sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
    T = Float32,
    ν = Float32(ν_cell),
    body = hull,
    ϵ = 1,
    perdir = (2, 3),
    exitBC = true,
    Δt = 0.25,
)

const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "300"))
@info "Running $(NSTEPS) mom_step!s on Re=$(Re)…"
for step in 1:NSTEPS
    WaterLily.sim_step!(sim; remeasure = false)
    if mod(step, 50) == 0 || step ≤ 3
        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        Ftot = Fp .+ Fv
        u_max = maximum(abs, sim.flow.u)
        CD_press = -2 * Float64(Fp[1]) / (1.0 * U∞^2 * A_wet)  # negate: force on body = -force on fluid; convert via /(½ρU²A)
        CD_visc  = -2 * Float64(Fv[1]) / (1.0 * U∞^2 * A_wet)
        @info @sprintf("step=%4d  Δt=%.2e  |u|=%.3e  F_p=(%.2f,%.2f,%.2f)  F_v=(%.3f,%.3f,%.3f)  CD_p=%.4f  CD_v=%.4f",
            step, sim.flow.Δt[end-1], u_max,
            Fp[1], Fp[2], Fp[3], Fv[1], Fv[2], Fv[3],
            CD_press, CD_visc)
    end
end

Fp = WaterLily.pressure_force(sim)
Fv = WaterLily.viscous_force(sim)
Ftot = Fp .+ Fv

println()
@printf "Final pressure_force (on fluid)  = (%+.3f, %+.3f, %+.3f)\n" Fp[1] Fp[2] Fp[3]
@printf "Final viscous_force  (on fluid)  = (%+.5f, %+.5f, %+.5f)\n" Fv[1] Fv[2] Fv[3]
F_drag_on_body = -Float64(Ftot[1])      # force on body in -x direction = drag
CD_total = 2 * F_drag_on_body / (1.0 * U∞^2 * A_wet)
@printf "Drag on body (x-component)       = %+.4f (cell-units)\n" F_drag_on_body
@printf "C_D total                        = %.4f\n" CD_total
@printf "  (split: C_D_pressure = %.4f, C_D_viscous = %.4f)\n" (
    -2 * Float64(Fp[1]) / (1.0 * U∞^2 * A_wet)) (-2 * Float64(Fv[1]) / (1.0 * U∞^2 * A_wet))

# Order-of-magnitude reference at Re=1000 (laminar regime), a streamlined
# 3D body has C_D ~ 0.01 - 0.05 (similar to a streamlined fuselage).
# This is a SMOKE TEST — not a calibrated validation.
