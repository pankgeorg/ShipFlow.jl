#!/usr/bin/env julia
#
# Same as wigley_resistance.jl but uses a TabulatedHull built by
# sampling the analytic Wigley SDF.  Validates the tabulated-SDF path
# that the DTC hull (and any other non-analytic ship-shape) will use.

using WaterLily
using ShipShapes
using ShipShapes: wigley_sdf, sample_sdf, tabulated_sdf, StaticArrays
const SVector = StaticArrays.SVector
using Printf

const NX = parse(Int, get(ENV, "WL_NX", "192"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const Re = parse(Float64, get(ENV, "WL_RE", "1000"))
const U∞ = 1.0f0

const L_c = parse(Float64, get(ENV, "WL_L", string(NX/2)))
const B_c = parse(Float64, get(ENV, "WL_B", string(NY/4)))
const T_c = parse(Float64, get(ENV, "WL_T", string(NZ/4)))
const xc = NX/2;  const yc = NY/2;  const zc = NZ/2

@printf "=== Wigley resistance via TabulatedHull (validates DTC SDF path) ===\n"
@printf "  Grid     = %d × %d × %d\n" NX NY NZ
@printf "  L,B,T    = %.1f, %.1f, %.1f (cells)\n" L_c B_c T_c
@printf "  Re_L     = %.0f → ν = %.3e\n" Re U∞ * L_c / Re

# Tabulate the SDF over the FULL simulation domain (in hull-frame).
# AutoBody needs a smooth SDF everywhere — sampling only the hull
# region leaves gradient=0 outside which gives NaN normals.
# Use 1-cell spacing matching the WaterLily grid for best accuracy.
const N_S = (NX + 2, NY + 2, NZ + 2)
# Map: x_world ↦ x_hull = x_world - centre, so the table covers
# x_hull ∈ [-xc, NX - xc] × [-yc, NY - yc] × [-zc, NZ - zc] in 1-cell steps.
origin = (-Float64(xc), -Float64(yc), -Float64(zc))
spacing = (1.0, 1.0, 1.0)

table = sample_sdf((x, t) -> wigley_sdf(x, L_c, B_c, T_c),
                   origin, spacing, N_S; T = Float32)
@printf "  Tabulated SDF: %d × %d × %d samples over full domain\n" N_S...

hull_map = (x, t) -> SVector(x[1] - xc, x[2] - yc, x[3] - zc)
hull     = tabulated_sdf(table; map = hull_map)

A_wet = (8/9) * L_c * (B_c + 4*T_c)
@printf "  A_wetted ≈ %.2f cells²\n\n" A_wet

ν_cell = U∞ * L_c / Re
sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0f0, 0f0), L_c;
    T = Float32, ν = Float32(ν_cell), body = hull, ϵ = 1,
    perdir = (2, 3), exitBC = true, Δt = 0.25,
)

const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "200"))
@info "Running $(NSTEPS) mom_step!s on Re=$Re (tabulated SDF)…"
for step in 1:NSTEPS
    WaterLily.sim_step!(sim; remeasure = false)
    if mod(step, 50) == 0 || step ≤ 3
        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        CD = -2 * (Float64(Fp[1]) + Float64(Fv[1])) / (1.0 * U∞^2 * A_wet)
        @info @sprintf("step=%4d  CD=%.4f  (CD_p=%.4f, CD_v=%.4f)",
            step, CD, -2*Float64(Fp[1])/A_wet, -2*Float64(Fv[1])/A_wet)
    end
end

Fp = WaterLily.pressure_force(sim)
Fv = WaterLily.viscous_force(sim)
CD = -2 * (Float64(Fp[1]) + Float64(Fv[1])) / (1.0 * U∞^2 * A_wet)
println()
@printf "Final C_D (tabulated) = %.4f\n" CD
@printf "Reference C_D (analytic, from RESULTS-wigley.md) = 0.0233\n"
δ = 100 * (CD - 0.0233) / 0.0233
@printf "Deviation from analytic-SDF run = %+5.2f %%\n" δ
println("\n(small deviations from interpolation error are expected.)")
