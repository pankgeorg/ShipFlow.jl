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

# Build the tabulated SDF.  Sample box ~ 1.4× hull extent in each axis.
const N_S = (96, 48, 48)   # sample points
const half_L = L_c * 0.7
const half_B = B_c * 0.7
const T_pad  = T_c * 0.2
spx = (2 * half_L) / (N_S[1] - 1)
spy = (2 * half_B) / (N_S[2] - 1)
spz = (T_c + 2 * T_pad) / (N_S[3] - 1)
origin = (-half_L, -half_B, -T_c - T_pad)

# We sample in HULL-FRAME coords (midship at 0, waterline at z=0).
table = sample_sdf((x, t) -> wigley_sdf(x, L_c, B_c, T_c),
                   origin, (spx, spy, spz), N_S; T = Float32)
@printf "  Tabulated SDF: %d × %d × %d samples, spacing (%.3f, %.3f, %.3f)\n" N_S... spx spy spz

# Hull map: world (xc, yc, zc) is hull-frame origin.  Sample box origin
# is relative to hull frame; the table's own coordinates are world-aligned
# after we feed `x_hull = x_world - centre`.
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
