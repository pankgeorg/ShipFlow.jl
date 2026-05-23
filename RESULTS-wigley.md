# RESULTS — Wigley hull resistance smoke-test

Sixth cross-validation: 3D Wigley hull (ShipShapes.jl) in uniform
stream, fully submerged (no free surface). Measures pressure + viscous
drag and compares to Blasius flat-plate friction theory.

This is the first integration test of WaterLily + ShipShapes
end-to-end. **No VoF, no propeller, no free surface** — just a hull
in laminar flow. It exists to flag wiring bugs before the full
DTC self-propulsion calculation in Phase 3.

## Headline

Re_L = 1000, hull 96 × 12 × 8 cells, domain 192 × 64 × 48.
Steady state reached by step ~50.

| Quantity              | Measured  | Reference                       | Notes |
|----------------------:|----------:|---------------------------------:|-------|
| C_D total             | 0.0233    | —                                | -     |
| C_D viscous           | 0.0213    | Blasius flat plate: 0.664/√Re = 0.021 | matches to 1% |
| C_D pressure (form)   | 0.0020    | —                                | ~10x smaller than friction, as expected for slender hull |

**Excellent**: the viscous drag matches Blasius flat-plate theory
(0.0213 vs 0.0210). Pressure drag is a small fraction, consistent with
a streamlined slender body.

## What was actually run

- Domain 192 × 64 × 48 cells, periodic in y/z, convective exit in x
- Wigley hull L = 96, B = 12, T = 8 cells (aspect L/B = 8, L/T = 12)
- Translated so midship sits at (NX/2, NY/2, NZ/2)
- Uniform inflow U∞ = 1, ν = U·L/Re = 0.096 (Re_L = 1000)
- BDIM kernel ϵ = 1
- 200 momentum steps (steady at ~50)

Driver: [`scripts/wigley_resistance.jl`](./scripts/wigley_resistance.jl).

## Caveats

- Wetted-area estimate uses the (8/9)L(B+4T) empirical formula — exact
  to ~5 % for the Wigley parabolic hull. A volumetric integration would
  be tighter.
- F_z = +10.6 is non-zero (≈ 25 % of F_x). The hull is symmetric in z
  about the midplane (waterline) but the grid is at NZ/2 = 24 — likely
  a half-cell asymmetry. Investigate before turning into a real test.
- Re = 1000 is laminar; the published Wigley resistance data is in the
  turbulent regime (Re > 10⁵). Layer-2 OF cross-validation will use
  turbulent Re with Turbulence.jl.

## Outstanding follow-ups

| # | Item                                                                |
|---|----------------------------------------------------------------------|
| 1 | F_z asymmetry — instrument with a finer waterline alignment          |
| 2 | Higher Re (10⁴, 10⁵) with Turbulence.jl Smagorinsky                  |
| 3 | Free-surface Wigley with VoF (Fr > 0) — wave-resistance calculation |
| 4 | Compare to published Wigley submerged-hull data (Bai & Webster 2003) |
| 5 | DTC hull resistance (Phase 3 release-blocker, requires tabulated SDF) |

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
JULIA_NUM_THREADS=auto WL_NX=192 WL_NY=64 WL_NZ=48 WL_L=96 WL_B=12 WL_T=8 \
    WL_RE=1000 WL_NSTEPS=200 \
    julia +1.12 --project=. scripts/wigley_resistance.jl
```
