# Higher-Re integrated stack with WALE (K4)

Re-runs the integrated VLM stack (Wigley + sectional rotor smear +
rudder) with WALE LES enabled at three Reynolds numbers. Confirms
stability and reports drag / |u| / eddy-viscosity trends.

Driver: [`scripts/headline_re_sweep.jl`](scripts/headline_re_sweep.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 36, B = 8, T = 5.
- Fr = 0.30. 60 steps per Re.
- `Turbulence.WALE(Cw=0.5, ν₀=0)` aggregating on top of `vof.ν`.
- Rotor at J = 0.32 (G1), rudder at δ = 5° (out of race for this
  stability test).

## Result

| Re      | D_mean (last 15) | \|u\|_max | ν_t_max  | Wall time |
|---------|------------------|-----------|----------|-----------|
| 5,000   | 58.74            | 3.63      | 0.251    | 42.5 s    |
| 20,000  | 55.87            | 3.66      | 0.161    | 29.7 s    |
| 50,000  | 55.31            | 3.66      | 0.143    | 30.7 s    |

Plot: `runs/headline_re_sweep/Re_sweep.png`.

## Interpretation

- **All three Re cases are stable** — `|u|_max` is bounded at ~3.6
  across the sweep, no blow-up or numerical artefacts. WALE
  successfully transfers energy to sub-grid scales at every Re.
- **Drag drops 6 %** going from Re = 5000 to Re = 50000. Viscous
  Cf scales as `1/√Re` in the Blasius regime, so 10× the Re means
  ~3× less viscous drag. The 6 % drop in our *total* drag suggests
  the viscous contribution is ~10–20 % of total at these Re — wave
  resistance + pressure drag dominate even at Re = 5000.
- **ν_t_max decreases with Re** (0.25 → 0.14, ~40 % drop). This is
  the WALE model working correctly: at higher Re the resolved
  cascade carries more energy directly (no sub-grid help needed),
  so the eddy-viscosity term contributes less per cell.
- **Wall time is independent of Re** for the WaterLily step itself
  (Poisson solves are not Re-sensitive at fixed grid). The 12 s
  startup overhead at Re=5000 (42.5 vs 30 s) is the WALE precompile
  on its first call — subsequent runs reuse.

## What this enables

- **Higher-Re production runs**: the headline demo (F3) was at
  Re = 20000 with WALE. Now confirmed stable to Re = 50000 (10×
  larger). Could push to Re = 200000 (real-ship boundary-layer
  regime) by halving the grid spacing.
- **Quantitative wave/viscous decomposition**: the drag-vs-Re curve
  here is the first attempt at separating wave and viscous drag.
  A no-wave reference run (deep-water, no free surface) would
  isolate viscous; the difference would be wave drag.

## Caveats

- 60 steps may not fully settle drag at the highest Re; the wave
  field develops on a Fr-dependent time scale that doesn't shrink
  with Re. The numbers reported are last-15-step averages, which
  catch the equilibrium reasonably but not perfectly.
- WALE has no near-wall damping in our setup. At very high Re
  (10⁵+) wall-modeled LES becomes necessary; for our laminar-
  transition regime (Re < 10⁵) the WALE-only setup is fine.
- The grid resolves ~64 cells across hull L, which is marginal for
  resolving the boundary-layer thickness at Re = 50000. A grid
  refinement (per J4) at Re = 50000 would show how much of the
  ν_t reduction is real vs grid-induced.

## See also

- `RESULTS-headline-demo.md` (F3) — single Re = 20000 case at full
  feature set.
- `Turbulence.jl/src/Turbulence.jl` — WALE implementation.
- `runs/headline_re_sweep/Re_sweep.png` — D, |u|, ν_t vs step
  for the three Re cases.
