# RESULTS — gentle sloshing (MULES + interface compression validation)

Follow-up 2″ of [RESULTS-damBreak.md](./RESULTS-damBreak.md): the
damBreak sweep showed compression keeps the interface sharp but that a
sharp algebraic interface at ρ-ratio 1000 destabilizes the u-advecting
coupling under *wall-slam*. This case probes compression in the regime
it exists for — a smooth standing wave that never breaks — where the
payoff is sharpness retention over many periods.

**Verdict: VALIDATED.** In the gentle regime MULES + c_α=1 is stable
over 5.5 periods at ρ-ratio 1000, reproduces the linear sloshing period
to 0.5 %, conserves mass to 1.3e-3 %, and holds the interface ~2.4×
sharper than either alternative — with no sign of the damBreak
instability.

## Setup

Closed 2D tank (no-slip walls), L = 0.585 m square, water depth
h = L/2, free surface tilted by the first standing mode
a·cos(πx/L) with a = 0.05·L. ρ_w/ρ_a = 1000, physical μ for both
phases, no surface tension. N = 128², 5 s ≈ 5.5 periods of the
analytic mode-1 (ω² = g·k·tanh(kh), k = π/L ⇒ T = 0.904 s).
Driver: `scripts/sloshing_mules.jl` (VoF.jl + WaterLily
`foam-integration`, `density_coefficient!` Poisson path, struct
tol = 1e-8 / itmx = 200).

## Results (N=128, 5 s, 547 steps each)

| config | survives | period vs analytic | mass drift | interface width (t≈1 s → 5 s) |
|---|---|---|---|---|
| MULES c_α=1 | ✓ | 0.909 s (+0.5 %) | 0.0013 % | **165 → 170** (≈1.3 cells/column, steady) |
| MULES c_α=0 | ✓ | 0.910 s (+0.7 %) | 0.0002 % | 403 → 405 (≈3.2 cells/column) |
| `step_vof!` clamp+repair | ✓ | 0.910 s (+0.7 %) | 0.0003 % | 403 → 405 |

Width = count(0.05 < α < 0.95) over the interior; the surface spans
128 columns, so sharp ≈ 1–3 cells per column. Elevation is the
mass-equivalent column sum, so the period estimate (zero-crossings of
η_left) is robust to smearing.

Notes:

- All three configs nail the linear period — the variable-density
  projection (`L = μ₀/ρ_face`) carries the gravity-wave dynamics
  correctly regardless of advection scheme.
- Without compression the width *plateaus* (~3 cells/col) rather than
  homogenizing as in damBreak: gentle straining + vanLeer's own
  steepening find an equilibrium width. Compression halves-plus it
  and keeps it flat — that margin is what matters for Kelvin-pattern
  fidelity where the wake strain field keeps trying to smear the
  surface.
- |u| stays ≈0.08 cell-units throughout — nowhere near the wall-slam
  regime (|u| ~ 50+) that killed the sharp damBreak runs.

## Recommendation (final, supersedes interim damBreak note)

| flow regime | scheme |
|---|---|
| violent (slamming, breaking, dam-break) | `step_vof!` + `mass_repair` (clamp keeps it sharp *and* stable) |
| gentle free surface (ship wave-making, seakeeping below breaking) | `step_vof_mules!` + `c_α = 1` |

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
WL_N=128 WL_TEND=5.0 WL_MULES=1 WL_CALPHA=1 WL_TAG=N128_calpha1 \
    julia +1.12 --project=. scripts/sloshing_mules.jl
WL_N=128 WL_TEND=5.0 WL_MULES=1 WL_CALPHA=0 WL_TAG=N128_calpha0 \
    julia +1.12 --project=. scripts/sloshing_mules.jl
WL_N=128 WL_TEND=5.0 WL_MULES=0 WL_TAG=N128_clamprepair \
    julia +1.12 --project=. scripts/sloshing_mules.jl
```

CSVs: `runs/sloshing_mules/slosh_N128_*.csv` (t, η_left, η_right,
m/m₀, width).
