# RESULTS — bare-hull resistance vs Froude number

Classical naval-architecture result: the wave-resistance hump of a
half-submerged Wigley hull, swept across Fr = 0.15 to 0.40 at Re=5000.

## Setup

Same hull + free-surface stack as the self-propulsion test, no
propeller:

- Wigley hull L × B × T = 48 × 10 × 6 cells, half-submerged (waterline = NZ/2)
- 96 × 48 × 48 grid, ρ_w/ρ_a = 10:1, Re_L = 5000
- WALE turbulence model on top of VoF per-cell ν
- 80 mom_step!s per Fr value; drag averaged over the last 25%
- Driver: [`scripts/wigley_resistance_vs_Fr.jl`](./scripts/wigley_resistance_vs_Fr.jl)

## Results

(Re-run with the corrected WaterLily Hook 1 array-ν reference
semantics; see "Note on the LES coupling fix" below.)

| Fr   | total drag | C_T (= 2D/(ρU²·A_wet)) |
|-----:|-----------:|-----------------------:|
| 0.15 |  27.5      | 0.076                  |
| 0.20 |  43.0      | 0.118                  |
| 0.25 |  47.7      | 0.132                  |
| 0.30 |  50.2      | 0.139                  |
| 0.35 |  50.6      | **0.140** *(peak)*     |
| 0.40 |  48.8      | 0.135                  |
| C_V  | ~3.5       | ~0.010                 |

C_T rises monotonically until **Fr ≈ 0.35**, then drops — the
classical wave-resistance hump driven by the Kelvin wave system. For
slender ship hulls (Wigley, Series 60, KCS) this hump typically sits
between Fr = 0.3 and Fr = 0.4 in published experimental and CFD data
(Bai & Webster 2003, Larsson & Stern 1999). Our peak at 0.35 lands in
the middle of that range.

Viscous coefficient C_V ≈ 0.010 holds across Fr (Cf depends on Re, not
Fr), close to laminar Blasius for the half-wetted area
(Cf_Blasius @ Re=5000 = 0.664/√Re ≈ 0.0094).

## Note on the LES coupling fix

The original scan reported a flat-low-Fr curve (Fr=0.15 → C_T=0.034)
and identically zero viscous component. Both were artefacts of a
WaterLily bug where `Flow` stored a COPY of the per-cell ν array
instead of the user's reference, so the turbulence model's updates
never reached the momentum equation. Fixed in WaterLily commit
`e4b8854`. With the fix:

* Low-Fr drag rose substantially (0.15: 12 → 28, +130 %): the LES
  now contributes proper SGS dissipation in the boundary layer.
* High-Fr drag fell (0.35: 59 → 51): wave-making is more efficient
  with a properly resolved boundary layer.
* C_V is finite (~0.010) for the first time.

The hump position (Fr ≈ 0.35) is robust across the bug fix.

## Caveats

- 80 steps is not enough to fully settle the wake at low Fr (Fr=0.15
  result is most affected; mean drag=12.1 but final-step pressure was
  −5.2 — the system was still in transient).
- Viscous component shows as ≈ 0 — likely the viscous_force kernel
  needs revisiting when fed a per-cell ν from VoFFlow (logged as a
  follow-up).
- 96×48×48 is too coarse to resolve the bow wave properly; the hump
  position is qualitatively right but the magnitude is unvalidated.

## What this validates

End-to-end, the full stack now produces a physically-interpretable
curve, not just a single number. Combined with the [self-propulsion
result](./RESULTS-selfprop.md), this closes the WiP "everything
works together" question — what's left is calibration and scale-up.

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
JULIA_NUM_THREADS=auto WL_RE=5000 WL_NSTEPS=80 \
    julia +1.12 --project=. scripts/wigley_resistance_vs_Fr.jl
```
