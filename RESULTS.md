# RESULTS — overview

The five-package Julia ship-CFD stack
(WaterLily + Turbulence + VoF + ShipShapes + Propellers) running on
`pankgeorg/*` produces the following validated and predicted
quantities. Each row links to a dedicated RESULTS doc with the full
methodology, data, and caveats.

## Validations (against external references)

| Case                                            | Result          | Reference / gate         | Doc                          |
|-------------------------------------------------|----------------:|--------------------------|------------------------------|
| Cylinder Re=100 — Cd, Cl_pp, St                 | within 5-10%   | Williamson 1996          | [RESULTS-cylinder.md]        |
| Channel Re_τ=395 — bulk u_x(y) RMS              | 0.028 (pre-fix) | OF Smagorinsky           | [RESULTS-channel.md]         |
| damBreak ρ=10:1 — front-position RMS            | **4.1 %**       | Martin-Moyce 1952, gate ±10% | [RESULTS-damBreak.md]    |
| Submerged-cylinder Archimedes — F_buoy          | within 2.4 %    | analytic ρ·g·V, gate ±5% | [RESULTS-archimedes.md]      |
| Actuator disk — U_disk, U_wake                  | -2.9 %, -4.4 %  | 1D Froude-Rankine        | [RESULTS-propeller.md]       |
| Submerged Wigley Re=1000 — C_D_viscous          | 0.0213          | Blasius 0.0210 (1 % match) | [RESULTS-wigley.md]        |
| Free-surface Wigley wave-resistance peak        | Fr ≈ 0.35       | Wigley literature 0.3–0.4 | [RESULTS-Fr-scan.md]        |

## Quantitative predictions from the stack

| Quantity                                                       | Value          | Doc                       |
|----------------------------------------------------------------|---------------:|---------------------------|
| Self-propulsion C_T (Wigley, half-submerged, Fr=0.25, high-Re) | **2.24 ± 0.03** (Re 10⁴–10⁵) | [RESULTS-selfprop.md] |
| Thrust deduction factor t at Fr=0.25                           | **0.36**       | [RESULTS-selfprop.md]     |
| t at Fr=0.20                                                   | 0.42           | "                         |
| t at Fr=0.35                                                   | 0.22           | "                         |
| Drag at self-propulsion vs bare hull                           | +55 %          | "                         |

## Integration milestones

- **Five-package headline simulation** runs end-to-end on a 96×48×48 grid in <1 min: see [RESULTS-headline.md]. Visualisations: `runs/wigley_snapshot_Fr025/*.png`, `runs/wigley_snapshot_Fr035/*.png`.
- **End-to-end integration test** in `test/runtests.jl` exercises all five packages with 9 finiteness + mass-conservation assertions in 27 s — CI-suitable.

## Headline plots

- `runs/wigley_selfprop_scan/self_propulsion.png` — T(C_T) crosses D(C_T) at C_T=2.33
- `runs/wigley_selfprop_scan/CT_vs_Re.png` — Re-scan asymptote at C_T≈2.24
- `runs/wigley_selfprop_FrScan/CT_t_vs_Fr.png` — Fr-scan; t falls with Fr
- `runs/wigley_resistance_Fr/drag_vs_Fr.png` — Wigley wave-resistance hump at Fr=0.35

## Major corrections found and fixed during validation

| Bug                                                                            | Resolution           |
|--------------------------------------------------------------------------------|----------------------|
| WaterLily `Flow` copied user's ν array — LES updates never reached `conv_diff!` | commit `e4b8854` — store reference, regression test added |
| WaterLily `viscous_force` not consistent with array-ν Hook 1                   | commit `dee8816` — `_ν(ν,I)` lookup |
| `Hook 3` Poisson tolerance wasn't routable through `mom_step!` for variable-ρ projection | commit `318f76d` |
| Channel395 "validation" was actually quasi-DNS, not LES                         | flagged in [RESULTS-channel.md]; re-calibration required |

## What's left for Phase-3 release

1. **DTC hull SDF** — `TabulatedHull` infra ready (6.2 % round-trip error on Wigley sample). Need DTC offsets data.
2. **Re=10⁶ self-propulsion** — would convert the Re-asymptote prediction (C_T≈2.24) into a calibratable result for ship-scale Cf.
3. **Channel395 LES re-calibration** — re-tune driver `g_x` now that LES is actually live.
4. **MULES α-redistribution** — drop damBreak mass loss from 5 % to the stated 0.1 % Phase-2 gate.
5. **OpenFOAM cross-validation** — blocked on arm64 (no qemu-user-static).

See [HANDOFF.md](../HANDOFF.md) at the foam root for the full resume snapshot and the per-package roadmaps.

[RESULTS-cylinder.md]:  ./RESULTS-cylinder.md
[RESULTS-channel.md]:   ./RESULTS-channel.md
[RESULTS-damBreak.md]:  ./RESULTS-damBreak.md
[RESULTS-archimedes.md]: ./RESULTS-archimedes.md
[RESULTS-propeller.md]: ./RESULTS-propeller.md
[RESULTS-wigley.md]:    ./RESULTS-wigley.md
[RESULTS-Fr-scan.md]:   ./RESULTS-Fr-scan.md
[RESULTS-selfprop.md]:  ./RESULTS-selfprop.md
[RESULTS-headline.md]:  ./RESULTS-headline.md
