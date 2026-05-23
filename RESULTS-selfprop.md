# RESULTS — Self-propulsion via C_T parameter scan

First quantitative self-propulsion result from the Julia WaterLily
stack: a 5-package simulation (WaterLily + Turbulence(WALE) + VoF +
ShipShapes(Wigley) + Propellers(ActuatorDisk)) operated at multiple
thrust coefficients to find the point where disk thrust equals hull
drag.

## Setup

- Domain 96 × 48 × 48 cells
- Wigley hull: L=48, B=10, T=6 cells (L/B=4.8, L/T=8)
- Half-submerged (waterline at z=NZ/2)
- Actuator disk: R = 1.5·T/2 = 4.5 cells (75% of draft, single-screw stern)
  at (hull stern + T/2, NY/2, waterline − T/2)
- Re_L = 5000, Fr = 0.25, ρ_w/ρ_a = 10:1
- WALE eddy-viscosity (Cw=0.5) added on top of vof.ν per cell
- 80 mom_step!s per C_T value; drag averaged over the last 25%
- Driver: [`scripts/wigley_self_propulsion_scan.jl`](./scripts/wigley_self_propulsion_scan.jl)

## Results

| C_T | thrust  | drag  | T − D     |
|----:|--------:|------:|----------:|
| 0.0 |   0.00  | 45.5  | −45.5     |
| 0.3 |   9.54  | 52.8  | −43.2     |
| 0.6 |  19.09  | 57.1  | −38.0     |
| 1.0 |  31.81  | 63.3  | −31.5     |
| 1.5 |  47.71  | 68.4  | **−20.7** |
| 2.5 |  79.52  | 73.3  | **+6.2**  |
| 4.0 | 127.23  | 78.3  | +49.0     |

Sign change between C_T = 1.5 and C_T = 2.5. Linear interpolation gives

**Self-propulsion C_T ≈ 2.27.**

At self-propulsion the disk thrust equals hull drag at ~70 cell-units.

## Physical observation: thrust deduction

Bare-hull drag (C_T = 0) is **45.5**. At self-propulsion (C_T ≈ 2.27)
the hull experiences drag ≈ **70** — a **+54%** increase.

This is the propeller-induced stern suction, commonly written as the
thrust deduction factor `t`:

```
t = (T − R_T) / T  where R_T is the bare-hull resistance
  ≈ (70 − 45.5) / 70  ≈ 0.35
```

Marine-engineering empirical correlations (Lurie & Taylor 1995;
Holtrop-Mennen 1982) give `t ≈ 0.18-0.25` for slender full-scale ships
at typical loadings. Our `t ≈ 0.35` is on the high side, expected
because (i) our Re=5000 is far below ship-scale Re~10⁹, (ii) the AD is
very heavily loaded (C_T=2.3 vs typical ship C_T~0.3), and (iii) the
disk is large relative to hull cross-section.

## What this validates

This is the first time the five Julia packages have produced a *number
that means something physically* together. The thrust-deduction
direction (drag rises with thrust) and order of magnitude match known
ship-CFD behaviour without us having calibrated anything.

## Caveats / not yet validated

- **Laminar regime.** Re=5000 is far below realistic ship-scale Re.
  Self-propulsion C_T will drop substantially at higher Re because hull
  Cf falls (~Re^−0.5 laminar, Cf~Re^−0.2 turbulent).
- **Wave-making dominant.** At Fr=0.25 with this hull, wave-making drag
  contributes much of the bare-hull 45.5 cell-units. Lowering Fr would
  drop both bare-hull drag and the absolute self-propulsion thrust.
- **Coarse grid.** 96×48×48 with hull L=48 cells gives ~6 cells across
  the hull beam — coarse for boundary-layer resolution.
- **Not benchmarked against OF.** OpenFOAM cross-validation is blocked
  on arm64 (no qemu-user-static).
- **Fixed thrust, not blade dynamics.** An actuator disk is a steady
  body force; a real propeller has time-varying loading.

## Path to release-blocker (Phase 3 release gate)

Per MASTER_PLAN §"What done means", the release gate is DTC
self-propulsion within ±15% of OpenFOAM and ±20% of el Moctar 2012
experimental data. Concretely:

1. Import DTC offsets → TabulatedHull (infra ready in ShipShapes.jl).
2. Run at Re ≥ 10⁶ with WALE (or upgrade to RANS k-ω SST).
3. Repeat this scan; compare interpolated C_T to OF.
4. Compare wake fraction `w` and thrust deduction `t` to el Moctar.

This pipeline is now end-to-end — only the inputs (DTC hull, higher Re)
are missing.

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
JULIA_NUM_THREADS=auto \
    WL_RE=5000 WL_NSTEPS_PER=80 WL_RPROP_FAC=1.5 \
    WL_CT_LIST="0.0,0.3,0.6,1.0,1.5,2.5,4.0" \
    julia +1.12 --project=. scripts/wigley_self_propulsion_scan.jl
```
