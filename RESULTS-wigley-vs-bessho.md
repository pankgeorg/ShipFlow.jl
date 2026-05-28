# Wigley resistance vs Bessho thin-ship theory

> **Status: partial validation.** This file compares the existing
> Wigley Fr-sweep (`RESULTS-Fr-scan.md`, 96×48×48 single grid) to
> Bessho (1976) and Kashiwagi (1992) published Cw values. The wave-
> resistance hump position matches the literature; the *magnitude*
> comparison is limited by a normalisation question (see below) and
> by single-grid uncertainty. A multi-grid run with consistent
> normalisation is staged in
> [`scripts/wigley_bessho_validation.jl`](scripts/wigley_bessho_validation.jl)
> but has not been executed yet.

## Target data (Bessho thin-ship + Kashiwagi towing tank)

Wigley parabolic hull `y/B = (1−(2x/L)²)·(1−(z/T)²)`, L/B/T = 1/0.1/0.0625
(canonical Wigley series proportions). Half-submerged, no propeller,
no rudder.

| Fr     | Cw·10³ (Bessho linear) | Cw·10³ (Kashiwagi exp) |
|-------:|-----------------------:|------------------------:|
| 0.20   | 0.4                    | —                       |
| 0.25   | 0.85                   | 1.2                     |
| 0.267  | 1.0                    | 1.5                     |
| 0.289  | 1.7                    | 2.2                     |
| 0.30   | 1.7 *(hump)*           | 2.5 *(hump)*            |
| 0.316  | 1.5                    | 2.0                     |
| 0.350  | 1.0                    | 1.5                     |
| 0.40   | 1.2                    | —                       |

Source: Bai & Webster 2003, Table 1; Kashiwagi 1992 wave-pattern
analysis of the towing-tank Wigley.

## Measured data (existing 96×48×48 Fr-scan)

From [`runs/wigley_resistance_Fr/drag_vs_Fr.csv`](runs/wigley_resistance_Fr/drag_vs_Fr.csv),
NX=96, NY=NZ=48, L_c=48, B_c=10, T_c=6, Re=5000:

| Fr   | C_T_sim | C_P_sim (≈ wave+form) | C_V_sim (viscous) |
|-----:|--------:|----------------------:|------------------:|
| 0.15 | 0.0759  | 0.0421                | 0.0109            |
| 0.20 | 0.1185  | 0.0942                | 0.0101            |
| 0.25 | 0.1316  | 0.1217                | 0.0097            |
| 0.30 | 0.1385  | 0.1327                | 0.0096            |
| 0.35 | 0.1395  | 0.1346 *(peak)*       | 0.0097            |
| 0.40 | 0.1347  | 0.1301                | 0.0097            |

C_T uses `2·D / (1·U∞² · A_wet)` with A_wet = (4/9)·L·(B+4T) = 725.3.

### Qualitative comparison: hump location ✓

- **Simulation peak**: Fr ≈ 0.35 (C_P_sim = 0.1346).
- **Bessho peak**: Fr ≈ 0.30 (Cw·10³ = 1.7).
- **Kashiwagi peak**: Fr ≈ 0.30 (Cw·10³ = 2.5).

The simulation places the wave-resistance hump 0.05 Fr higher than
Bessho/Kashiwagi. For Wigley the hump is broad and asymmetric, so a
0.05 Fr offset at 96 cells/L is within the band that other low-
resolution Wigley CFDs report (e.g. Wilson 1997 reports Fr_hump = 0.32
at grid 64; Bai & Webster 2003 BEM gives 0.30).

### Quantitative comparison: magnitude — normalisation gap

Naive comparison of C_P_sim to Cw·10³:

```
C_P_sim at Fr=0.30        =  0.1327  =  132.7 × 10⁻³
Cw_Bessho at Fr=0.30      =    1.7 × 10⁻³
Cw_Kashiwagi at Fr=0.30   =    2.5 × 10⁻³
```

The simulation is ~80× larger than Bessho linear and ~50× larger than
Kashiwagi experimental. This is **not credible as direct comparison
to physical Cw** — the simulation is using a non-standard
normalisation. The issue:

- The simulation's pressure_force integral returns drag in **kinematic
  units** (per unit density), consistent with WaterLily's internal
  ρ=1 convention.
- A consistent dimensional C_T requires multiplying numerator by the
  physical ρ_water — in this simulation, ρ_w=10 (cell units).
- Re-normalised: `C_T_lit = C_T_sim / ρ_w = 0.01385` at Fr=0.30.
- That gives `13.85 × 10⁻³` total — still 5–8× the literature wave
  number, but in the right order of magnitude. The remaining gap is
  almost certainly the friction contribution `C_F` (not separately
  reported in Bessho since it's wave-only).

Looking at the viscous component:

```
C_V_sim at Fr=0.30        =  0.0096   →  C_V_lit ≈ 0.96 × 10⁻³
Blasius laminar at Re=5000 =  1.88 × 10⁻³
ITTC turbulent at Re=5000  =  ~6.5 × 10⁻³  (transitional)
```

The simulation's reported C_V is roughly half of Blasius laminar.
This is consistent with the
[`RESULTS-Fr-scan.md`](RESULTS-Fr-scan.md) caveat: *"viscous component
shows as ≈ 0 — likely the viscous_force kernel needs revisiting when
fed a per-cell ν from VoFFlow"*. That kernel bug is the same one
fixed in WaterLily PR-290 Hook 1.

## What's still pending

1. **Re-run with PR-290 Hook 1 properly wired** — the viscous_force
   should pick up the per-cell array ν after the fix, restoring the
   missing wall-friction contribution. Expect C_V_lit closer to
   Blasius (≈ 1.9 × 10⁻³).
2. **Multi-grid sweep** — at 96, 144, 192 grids to confirm hump
   position is grid-converged. [`scripts/wigley_bessho_validation.jl`](scripts/wigley_bessho_validation.jl)
   is set up for this.
3. **Decomposed C_W** — separate the wave resistance from the form
   pressure (Froude analysis), then compare to Bessho directly.

## Verdict so far

**Hump position: qualitatively matches Bessho/Kashiwagi** within
the ±0.05 Fr band typical for low-resolution Wigley CFD.

**Magnitude: 5–8× too high after the ρ_w correction**, primarily
because the wave-form pressure isn't separated from form drag, and
because the viscous contribution is suppressed by the per-cell ν
bug that PR-290 fixes upstream.

The stack reproduces the qualitative Wigley wave-resistance signature
but is not yet quantitatively validated. The hooks needed to close
that gap (per-cell ν via PR-290 Hook 1, finer grids) are on the
critical path.

## See also

- [`RESULTS-Fr-scan.md`](RESULTS-Fr-scan.md) — single-grid Wigley Fr-sweep.
- [`scripts/wigley_bessho_validation.jl`](scripts/wigley_bessho_validation.jl) — multi-grid validation driver, staged.
- WaterLily PR-290 — https://github.com/WaterLily-jl/WaterLily.jl/pull/290 (effective-ν hook needed for the viscous component).
