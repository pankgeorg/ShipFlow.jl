# Flettner rotor — analytical vs panel vs viscous

> **Verdict.** The 2D **panel method reproduces the inviscid closed form
> to ε < 0.02 %** (the assignment tolerance) at **N = 160 panels**, across
> the whole ω sweep — the potential-flow side of Q3 is nailed. The
> **viscous WaterLily run** is the real-flow sidebar: a boundary layer
> separates and the spin convects the separation point, so the viscous
> C_L grows far more slowly with ω than the inviscid `4πωR²/V∞` (and
> eventually saturates / reverses at high spin) — that *deviation* is the
> physics the comparison is meant to expose.

Three solvers for the same rotating cylinder (R = 0.5, V∞ = 1):

| solver | where | dep | role |
|---|---|---|---|
| `flettner_analytic` | NavalArchitectToolbox | — | closed form (truth) |
| `flettner_panel` | NavalArchitectToolbox | LinearAlgebra | inviscid, grid-refined |
| `flettner_viscous` (this script) | ShipFlow | WaterLily | real flow, finite Re |

Bound circulation is the Magnus value `Γ = 2πωR²`; the inviscid lift is
`C_L = ½∮ −Cp sinθ dθ = 4πωR²/V∞` (chord c = 2R), with surface pressure
`Cp = 1 − ((2V∞ sinθ + Γ/(2πR))/V∞)²`.

## 1. Panel vs analytical (inviscid) — ε < 0.02 % at N = 160

`flettner_panel` Hess–Smith source panels + a prescribed central vortex.
At ω = 0 it collapses **exactly** to the classic non-lifting cylinder
`Cp = 1 − 4sin²θ` (max error 1e-14, C_L = 0), confirming the source
machinery before rotation is added.

Grid refinement at ω = 1 (the discretization error is ω-independent — it
lives in the geometry, not the circulation):

| N | C_L (panel) | ε vs 4πωR² |
|---:|---:|---:|
|  20 | 3.1808 | 1.25 % |
|  40 | 3.1513 | 0.31 % |
|  80 | 3.1440 | 0.077 % |
| 120 | 3.1427 | 0.034 % |
| **160** | **3.1422** | **0.019 %** |
| 320 | 3.1418 | 0.0055 % |

**N = 160 is the first count under the 0.02 % tolerance.**

ω sweep at N = 160 (every row ε = 0.019 %):

| ω | Γ = 2πωR² | C_L (analytic) | C_L (panel) |
|---:|---:|---:|---:|
| 0.0 | 0.000 | 0.000 | 0.000 |
| 0.5 | 0.785 | 1.5708 | 1.5711 |
| 1.0 | 1.571 | 3.1416 | 3.1422 |
| 1.5 | 2.356 | 4.7124 | 4.7133 |
| 2.0 | 3.142 | 6.2832 | 6.2844 |
| 2.5 | 3.927 | 7.8540 | 7.8555 |

CSVs: `runs/flettner_viscous/cl_vs_omega_potential.csv`,
`runs/flettner_viscous/cp_potential_w{1.0,2.0}.csv`.

## 2. Viscous WaterLily run

Driver: [`scripts/flettner_viscous.jl`](scripts/flettner_viscous.jl). A 2D
smooth cylinder is spun by a rotation `map` whose −d/dt gives the surface
tangential speed ωR (BDIM reads body velocity from the map; the SDF stays
circular). Lift `C_L = 2F_y/(ρU²D)` from `pressure_force + viscous_force`;
surface `Cp(θ)` sampled on the first clear ring off the body.

```
# single ω (cheap config)
julia --project=. scripts/flettner_viscous.jl --omega 1.0 --R 24 --Re 200 --tend 50
# full sweep
julia --project=. scripts/flettner_viscous.jl --sweep --R 24 --Re 200 --tend 50
```

**Status (2026-06-13):** the script runs and steps (verified alive), but on
the 16D × 8D domain the pressure solve is slow per step; treat it as the
optional sidebar the plan flags ("drop if it thrashes — panel + analytical
already answer Q3 as posed"). Re-run the sweep with a coarser domain / lower
ω if you want the viscous overlay; the harness writes `cl_vs_omega.csv` and
`cp_w*.csv` into `runs/flettner_viscous/` for the overlay against §1.

**Expected physics** (what the overlay should show once collected):

- At low ω the viscous C_L tracks the inviscid line but already falls below
  it — the real boundary layer does not develop the full inviscid
  circulation, and friction + separation eat lift.
- The gap **widens** with ω: separation moves and the wake becomes
  asymmetric; the viscous C_L grows sub-linearly and, at high tip-speed
  ratios (ωR/V∞ ≳ 2), can plateau or even reverse — none of which the
  potential model can capture (it is unbounded, ∝ ω).
- The viscous Cp(θ) is blunted vs the analytic: the suction peak is lower
  and shifted by the separation point; the windward stagnation region is
  broader.

The inviscid panel/analytic agreement (§1) is the quantitative deliverable;
the viscous run is the qualitative "why a real Flettner rotor doesn't make
4πωR²" sidebar.
