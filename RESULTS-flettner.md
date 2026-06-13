# Flettner rotor — analytical vs panel vs viscous

> **Verdict.** The 2D **panel method reproduces the inviscid closed form
> to ε < 0.02 %** (the assignment tolerance) at **N = 160 panels**, across
> the whole ω sweep — the potential-flow side of Q3 is nailed and is the
> quantitative deliverable. The **viscous WaterLily run is PARTIAL / not
> validated**: the static (ω=0) case is correct, but the rotating-body
> force extraction returns nonphysical lift (wrong sign, ~10× too large) —
> see §2. Per the plan it is the optional sidebar to drop if it thrashes;
> it is kept here as a documented starting point, not as a result.

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

**Status (2026-06-13): PARTIAL — not validated; the sidebar the plan says
to drop if it thrashes.** The script runs end-to-end and the *static* case
is correct — at ω = 0 it gives C_L ≈ 0 (−0.0007), confirming the
pressure/viscous force machinery and the Cp sampling. But with the spin
turned on the extracted lift is nonphysical:

| ω | C_L viscous (this script, R=8, Re=150) | C_L inviscid |
|---:|---:|---:|
| 0.0 | −0.001 | 0.000 |
| 1.0 | **−32.4** | 3.142 |
| 2.0 | **−58.4** | 6.283 |

Both the **sign and the magnitude are wrong** (a real Flettner C_L is a
modest positive fraction of the inviscid value, not ~10× larger and
negative). The fault is in the rotating-body force extraction, not the
flow solve: spinning the cylinder through a rotation `map` makes
`pressure_force`/`viscous_force` pick up the unsteady and rotational
pressure contributions on the immersed boundary, which the simple
`F = −(fp+fv)` reduction does not separate from the true aerodynamic lift.
On a smooth (rotation-symmetric SDF) body the spin also has to be imposed
purely through the BDIM velocity, and at R ≈ 8–16 cells that surface
velocity is under-resolved.

**This is left as a documented starting point, not a result.** Fixing it
needs one of: (a) a no-slip *surface-velocity* boundary condition that does
not move the SDF (so the gauge pressure stays clean), (b) subtracting the
rigid-body/added-mass pressure contribution before integrating, and (c) a
finer cylinder (R ≳ 32) once (a)/(b) are right. The inviscid panel/analytic
agreement (§1) is the quantitative Q3 deliverable; this viscous overlay is
the optional qualitative sidebar and is **not** trustworthy as written.

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
