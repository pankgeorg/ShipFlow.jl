# RESULTS — Full-stack headline: ship + propeller + waves + turbulence

This is the project's headline integration: every package in the stack
running together on one Julia-WaterLily simulation. Compares the
ship-CFD goal sketched in [MASTER_PLAN.md](./MASTER_PLAN.md) against
what now actually runs.

## What runs

`scripts/wigley_propeller.jl`, 96 × 48 × 48 grid:

| Component        | Source                                  |
|------------------|------------------------------------------|
| Time integration | `WaterLily.mom_step!` with Hook 3 (`pois_tol=1e-8, pois_itmx=100`) |
| Hull body        | `ShipShapes.Wigley`  (analytic SDF, AutoBody) |
| Free surface     | `VoF.VoFFlow` (variable-density Poisson via Hook 1) |
| Turbulence       | `Turbulence.WALE` with per-cell ν₀ = vof.ν (Hook 1) |
| Propulsion       | `Propellers.ActuatorDisk` via WaterLily `udf` |
| Pressure         | Variable-coefficient multigrid (WaterLily MultiLevelPoisson) |

Numerical settings (smoke test):
- ρ_w/ρ_a = 10 (eased ratio — the gates are documented for 1000 too)
- Fr = 0.25, Re = 2000, U∞ = 1
- 60 momentum steps; total physical t ≈ 37 cell-time-units
- C_T = 0.3 (constant, no thrust-balancing loop)

## Numerical state at step 60

```
hull drag (pressure + viscous) = 51.3  cell-units
prescribed disk thrust          =  1.5
imbalance T - D                 = -49.8   (not self-propelling)
max wave amplitude              =  1.5 cells
|u|_max                         = 1.18 (slightly super-inflow)
```

The hull-drag / disk-thrust imbalance is BY DESIGN at this stage —
the smoke test runs with a fixed C_T = 0.3, far below the value
needed to balance hull drag at Re=2000. A self-propulsion calculation
(thrust ↔ drag balance) is the next gate; see HANDOFF.md §"Resume here".

## Headline meaning

This is the first time all five packages run together:
- **WaterLily** (the substrate)
- **Turbulence.jl** providing eddy ν₀
- **VoF.jl** providing the free surface
- **ShipShapes.jl** providing the hull
- **Propellers.jl** providing the disk

Each was validated independently against an external reference before
the integration:

| Package        | Standalone validation                                          | Headline RMS / max |
|---------------:|----------------------------------------------------------------|-------------------:|
| WaterLily      | cylinder Re=100 vs Williamson 1996, channel Re_τ=395 vs OF      | (substrate)        |
| Turbulence.jl  | Smagorinsky channel395: u_max within 0.5% of OF, RMS 0.028     | (LES eddy)         |
| VoF.jl         | damBreak vs Martin-Moyce 1952: RMS 3.8% at ρ=10:1              | (free surface)     |
| ShipShapes.jl  | Wigley resistance C_D_viscous = 0.0213 ≈ Blasius 0.0210         | (hull body)        |
| Propellers.jl  | Actuator disk vs 1D theory: U_disk -2.9%, U_wake -4.4%          | (propulsion)       |

## Outstanding work to Phase 3 / 0.1.0 release

1. **Thrust-balanced self-propulsion**: add a PI controller on C_T to
   match instantaneous hull drag.
2. **DTC hull**: tabulate offsets, build `TabulatedHull` (infra ready),
   swap in for Wigley.
3. **Higher-Re Layer-2**: WALE channel395 at N_HC=32, full t=400.
4. **MULES α-redistribution**: drop the 5% damBreak mass loss to 0.1%.
5. **OF Layer-2 cross-validation**: blocked on arm64 hosts.

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
JULIA_NUM_THREADS=auto \
    WL_NX=96 WL_NY=48 WL_NZ=48 \
    WL_L=48 WL_B=10 WL_T=6 \
    WL_NSTEPS=60 WL_CT=0.3 \
    julia +1.12 --project=. scripts/wigley_propeller.jl
```

## Update — self-propulsion controller trajectory

The PI controller (Kp=0.05, Ki=0.002, integral cap = 30) runs on a
96×48×48 Wigley smoke test (Re=2000, Fr=0.25, ρ=10:1) with 80-step
warmup. Trajectory over 300 steps:

| step | drag (smoothed) | thrust | err/T  |
|-----:|----------------:|-------:|-------:|
|   80 | 58.5 (warmup end, T seeded)         | 58.5  |  0.000 |
|  108 | 66.9            | 64.8   | +0.031 |
|  120 | 67.2            | 66.4   | +0.012 |
|  132 | 66.6            | 67.2   | **-0.009** |
|  168 | 56.7            | 62.4   | -0.099 |
|  216 | 31.7            | 40.6   | -0.279 |
|  252 | 12.5            | 21.0   | -0.684 |
|  300 |  7.8            |  6.8   | +0.133 |

**At step 132 the system reaches transient self-propulsion balance**
(err = −0.9%, within the 10% gate). However, the balance does not
persist: as the propeller wake speeds up the flow behind the hull,
the apparent hull drag *falls* (this is the marine-engineering
"thrust deduction" effect — flow acceleration by the propeller
reduces the effective hull resistance). The controller then chases
the falling drag, the system overshoots, and after a long swing
returns to a different operating point.

**This is a physical feature, not a controller bug.** Quantitative
self-propulsion needs:

- A higher Reynolds number where hull drag is large enough that the
  propeller-wake feedback is a small perturbation, *or*
- Time-averaging the controller signal over several wake-flow
  turnover times, *or*
- A parameter scan: run multiple fixed-C_T simulations to convergence
  and find the C_T where drag = thrust.

Documented; will revisit at Re=10⁴ with Turbulence.jl Smagorinsky.

## Update — self-propulsion via C_T parameter scan

The PI controller hypothesised that the propeller wake reduces apparent
hull drag ("thrust deduction"). A clean 7-point parameter scan at
Re=1000, Fr=0.25, ρ=10:1 settles the question:

| C_T | thrust | drag (steady) | T − D  |
|----:|-------:|--------------:|-------:|
| 0.0 |  0.00  | 53.00         | −53.00 |
| 0.1 |  0.90  | 54.25         | −53.35 |
| 0.2 |  1.81  | 55.29         | −53.48 |
| 0.3 |  2.71  | 56.16         | −53.45 |
| 0.5 |  4.52  | 57.67         | −53.15 |
| 0.8 |  7.24  | 59.13         | −51.89 |
| 1.2 | 10.86  | 59.49         | −48.63 |

**Drag rises monotonically with thrust** (the classic propeller-induced
suction at the stern). The PI controller's "drag fell with thrust"
observation was a transient artefact of incomplete steady-state.

**Self-propulsion would require C_T ≈ 6** at this hull / Re combination
— far outside the realistic AD range. The blocker is Re: at Re=1000
the friction coefficient is ~Cf=0.1, ~30× the full-scale Cf=0.003.
A meaningful self-propulsion run needs Re ≥ 10⁵ with the Turbulence.jl
LES wall model. That is the next prerequisite.

Scan committed at `runs/wigley_selfprop_scan/scan.csv`.
