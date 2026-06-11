# MASTER PLAN — Ship-flow CFD on WaterLily.jl

The ultimate target is *flow around a ship with a rotating propeller in
water* — hull resistance, wave-making, and propeller wake on a single
unified solver. This document sequences the five workloads that get us
there, defines the decision gates between them, and points to the
shared validation harness in [HARNESS.md](./HARNESS.md).

## The five workloads

| # | Repo                                              | Plan file                                                                                       |
|---|---------------------------------------------------|-------------------------------------------------------------------------------------------------|
| 1 | `pankgeorg/WaterLily.jl` (branch `plan/upstream-hooks`) | [PLAN-upstream-hooks.md](https://github.com/pankgeorg/WaterLily.jl/blob/plan/upstream-hooks/PLAN-upstream-hooks.md) |
| 2 | `pankgeorg/Turbulence.jl`                         | [PLAN.md](https://github.com/pankgeorg/Turbulence.jl/blob/main/PLAN.md)                          |
| 3 | `pankgeorg/VoF.jl`                                | [PLAN.md](https://github.com/pankgeorg/VoF.jl/blob/main/PLAN.md)                                 |
| 4 | `pankgeorg/Propellers.jl`                         | [PLAN.md](https://github.com/pankgeorg/Propellers.jl/blob/main/PLAN.md)                          |
| 5 | `pankgeorg/ShipShapes.jl`                         | [PLAN.md](https://github.com/pankgeorg/ShipShapes.jl/blob/main/PLAN.md)                          |

## Dependency graph

```
                            WaterLily.jl  (upstream)
                                  ▲
                                  │
                      ┌───────────┴──────────────┐
                      │  PLAN 1: upstream hooks  │
                      └───────────┬──────────────┘
                                  │
        ┌──────────────┬──────────┴────────┬──────────────┐
        ▼              ▼                   ▼              ▼
  Turbulence.jl     VoF.jl           Propellers.jl   ShipShapes.jl
  (needs Hook 1+2)  (needs 1+2+3)    (needs 3, opt)  (needs nothing)
        │              │                   │              │
        └──────────────┴───────┬───────────┴──────────────┘
                               ▼
                         ShipFlow.jl
                  (this repo — application + harness)
```

ShipShapes.jl and Propellers.jl (using `udf`) are independent of the
upstream hook PR — they can ship in parallel as Day-1 work.

## Sequencing — phases, not milestones

### Phase 0 — the foundation (weeks 1–2)

In parallel:

- Spike ShipShapes.jl Wigley analytic SDF + L1 tests. Trivially shippable;
  unblocks every later flow test.
- Spike Propellers.jl uniform actuator disk via `udf`. Validate against
  1D actuator-disk theory in a uniform stream.
- Set up the validation harness (see [HARNESS.md](./HARNESS.md)): the
  OpenFOAM Docker image, the reference-data repo, the comparison script
  skeleton.

**Decision gate.** Are the analytic checks passing? Is the OpenFOAM
container reproducible? If yes, continue. If the harness is the
bottleneck, fix that before writing any more physics.

### Phase 1 — turbulence (weeks 3–8)

- Implement WaterLily upstream Hook 1 (effective viscosity) on the
  `plan/upstream-hooks` branch. Run WaterLily's own tests with array-`ν`
  set to a scalar; bit-identical results required.
- Build Smagorinsky LES in Turbulence.jl against the branch.
- Hit Layer 1 (decaying-HIT spectrum) and Layer 2 (channel395) for
  Smagorinsky.
- Add WALE.
- Implement Hook 2 (`transport!`) upstream.
- Add Spalart–Allmaras → k–ω SST in Turbulence.jl.

**Decision gate.** Does channel395 RANS+LES match OpenFOAM within ±5%
on u⁺? If yes, open the upstream WaterLily PR (Hook 1 + Hook 2 only,
not Hook 3 — see below). If not, decide whether the BDIM-wall
treatment is salvageable or whether the wall layer needs a different
approach (e.g., immersed boundary with explicit wall stress).

### Phase 2 — free surface (weeks 9–18)

- Build VoF.jl on the upstream-hooks branch (Hook 2 is what α
  advection needs).
- Hit Layer 1 (Zalesak, dam break, Hysing bubble).
- Hit Layer 2 (OpenFOAM `damBreak`) within ±10%.
- Add ShipShapes KCS/DTC tabulated SDFs.
- Combined VoF + ShipShape integration test: static floating box
  reaches Archimedes equilibrium.

**Decision gate.** Does the Poisson solver still converge at ρ-ratio
1000? Does mass drift stay under 0.1% over 5 seconds? If both yes,
proceed. If not, decision: switch to a pressure decomposition
(`p_rgh`) or stop and fix the multigrid coarsening at the interface
before adding propellers.

> **GATE PASSED — 2026-06-11.** damBreak ρ=1000, N=128, 5 s through the
> `density_coefficient!` Poisson path: mass m/m₀ = 1.0000 (both
> clamp+mass_repair and MULES), Poisson at ≈2 V-cycles/projection
> throughout. No `p_rgh` needed. Caveat discovered en route: MULES
> without an interface-compression flux homogenizes α over long runs —
> use `step_vof!`+`mass_repair` until the `cAlpha`-style term lands.
> Details: [RESULTS-damBreak.md](./RESULTS-damBreak.md) §Phase-2 gate.

### Phase 3 — propulsion (weeks 19–24)

- Finalize Propellers.jl actuator-disk + actuator-line.
- Hit Layer 2 (`incompressibleFluid/propeller`) within ±10%.
- Hit Layer 3 (DTCHullProp) within ±15% — *this is the headline result*.
- Decide on Hook 3 (composable body forces): if the user-facing API
  needs it, upstream it; otherwise stay with `udf` and keep the PR
  small.

**Decision gate.** Is the self-propulsion point predicted within ±15%
of OpenFOAM and ±20% of experiment? If yes, ship a v0.1 of every
package and open the upstream PR. If not, debug — most likely
suspects: actuator-disk smearing width, BDIM wall treatment under the
hull, or a missing free-surface refinement around the wake.

### Phase 4 — public release (week 25+)

- Make the five repos public, possibly request transfer to the
  WaterLily-jl org.
- Land the upstream WaterLily PR.
- Tag 0.1.0 on all five.
- Write the paper.

## Cross-cutting decisions (made once, applied everywhere)

| Decision                          | Choice                                                                  |
|-----------------------------------|-------------------------------------------------------------------------|
| Default numeric precision         | Float32 for fields; Float64 only for `ω` if range demands it            |
| Backend matrix                    | CPU (threads + SIMD) and CUDA mandatory; AMDGPU best-effort             |
| Reference data hosting            | `pankgeorg/cerulean-reference-data` (private Git LFS repo, created lazily) |
| OpenFOAM version                  | OpenFOAM 13 / OpenFOAM-dev (matches what we cloned)                    |
| OpenFOAM image                    | `opencfd/openfoam-default:2406` *or* build from the cloned `OpenFOAM/` |
| Per-PR CI budget                  | < 10 minutes total across all packages                                  |
| Nightly OpenFOAM CI budget        | < 4 hours total                                                         |
| Coordinate convention             | x forward, y starboard, z up (right-handed)                             |
| Force / moment sign convention    | Force on fluid is positive; force on body is the negative               |

## Risk register

| Risk                                                            | Mitigation                                                  |
|-----------------------------------------------------------------|-------------------------------------------------------------|
| Upstream rejects the hook PR                                    | Maintain the hooks on `pankgeorg/WaterLily.jl` indefinitely; downstream packages depend on the fork |
| BDIM wall treatment is inadequate for high-Re turbulent flows   | Plan B: switch to a true immersed boundary with explicit wall-stress (Kempe, Fröhlich 2012) — adds ~weeks |
| Uniform Cartesian grid can't resolve ship + wake + propeller    | Plan B: limit Phase 3 scope to bare hull + actuator disk; defer blade-resolved propeller indefinitely |
| Multigrid stalls at large ρ-ratio                               | Plan B: GMRES-preconditioned-by-multigrid; or 2-phase-aware coarsening |
| OpenFOAM container is non-reproducible                          | Pin a single image SHA; build from source as fallback       |

## What "done" means at the milestone

A user runs:

```julia
using WaterLily, ShipShapes, Turbulence, VoF, Propellers, ShipFlow

ShipFlow.dtc_self_propulsion(
    scale = 1/59.4,
    Fr    = 0.218,
    n     = 9.5,    # rev/s
    turb  = :komega_sst,
)
```

…and gets resistance, thrust, wake fraction, and a wave-pattern PNG, in
under an hour on a single GPU, with results within ±15% of OpenFOAM and
±20% of the el Moctar 2012 experimental data.

That's the bar.

## Iteration cadence — how to actually work this plan

- One workload "in motion" at a time. The user (me) is one person; pretending parallel work helps just creates context-switch cost.
- When a workload's PLAN milestone goes green, commit + push, then move
  to the next workload. The plan files are living documents — update
  them in place when reality disagrees.
- Don't open the upstream PR until Phase 2 ends. The "don't lick the
  cookie" rule: the upstream issue must arrive with code, tests, and
  two real downstream consumers (Turbulence.jl and VoF.jl) already
  using it.
- A weekly checkpoint (15 minutes) just reads the milestone table in
  each PLAN and asks: did the most recent gate pass? If not, what's
  blocking? No rewrite of plans unless a gate flips.
