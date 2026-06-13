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

> **Progress note (2026-06-12) — DTC bare-hull resistance, step 1 of Layer 3.**
> First resistance run on the *real* DTC geometry (full five-package
> stack: WaterLily+VoF+Turbulence+ShipShapes, tabulated SDF). Stack
> integrates and is stable/mass-conserving, but the calm-water
> resistance at Fr=0.218 **does not pass** the ±15 % gate: C_P,sim ≈
> 3.6e-3 (NL=128) vs reference residuary C_R = 0.62e-3 (+475 %). The
> resolved force is dominated by a trapped longitudinal sloshing mode
> (C_T ±56 % std, never settles) and BDIM form-drag bias; the Froude
> subtraction C_T−C_F,ITTC is invalid because BDIM develops ~no wall
> friction (no y⁺ wall function — exactly the risk-register item). The
> *fully-submerged* DTC drag IS grid-converged (C_P,deep=5.25e-3,
> split-half 1.6 %), confirming the geometry/solver integration is
> sound. Details + next steps: [RESULTS-dtc-resistance.md]. Bare-hull
> resistance must settle (kill sloshing, add near-wall stress model)
> before the DTCHullProp self-propulsion gate below is reachable.
>
> **Progress note (2026-06-13) — DTC resistance round 2.** Acted on the
> three round-1 root causes. **Fix 1 (sloshing): PARTIAL.** Built a
> sponge/wave-damping zone; the key result is a numerics finding — a
> Rayleigh *velocity* body force in `mom_step!`'s `udf` stalls the
> multigrid Poisson (nit 50 vs 2–3) on the two-phase 100:1 system, so the
> velocity relaxation must be applied **post-projection** (`post_sponge!`).
> That + α (surface) damping reduces the swing and is Poisson-safe, but
> the confined-box longitudinal **seiche** persists (split-half ≈10 %,
> not the <2 % target) — reported as a windowed mean ± std. **Fix 2
> (Froude walk): BLOCKED** on a settled C_T; also corrected the
> reference — the el Moctar 2012 towing-tank envelope ends at **Fr=0.218**
> (no 0.28/0.33 truth), so the grounded walk is down the tested ladder.
> **Fix 3 (near-wall stress): WIRED + smoke-tested** — Turbulence.jl's
> Spalding BDIM wall function (SA path, water-gated) runs end-to-end from
> the DTC SDF but is not yet exercised at scale (gates on Fix 1).
> Next: settle the seiche (longer upstream fetch for a real inlet sponge,
> or an integer-period averaging window / steady local-time-stepping),
> then exercise Fix 3 and walk the Froude ladder. Details:
> [RESULTS-dtc-resistance.md] §Round 2.
>
> **Progress note (2026-06-13) — propeller actuator-disk Layer-2 PASS
> (step 2, "option D").** Chosen because it is *independent* of the
> unsolved DTC free-surface seiche above. The in-grid actuator disk now
> round-trips the validated DTMB-4382 open-water VLM: feed the VLM's
> radial loading (dT/dr, dQ/dr) into a new `Propellers.GradedDisk`, run
> single-phase in WaterLily, integrate the resolved thrust+torque, and
> recover KT/KQ/η within **<0.3 %** at J = 0.6, 0.889, 1.1 — well inside
> the ±10 % Layer-2 gate. The tightness is partly structural (in a
> confined periodic duct the control-volume momentum/angular-momentum
> balance is a conservation identity), so the headline is that the
> *whole pipeline is self-consistent* — VLM coefficients → radial loading
> → cell-unit (J,n,D,ρ) conversion → body-force deposition → resolved-flux
> recovery — and the in-grid axial/swirl profiles reproduce the VLM's
> bell-shaped radial loading. Ladder 1 (analytic): the uniform disk
> matches Froude–Rankine induced velocity to ~2 % across C_T=0.2–1.5 and
> shows the slipstream contraction. Ladder 3 (OpenFOAM `propeller`
> tutorial): identified as a *resolved rotating-mesh* propeller (1500 rpm
> solidBody, snappyHexMesh, createBaffles), **not** an actuator disk and
> with no published reference KT/KQ — so no like-for-like and no
> documented-value compare; a rerun is infeasible on this aarch64 box
> (OpenFOAM images are amd64-only, the Phase-0 qemu blocker). The
> defensible external anchors for the disk are therefore momentum theory
> + the experiment-matched VLM. This **unblocks** the actuator-disk side
> of the DTCHullProp self-propulsion gate below: the propeller model is
> validated; what remains is the hull resistance (the seiche) and wiring
> the disk behind the resolved DTC. Details:
> [RESULTS-propeller-layer2.md]; `Propellers.jl` PLAN milestone 2.

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

## Aside — 8401 lifting-device tools (capability additions)

Outside the five-package ship-CFD arc, the stack also carries reusable
lifting/circulation tools for the NTUA 8401 themes (see
`PLAN-lifting-devices.md`): `LiftingSurfaces.Wing` (finite-wing VLM,
validated vs lifting-line), `NavalArchitectToolbox.flettner_panel` /
`flettner_analytic` (2D rotating-cylinder potential flow, ε<0.02 % vs the
closed form at N=160), and `scripts/flettner_viscous.jl` here (the
WaterLily real-flow comparison; see `RESULTS-flettner.md`). NAT is the
unified surface (`using NavalArchitectToolbox` → propeller VLM + wing VLM +
Flettner panel); the viscous run stays in ShipFlow so NAT keeps no WaterLily
dep.
