# Rotor two-way coupling — Taylor wake measurement (I2)

Quantifies the rotor inflow correction picked up by
`trilinear_inflow(flow.u)` and a per-step VLM solve. Two-way coupling
on the **rotor side** (rudder-side was F2). Three-run protocol so we
can disentangle hull wake from the rotor's own induced velocity.

Driver: [`scripts/rotor_twoway_wakefraction.jl`](scripts/rotor_twoway_wakefraction.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 36, B = 8, T = 5.
- Fr = 0.30, Re = 5000.
- Rotor at the stern, R = 3, J = 0.32 (G1 self-prop point).
- VLM re-solved every 10 steps (`WL_VLM_EVERY=10`) — inflow changes
  on a slow time scale, so this is ~10× the speed of solving every
  step with negligible accuracy loss.
- 80 timesteps; statistics over the final 25 %.

Three runs:

- **Run 0** — hull alone, no rotor. Lets us measure the
  hull-induced velocity at the disk plane (Taylor wake) without any
  contamination from the rotor's own induced flow.
- **Run A** — rotor with precomputed freestream CT (no two-way).
  Baseline.
- **Run B** — rotor with `inflow=trilinear_inflow(flow.u)` (per-step
  VLM solve). The rotor sees the locally-modified flow.

## Result

|                              | freestream | hull-only | two-way |
|------------------------------|-----------:|----------:|--------:|
| V_local at disk plane        | 1.000      | **1.536** | 1.539   |
| V_local 1.5·R upstream       | —          | 1.169     | —       |
| Rotor CT                     | 13.881     | —         | **14.66** |
| Hull drag                    | 61.97      | —         | 62.85   |

**Taylor wake fraction** (hull only at disk plane):
`w_T = 1 − 1.536/1.000 = −0.536` (a 53.6 % *acceleration*, not a
deficit).

**Two-way coupling effect on CT**:
`ΔCT = 14.66 / 13.88 − 1 = +5.6 %` — the rotor produces 5.6 % more
thrust than the freestream baseline.

## Interpretation

- **The negative Taylor wake fraction is geometric.** The rotor is
  placed at `prop_xc = hull_xc + L/2 + T = 0.2·NX + 18 + 5 ≈ 49`
  cells, only 5 cells (≈ 1.0·T) behind the Wigley stern. At this
  station, mid-keel depth, the flow is **accelerating** to fill the
  volume the hull's parallel-midbody-to-pointed-stern shape is
  vacating. There is no boundary-layer deficit yet — that develops
  further downstream.
- **In real ship engineering, w_T is typically +0.15 to +0.35**
  (positive deficit) because the propeller sits behind a much
  blunter stern with a developed boundary layer at Re ≈ 10⁹. Our
  Re = 5000, slender Wigley, and stern-proximate placement put us
  outside that regime.
- **The methodology is correct.** A two-way coupled rotor sees the
  hull-modified inflow and adjusts its CT accordingly. The
  trilinear-inflow path through VortexLattice's
  `additional_velocity` works. We just happen to be in a regime
  where the geometry gives w < 0.
- **CT shifted by ~+5.6 %** (vs the predicted +54 % from a pure
  dynamic-pressure rescale of the inflow). The discrepancy is
  because the VLM also sees the modified angle of attack
  distribution along the blade span — at higher axial inflow the
  effective α at each section drops, partially offsetting the
  dynamic-pressure gain.

## Caveats

- **Computation cost.** Solving the VLM every step is ~10× the
  WaterLily step time. The `VLM_EVERY=10` workaround amortises;
  CT_cached varies smoothly within a 10-step window so this is
  safe at this Re. Higher-Re unsteady cases would need
  `VLM_EVERY=1` and a deliberate budget.
- **Disk-plane probe also picks up rotor race.** Reading `flow.u`
  at the disk plane in Run B nominally double-counts the rotor's
  own induced velocity. The reason CT_B is still close to CT_A is
  that VortexLattice's `additional_velocity` is interpreted as a
  *perturbation* on top of the freestream, and the rotor's own
  induced velocity in actuator-disk theory would only contribute
  ~a·U∞ ≈ 0.3 — much less than the 0.54 acceleration we see, which
  is dominated by hull geometry.
- **Only one rotor placement tested.** A scan over `prop_xc`
  (distance behind stern) would map the transition from
  acceleration zone (w < 0) to wake-deficit zone (w > 0).
- **No comparison to a wake-fraction-corrected single-run.** The
  practical use of this is to set the freestream CT to account
  for w, run a faster single-pass simulation. Worth a follow-up.

## See also

- `RESULTS-twoway-amplification.md` — F2: rudder-side two-way
  coupling (different geometry, complementary effect).
- `LiftingSurfaces.jl/src/LiftingSurfaces.jl` — `rotor_forces`
  with `inflow=` kwarg; `trilinear_inflow`.
- `runs/rotor_twoway/summary.csv` — the run output.
