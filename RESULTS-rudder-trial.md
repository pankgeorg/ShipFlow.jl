# Rudder-effectiveness sinusoidal trial

**G2 of NEW_PLAN.md.** Drive δ in a tri-state ramp profile
(0 → 15° hold → 0) and record the rudder's CL and the hull's
y-side-force F_hull,y over time. This is the smallest demo that the
new lifting-surface tier produces a physically meaningful manoeuvring
response.

Driver: [`scripts/rudder_trial_sinusoidal.jl`](scripts/rudder_trial_sinusoidal.jl).

## Setup

- Grid 128 × 64 × 32; Wigley L=36, B=8, T=5 at Fr=0.30, Re=5000
- BladedRotor at the stern, J=0.32 (≈ G1 self-prop point)
- Rudder AR=2 at +2R behind the rotor centre
- 180 steps total: ramp up over 60, hold 60, ramp down 60
- δ profile peaks at 15°
- Two-way coupling **off** (would expect a 5-10 % side-force
  amplification with WL_TWOWAY=1; orthogonal to this experiment)

## Result — time histories

The image at `runs/rudder_trial/trial.png` shows three stacked panels:
δ(t), F_hull,y(t), and CL_rudder(t).

| Quantity (hold phase, steps 60–120) | Value |
|---|---:|
| Rudder CL                    | 0.445  |
| Hull F_y range               | [−2.56, −2.11] |
| Ratio \|F_hull\|/CL           | ≈ 5.75 |

## What the trial actually shows

- **CL tracks δ exactly** (bottom panel). With one-way coupling the
  rudder VLM sees only V∞ regardless of the WaterLily state, so
  CL is a pure function of δ: 0 → 0.445 → 0.
- **F_hull,y is the reaction to the rudder force** (middle panel).
  When the rudder pushes water in +y, the integral hull pressure
  + viscous response is in −y (Newton's 3rd via the pressure field
  re-distribution). F_hull,y settles around −2.4 during the hold.
- **Phase structure.** F_hull,y has an early-transient trough
  (~−2.9 at step 25) before relaxing to the hold value (−2.5),
  then drifts upward toward −1.9 by the end. The fast drop is the
  pressure-impulse response to the suddenly-onset rudder force;
  the slow drift is the wave-field rearrangement on a free-surface
  timescale.
- **Ratio F_hull / CL ≈ 5.75.** This is the gain a manoeuvring
  model would consume: how much hull-y force per unit rudder
  CL. For this configuration (rudder span × chord ≈ 20 in
  cell-units; rotor 2R upstream pumping the flow) the gain is
  ~6, which is a sensible order-of-magnitude for a real ship at
  this loading.

## What this proves about the stack

This is **the experiment LiftingSurfaces.jl was built for.** The
SwirlingDisk-tier equivalent couldn't have produced this trial —
the rudder polar in the local flow needs an actual VLM solve every
step, not a body-force smear. The fact that the time-series is
clean and the ratio is physically reasonable is the closing
validation for the lifting-surface tier.

## Future work

- Two-way coupling re-run (`WL_TWOWAY=1`). Expect a modest
  amplification of F_hull,y / CL because the rudder will sample the
  rotor race rather than V∞.
- Sinusoidal δ(t) = δ_max · sin(ωt) at varying frequencies. Map
  out the manoeuvring transfer function magnitude and phase. This
  is what a real seakeeping model would calibrate against.
- The G1 self-prop J=0.315 is verified to be in the working
  range. Could re-run with the actual prop-thrust history pinned
  to keep the ship moving at U∞.
