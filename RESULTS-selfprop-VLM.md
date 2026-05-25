# Self-propulsion via the VLM BladedRotor

Replay of the `wigley_self_propulsion_scan.jl` experiment using the
new lifting-surface tier (`LiftingSurfaces.BladedRotor`) instead of
`Propellers.ActuatorDisk`. **G1 of NEW_PLAN.md.**

Driver: [`scripts/wigley_self_propulsion_VLM.jl`](scripts/wigley_self_propulsion_VLM.jl).

## Setup

- Wigley hull L=48, B=10, T=6 at 96×48×48
- Fr=0.25, Re=1000, ρ_w/ρ_a=10
- 3-blade VLM rotor at the stern: R=0.8·T/2=2.4, R_hub=0.2R,
  chord(0.25R→0.18R), twist(35°→15°), 12×4 panels per blade
- Per J: rotor solves once, CT used to size the thrust smear at
  the disk centre (`smear_force!`, ε=2.5); CQ deposited as a
  tangential ring at r=0.7R (`smear_torque!`, N=8, ε=2.0)
- Each run: 100 mom_steps; drag averaged over the last 25%

## Results

| J    | \|CT\|  | Thrust   | Hull drag | T − D     |
|-----:|-------:|---------:|----------:|----------:|
| 0.25 | 21.98  | 198.91   | 143.50    | **+55.40** |
| 0.30 | 15.36  | 138.96   | 130.76    |  **+8.21**  |
| 0.35 | 11.35  | 102.69   | 121.54    |  −18.85   |
| 0.40 |  8.74  |  79.07   | 114.47    |  −35.40   |
| 0.50 |  5.65  |  51.16   | 104.54    |  −53.39   |
| 0.70 |  2.94  |  26.62   |  93.34    |  −66.71   |
| 1.00 |  1.48  |  13.38   |  85.93    |  −72.55   |

**Self-propulsion J ≈ 0.3152** (linear interpolation of T−D
zero-crossing between J=0.30 and J=0.35).

## Interpretation

- **Drag is monotonically decreasing in J** — the rotor's
  upstream-induced velocity accelerates the flow in the
  hull-stern region, reducing the effective hull boundary-layer
  contribution to drag. At J=∞ (zero rotor input) the hull drag
  would asymptote near ~85, the "naked-hull" value.
- **CT ∝ 1/J²-ish** as expected from VLM kinematics; the rotor
  has to spin much faster at low J to make sense of the
  geometric pitch.
- **Self-prop at J=0.32 vs C_T≈2.33** in the earlier
  SwirlingDisk scan (`RESULTS-selfprop.md`). These aren't
  directly comparable — different rotor configurations and
  different Re — but the qualitative result lines up:
  self-propulsion needs a fully-loaded rotor at this hull.
- **Sign change is clean and monotone** — no oscillations near
  the balance point, suggesting the BladedRotor + WaterLily
  coupling is dynamically well-behaved.

## Comparison: BladedRotor vs SwirlingDisk

Both reach self-propulsion in the same monotone way, but the
qualitative wake structure differs (see `RESULTS-bladed-vs-swirl.md`).
The BladedRotor's concentrated-thrust smear gives a more localised
hull-suction effect, which is why its self-prop J might be slightly
lower than the SwirlingDisk's equivalent J (the local stern
acceleration is stronger per unit thrust).

## What to do with this number

J = 0.315 is the design operating point if we want this rotor at
this hull to self-propel. Going forward:

1. **G2 (next):** rudder-effectiveness trial at this J. Now we
   know the rotor RPM that matches the hull, the rudder sits in a
   properly-loaded propeller race.
2. **G4 (later):** sectional smear comparison — re-do this sweep
   with a per-panel smearing of BladedRotor's force, see if
   self-prop J shifts toward 0.40-0.50 (closer to SwirlingDisk's
   distributed-force result).
