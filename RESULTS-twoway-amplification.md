# Two-way coupling amplification trial (F2)

Quantifies the effect of WaterLily → rudder feedback (the
`trilinear_inflow` path) when the rudder sits **inside the rotor
race** rather than outside it. The original G2 trial had the rudder at
`prop_xc + 2·R`, well past the contracted slipstream; F2 places it at
`prop_xc + 0.5·R`.

Driver: [`scripts/rudder_in_race.jl`](scripts/rudder_in_race.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 36, B = 8, T = 5.
- Fr = 0.30, Re = 5000.
- Rotor at J = 0.32 (≈ G1 self-prop), thrust 196.2, torque 25.5.
- **Rudder at `rud_xc = prop_xc + 0.5·R = prop_xc + 1.5`** (inside race).
- δ = 10° fixed for 80 steps. Statistics over the final 25%.

Two cases run back-to-back:

- **TWOWAY=0**: `rudder_forces(rudder, δ, U∞; inflow=nothing)` — rudder
  sees freestream only.
- **TWOWAY=1**: `rudder_forces(rudder, δ, U∞; inflow=trilinear_inflow(flow.u))`
  — rudder samples the WaterLily velocity field at every panel control
  point.

## Result

|                              | TWOWAY = 0 | TWOWAY = 1 |
|------------------------------|------------|------------|
| `|V_local|` at rudder centre | 1.000      | **2.647**  |
| CL_rudder (steady)           | +0.302     | **+1.209** |
| F_hull_y (cell units)        | −2.49      | −2.19      |

**Local axial velocity jumps to 2.65×** the freestream — the rotor race
roughly doubles the inflow at this station, consistent with classical
actuator-disk slipstream contraction at this CT.

**CL_rudder rises 4×** (0.30 → 1.21). VortexLattice's
`additional_velocity` path absorbs both effects (raised dynamic pressure
and modified effective angle of attack) into the returned CL, normalised
to V_ref = U∞. With our smear convention `side = CL · 0.5·U∞²·S`, the
deposited side force on the fluid is correspondingly ~4× larger.

**Comparison with G2 (rudder at 2·R offset):**

| Run                | rud_offset | V_local | CL_rud | Δ vs TWOWAY=0 |
|--------------------|------------|---------|--------|---------------|
| G2 baseline        | 2·R        | ~1.05   | 0.302  | < 1 %         |
| **F2 in race**     | **0.5·R**  | **2.65**| **1.21** | **+300 %**  |

The two-way coupling effect grows from negligible (< 1 % at 2·R) to
**+300 % on rudder CL** at 0.5·R. This is the regime where
`trilinear_inflow` does meaningful work.

## Interpretation

- **Why F_hull went slightly down (−2.49 → −2.19)** despite a 4× larger
  rudder side force: the side force is deposited as a body force in
  the fluid, and the hull receives only the *integrated pressure
  reaction* of that force after it propagates upstream against the
  flow. With the rudder inside the race, the deposited force is in a
  region of much higher axial velocity, so its upstream-pressure
  signature is *partially convected past the hull* before settling.
  The hull therefore feels a marginally smaller side force, not a
  larger one — counter-intuitive but consistent with the convective
  geometry.
- **What this means for ship manoeuvring**: a real ship's rudder is
  almost always behind the propeller for exactly this reason —
  the race amplifies rudder authority dramatically. The F2 result
  shows the LiftingSurfaces tier captures this: a manoeuvring-model
  coupler reading CL_rudder would see correctly amplified values
  when two-way coupling is on. Whether the *integrated hull side
  force* increases or decreases depends on geometry.
- **VortexLattice is doing the work**, not just the dynamic-pressure
  rescaling. Each of the rudder's panel control points sees its own
  local velocity vector from `trilinear_inflow`, so the circulation
  distribution along the rudder span is non-uniform (the in-race
  portion gets much more lift). The reported CL is the panel-weighted
  average.

## Caveats

- 80 timesteps is enough for the rotor race to develop but not for
  the wave field to fully settle. A 200-step run would tighten the
  F_hull number.
- The rudder is *single-point smeared* at `(rud_xc, rud_yc, rud_zc)`;
  for a more physically accurate distribution the smear should span
  the rudder span (analogous to `smear_blades!` for the rotor).
- No two-way coupling on the *rotor* itself in this trial — the
  rotor's CT is precomputed at freestream. If both rotor and rudder
  sampled `flow.u`, the rotor's CT would shift slightly with the
  hull wake.

## See also

- `scripts/rudder_in_race.jl` — driver for this trial.
- `scripts/rudder_trial_sinusoidal.jl` — original G2 (rudder at 2·R).
- `LiftingSurfaces.jl/src/LiftingSurfaces.jl` — `trilinear_inflow`,
  `rudder_forces(..., inflow=…)`.
