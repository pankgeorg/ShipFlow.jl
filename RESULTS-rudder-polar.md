# Rudder polar via VLM (LiftingSurfaces.jl)

Sweep of rudder angle δ from -25° to +25° in 5° increments, all via
`LiftingSurfaces.rudder_forces`. Rectangular flat-plate rudder,
AR = 2 (chord = 1, span = 2), 16 spanwise × 8 chordwise panels.

Driver: [`scripts/rudder_polar_via_VLM.jl`](scripts/rudder_polar_via_VLM.jl).

## Results

| δ (°) |       CL |       CD |    L/D |
|-----:|---------:|---------:|-------:|
| −25  |  −1.0042 |  0.16841 |  −5.96 |
| −20  |  −0.8243 |  0.11163 |  −7.38 |
| −15  |  −0.6308 |  0.06452 |  −9.78 |
| −10  |  −0.4267 |  0.02923 | −14.60 |
|  −5  |  −0.2152 |  0.00739 | −29.11 |
|   0  |   0.0000 |  0.00000 |    Inf |
|  +5  |  +0.2152 |  0.00739 | +29.11 |
| +10  |  +0.4267 |  0.02923 | +14.60 |
| +15  |  +0.6308 |  0.06452 |  +9.78 |
| +20  |  +0.8243 |  0.11163 |  +7.38 |
| +25  |  +1.0042 |  0.16841 |  +5.96 |

**dCL/dα at linear range (±5°)** = 2.466 / rad.
Prandtl-LLT elliptic-loading prediction for AR = 2 is 2π·AR/(AR+2) =
3.142 / rad. The VLM/LLT ratio is 0.785, which agrees with the
classical rectangular-planform correction factor of ≈ 0.78.

CD is exactly zero at δ = 0 and grows symmetrically as CL² beyond
that — induced drag, as expected for an inviscid lifting-surface
solution (no skin friction in VLM).

## Interpretation

- VLM is solving correctly: lift is monotonic in δ, drag is even, both
  go to zero at δ = 0, and the linear-range slope matches
  rectangular-planform theory to within numerical noise.
- The `Rudder` API in `LiftingSurfaces.jl` correctly invokes
  VortexLattice.jl: it accepts kwargs `δ` and `V∞`, returns a
  named-tuple with `CL`, `CD`, `CY`, `CM` in the right frames.
- `CY = 0` everywhere because the rudder is built with span along
  +y (body-frame), so the lift acts in body-z, not body-y. The
  user is responsible for mapping VLM body-frame to the ship's
  world frame.

## Use case

This polar is the **lookup table** a maneuvering-model coupler would
need: given a commanded δ, return the resulting side force and
yaw moment on the ship. In the full coupled run
(`scripts/rudder_in_freestream.jl`), the same call is made every
WaterLily step with VortexLattice's `additional_velocity` keyword
set to a closure that samples the local WaterLily velocity at each
panel control point — so VLM sees the actual propeller-race +
hull-wake inflow rather than just V∞.
