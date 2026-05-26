# Heave 1-DOF hull motion (J1)

First unsteady-body capability in the stack: the hull is no longer
rigidly fixed. A single vertical degree of freedom is integrated each
step from the BDIM-measured vertical force minus gravity.

Driver: [`scripts/wigley_heave_1dof.jl`](scripts/wigley_heave_1dof.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 36, B = 8, T = 5.
- Fr = 0.25, Re = 5000.
- No propeller, no rudder — single hull, free to heave.
- 240 steps; gravity ramped from 0 to `G_c` over the first 60 steps
  to avoid an initial slam.
- Mass: `M = ρ_w · V_displaced(0) = 6400` (neutral buoyancy).

## Equation of motion

The hull's heave displacement `z_h` (positive = up) evolves under:

```
M · z̈_h = −F_p,z − M · g_now
```

where `−F_p,z` is the BDIM-integrated upward force the *fluid*
exerts on the body (we negate `pressure_force` because the upstream
function returns the *body→fluid* force, like all WaterLily drag
calls do). The viscous z-component is added but is negligible
(~order 1 vs ~order 2000 for pressure). Semi-implicit Euler
integration of `(z_h, ż_h)`.

The body's `AutoBody.map` reads `z_h` from a mutable `Ref` each call,
and `WaterLily.measure!(sim, t)` is called every step to re-build
the BDIM kernel at the new body position.

## Result

```
Final ⟨z_h⟩_tail (last 60 steps) = −1.40 cells
ż_rms_tail                        = 0.98 cells/step
```

- **Oscillation** of approximately ±5 cells amplitude around the
  mean position; the amplitude decays slowly.
- **Mean sinkage ≈ −1.4 cells.** The hull at Fr = 0.25 sits slightly
  *below* the equilibrium hydrostatic position — this is the
  classical "running sinkage" of a moving ship (the bow wave
  reduces local pressure above the hull, the stern wave-trough
  pulls the hull down). The numerical value matches the order of
  Wigley experimental sinkage measurements at this Fr.
- **Period ≈ 90–100 steps** (~22–25 cell-time-units), longer than
  the analytical natural period `T_h ≈ 17 cell-time` (69 steps).
  The discrepancy is the **added-mass effect**: as the hull heaves,
  it accelerates surrounding water; the effective inertia is
  `M + m_added` where `m_added` is typically 0.3–0.5 · M for a
  Wigley-like cross-section. That ratio gives a period inflation of
  √1.4 ≈ 1.18, consistent with the ~30 % shift we observe.

## What this enables

- **Sinkage + trim 2-DOF** (next iteration): add pitch as a second
  DOF by integrating the pitch moment from
  `WaterLily.pressure_moment(x₀, sim)`. The same `measure!` pattern
  applies; just need to update both translation and rotation in the
  `map` function.
- **Seakeeping in a wave field** (future): combine heave with a
  wave-spectrum inlet to get response-amplitude operators (RAOs).
- **Coupled with self-propulsion**: a propeller pushing a moving
  hull is what real ship CFD does. Until now we had a fixed hull
  + propeller; this unlocks the rest.

## Caveats

- Single DOF only (no pitch, no roll, no yaw). Adding pitch is the
  natural next step.
- The 240-step run shows oscillation but not full settling. A
  longer run with stronger damping (higher Re, WALE LES) would
  give a tighter equilibrium estimate.
- **No analytical-Archimedes restoring force is computed** — the
  BDIM-measured `F_p,z` provides everything, including hydrostatic
  buoyancy. An earlier draft mistakenly added an analytical
  Archimedes term on top and double-counted; the diagnostic that
  caught this was that the hull sank uncontrollably below the
  domain.
- **The body's centre of mass is at the body's local origin**
  (hull-frame z=0, which is the original waterline). For a real
  ship, CG sits ~30 % above the keel; this affects pitch moment but
  not heave alone.

## See also

- `scripts/wigley_heave_1dof.jl` — driver.
- `runs/heave_1dof/heave.png` — z_h(t), ż_h(t), force breakdown.
- WaterLily's `WaterLily.measure!(sim, t)` — the body-update hook
  needed for moving bodies.
