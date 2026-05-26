# Wigley above-waterline deck SDF (L1)

Addresses the K1 finding (`RESULTS-heave-pitch-2dof.md`) that the
Wigley pitches unstably because the SDF is truncated at body z = 0,
leaving no above-waterline body for the lifted bow to push water
with. Extends `wigley_sdf` (and the `Wigley` constructor) to support
a `deck_h` parameter.

## Implementation

```julia
Wigley(; L, B, T, deck_h=0, map=(x,t)->x)
```

`deck_h > 0` extends the hull above the waterline with a vertical-
walled prism. The deck cross-section is identical to the parabolic
waterline cross-section, so the body is C0 + C1 continuous at z = 0
(both the below-waterline parabolic surface and the deck wall are
vertical at z = 0).

**Bug in v1**: the deck interior was initially treated with the
same Eikonal formula as the parabolic body, which returns the
y-side-wall distance only. That over-extends the BDIM kernel
*upward* (the SDF claims the body is much taller than it is) and
produces spurious pressure forces. K1 with the buggy deck gave
worse instability than no deck (M_y = 10500 vs 7000).

**v2 (committed in `ShipShapes.jl/8ece26c`)**: the deck interior
returns `−min(d_top, d_side)` — true minimum-distance to either
the top deck face (`deck_h − z`) or the parabolic side wall
(Eikonal-corrected). BDIM now sees the body's actual extent.

## ShipShapes tests

7 new tests covering:
- bit-identical behaviour at `deck_h = 0` (back-compat)
- interior points inside the deck region return negative SDF
- side and above-deck points return positive SDF
- `Wigley(; deck_h=h)` factory propagates the parameter

Total: 34 → 41 tests, all pass.

## Effect on K1 (heave + pitch 2-DOF)

| Setup                          | ⟨z_h⟩  | ⟨θ⟩      | M_y (saturated) |
|--------------------------------|--------|----------|-----------------|
| No deck (truncated, original)  | −2.4   | clamp +20° | +7000           |
| Broken deck (v1)               | −2.4   | clamp +20° | +10500          |
| **Fixed deck (v2)**            | **−1.6** | clamp +20° | **+5800**       |

- Mean sinkage shifted from −2.4 to −1.6 with the fixed deck —
  the deck adds buoyancy as the hull bobs, reducing the net
  sink-into-water displacement.
- Pitch still saturates at the clamp. The fixed deck reduces M_y
  by ~30 %, but it's still strong enough to overcome the metacentric
  restoring + my linear damping at the simple-Euler integration
  scheme.

## Why pitch is still saturated (open issue)

Linear analysis: K_pitch = ρ_w · g · B · L³ / 20 = 82 000 (parabolic
waterplane formula). With M_y_wave ≈ 5800 (from the BDIM
measurement), small-angle equilibrium θ_eq = M_y_wave / K_pitch ≈
0.07 rad = 4°. The simulation should settle around this value.

It doesn't, because:

1. **Added-mass inertia** is missing. The slender-body `I = M·L²/12`
   underestimates the actual pitch inertia by ~30–50 %, but worse,
   the *frequency* of the natural pitch mode is what the explicit
   Euler integrator's stability depends on. With under-estimated
   inertia, the period is too short, and the integrator overshoots.
2. **Linear damping coefficient** is a fudge — set as 0.2/dt
   (decay time ~5 steps). Adequate for one cycle but the wave-
   moment forcing is sustained, not transient.
3. **Implicit integration** would absorb both issues; explicit Euler
   on stiff systems is fragile.

## Next steps (not implemented)

1. Switch to a Newmark-β or generalized-α integrator for the rigid-
   body EOM. ~50 lines of code, would fix the explicit-Euler
   stability issue.
2. Add added-mass correction: `I_pitch_eff = I_pitch · (1 + k_22)`
   where `k_22 ≈ 0.4` is the slender-body added-mass coefficient.
3. Reduce the wave-induced M_y at θ = 0 by using the Containership
   instead of Wigley (vertical walls have a different pressure
   distribution near the stern); previous K2 work suggests Cb > 0.7
   hulls are more stable in pitch.

## See also

- `RESULTS-heave-pitch-2dof.md` (K1) — original 2-DOF result.
- `RESULTS-heave-1dof.md` (J1) — stable 1-DOF heave.
- `ShipShapes.jl/src/ShipShapes.jl` — `wigley_sdf` with deck.
- `ShipShapes.jl` commit `8ece26c`.
