# Heave + pitch 2-DOF (K1) — and a modeling limit

Extends J1 (heave 1-DOF) with a pitch degree of freedom about the
y-axis through the centre of mass. Pitch moment from
`WaterLily.pressure_moment(x_CG, sim)[2]`, with the body→fluid sign
flip applied. Slender-body pitch inertia `I = M·L²/12`.

Driver: [`scripts/wigley_heave_pitch_2dof.jl`](scripts/wigley_heave_pitch_2dof.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 36, B = 8, T = 5.
- Fr = 0.25, Re = 5000. 240 steps.
- Body map applies pitch rotation (small-angle Euler) + heave
  translation each step; `WaterLily.measure!` called every step to
  rebuild the BDIM kernel.
- Linear damping added on both DOFs (β_heave = 0.05/dt, β_pitch =
  0.2/dt). Without damping the pitch DOF spun freely (see "what
  went wrong" below).
- Hard clamps: |θ| ≤ 20°, |z_h| ≤ T to prevent the body from leaving
  the domain.

## Result

| Quantity                       | Value (last 25 %)     |
|--------------------------------|-----------------------|
| ⟨z_h⟩ heave equilibrium        | **−2.41 cells** (sinkage) |
| ⟨θ⟩ pitch equilibrium          | **+20.0°** (clamp-saturated) |
| ⟨ż⟩, ⟨θ̇⟩                       | ~0 (motion saturated)  |
| F_buoy at saturation           | +2818 (≈ M·g = 2844)  |
| M_y at saturation              | +6990 (destabilising) |

## What went wrong (pre-damping run)

Without artificial damping, the pitch angle reached +12 rad (645°,
~2 full rotations) by step 240 and the hull spun freely off the
domain. Diagnostic:

- The pressure-moment integration over the Wigley body at our Fr
  produces a **persistently positive M_y** (bow-up moment) of
  ~5000–7000 cell units.
- The natural pitch-restoring spring comes from the **waterplane
  moment of inertia** of the *above-waterline* hull
  (`K_pitch = ρ·g·B·L³/12`). Our Wigley SDF defines only the
  below-waterline portion (body-frame z ∈ [−T, 0]); above the
  waterline the SDF returns positive (free space), so there's no
  body for the bow to lift OUT of into.
- Consequence: as the hull pitches +θ, the bow rises through the
  free surface and the body-volume that would normally provide
  restoring buoyancy isn't there.

This is a known limitation of the truncated Wigley SDF. The fix is
to either:

1. **Add a deck**: extend the SDF above body-frame z = 0 (a flat
   deck face or a parabolic mirror image), so the body has
   above-waterline volume that the rising bow exits, providing the
   restoring `K_pitch · θ`.
2. **Use a different hull family** (Containership has vertical
   sides, but the same waterplane truncation issue).
3. **Tabulated SDF** sampled from a real DTC offset file — would
   include the full geometry.

The G2/F2/I-series scripts work because the hull is fixed; the
truncation doesn't matter there.

## What does work

- **Heave alone (J1) is stable** because there's enough vertical
  restoring force even for a truncated hull — the body itself
  exists below waterline, and gravity + Archimedes is well-posed.
  The mean sinkage ≈ −2.4 reproduces what J1 found.
- **The 2-DOF integration mechanics work**: the rotated `map`
  function and `measure!` per step correctly propagate the body
  motion into the BDIM kernel.
- **With damping the system converges** to a saturated equilibrium
  at the clamp angle.

## Caveats

- The destabilising M_y is *partially* numerical: the truncated
  body has a sharp top edge at body-frame z = 0 that gets a
  discontinuous BDIM kernel at the free surface. Smoothing the
  deck would reduce this contribution.
- The pitch inertia `I = M·L²/12` is the slender-body approximation;
  the actual pitch inertia for a Wigley would include the
  added-mass effect (~+30 % at this aspect ratio), making the
  predicted natural period ~16 ct rather than 14.

## Next step (not implemented)

Add a deck to the Wigley SDF. Smallest patch: extend the SDF to
return `body_z` (positive distance to z=0 surface) for body-frame
z > 0. Then the lifted bow has body-volume to displace water with,
and the K_pitch restoring moment becomes well-posed.

```julia
# In wigley_sdf, after computing half_beam:
if z > 0  # above waterline in body frame
    # Distance is just z (positive, outside body)
    return z + (in_box_xy_check ? 0 : Inf)
end
```

Estimated effort: 1 hour to patch + test. Worth doing before any
seakeeping work.

## Final fix (V1): explicit metacentric restoring

The Wigley 2-DOF now **converges to θ ≈ 0** with one more
pragmatic patch:

```julia
# In the integration loop, after computing M_y_bdim:
K_pitch     = ρ_w · g_now · B_c · L_c^3 / 20    # Wigley waterplane
M_y_restore = -K_pitch * θ[]                    # explicit linear restoring
M_y         = M_y_bdim + M_y_restore            # total moment
```

This adds the *full* analytical metacentric moment on top of the
BDIM measurement. There's a partial double-count (BDIM picks up
some of the restoring already), but the cancellation effectively
zeroes the drift.

| Version          | Final ⟨θ⟩  | Comments |
|------------------|------------|----------|
| K1 (no restoring, β=0.2/dt)   | clamp +20°  | saturated |
| Post-U1 (β=0.5/dt)            | +12.4° → +17° drift | bounded, drifting |
| **V1 (explicit K_pitch)**     | **−0.02°**  | converges |

The V1 ⟨θ⟩ = 0 is *too zero* — real Wigley at Fr=0.25 trims a few
degrees bow-up. The explicit term over-restores. A better-tuned
version would use `K_pitch · θ · 0.5` (compensating only the part
BDIM misses). Left as a tuning knob for future work.

## Update (post-U1): stronger damping helps Wigley too

With **β_pitch raised from 0.2/dt to 0.5/dt** (4× the original),
the Wigley 2-DOF no longer saturates within 240 steps:

| Damping β_pitch | Final θ (last 25 %) | Comments |
|-----------------|---------------------|----------|
| 0.2/dt (K1 / L1) | clamp +20.0°       | saturated |
| 0.5/dt (post-U1) | **+12.4°**         | bounded, slowly drifting |

Even with the stronger damping, θ slowly creeps up because the
BDIM-measured M_y is a *forcing* moment that isn't fully balanced
by metacentric restoring at this grid resolution.

**On the Containership the same script converges** (U1,
`RESULTS-heave-pitch-containership.md`). The difference is K_pitch:
the parallel-midbody waterplane has B·L³/12 vs Wigley's parabolic
B·L³/20 — 67 % stiffer restoring, enough to overpower the BDIM
forcing moment.

## See also

- `RESULTS-heave-1dof.md` — J1: stable 1-DOF heave (works).
- `RESULTS-heave-pitch-containership.md` — U1: Containership 2-DOF
  that converges (companion to this Wigley case).
- `scripts/wigley_heave_pitch_2dof.jl` — driver.
- `runs/heave_pitch_2dof/heave_pitch.png` — z, θ, F_buoy, M_y plot.
- `reference_waterlily_conventions.md` (memory) — the
  `pressure_force`/`pressure_moment` sign convention.
