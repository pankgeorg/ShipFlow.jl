# Heave 1-DOF on the Containership hull (R1)

Validates the J1 heave mechanics on a second hull family. Result:
the Containership is **less BDIM-friendly** than the Wigley — the
measured Archimedes is ~37 % of analytical at static, not 88 %
(per L3 hydrostatic check on Wigley).

Driver: [`scripts/containership_heave_1dof.jl`](scripts/containership_heave_1dof.jl).

## Setup

Identical to J1 except `ShipShapes.Containership` (par_frac=0.5,
Cb=0.75) instead of `Wigley`. M_ship = ρ_w · V₀ where
V₀ = 0.75 · L · B · T = 1080 cell³ (vs Wigley's 640).

## Result

The hull sank monotonically. At step 60 it had already gone below
z_h = −19 cells (well below the original keel at z_h = −5), and the
BDIM kernel reports F_hyz = 0 there (the body has left the
domain).

Diagnostic from the early steps:
- M·g = 10 · 0.444 · 1080 = **4800** (analytical weight)
- F_hyz at step 1–5 = **1700 – 1800** (measured)
- Ratio = 37 %

For comparison, the Wigley hydrostatic check (L3) gave F_hyz /
(ρ·g·V) ≈ 88 %. The Containership's bias is *three times worse*.

## Why

- **Sharp corners**: the Containership SDF has C0-continuous but not
  C1-continuous corners at the bow/stern parallel-taper joins (at
  `|2x/L| = par_frac`). The BDIM kernel μ₀ has finite width
  (~1 cell), so it smears these corners over a region of the
  pressure field. The integrated pressure × normal × μ₀ then
  under-counts the corner contributions.
- **Vertical sides + flat keel**: the Containership has a flat
  bottom (z = −T constant for `|2x/L| ≤ par_frac`). The keel
  surface normal is straight downward, and the BDIM kernel
  there shows the same finite-width smearing — the result is
  body-shape-dependent magnitude of the bias.
- **Higher Cb amplifies**: the Containership has 1.7× the Wigley's
  volume in the same domain, so the absolute deficit scales
  proportionally. The Wigley's 12 % deficit at 640 cell³ becomes
  ~21 % in the same units; on top of that, the corner contribution
  adds another ~40 %.

## What this means for the project

- **Heave/pitch 6-DOF on the Containership is currently broken.**
  Without a corrected M_effective (or a smoothed SDF), the hull
  doesn't float. This is a hull-specific calibration step.
- **Practical fix**: pre-run a hydrostatic-sanity test like L3 for
  each new hull, measure the BDIM bias, and use
  `M_effective = (F_hyz_measured / (ρ·g·V_analytical)) · M_analytical`
  for the floating equilibrium. The Wigley needs M·0.88, the
  Containership needs M·0.37. Plug into the heave script.
- **Long-term fix**: smooth the Containership SDF with a small
  rounding radius at the bow/stern joins. ~10 lines of code change
  to `containership_sdf`. Reduces the corner-smearing artefact.

## Caveats

- The 240-step run shows the hull leaving the domain, then F_hyz =
  0 because BDIM finds no body. After that point ż_h diverges
  linearly. The "tail-25 %" mean is meaningless.
- L3's 88 % bias on Wigley may also be smaller than expected if we
  had a finer grid or larger hull-to-cell ratio. The Containership's
  37 % is at the same grid; refining would help both.

## Follow-up: bias calibration alone is insufficient

Tried setting `M_ship = 0.37 · M_analytical` (i.e. M matches the
BDIM-measured Archimedes at static). Result: hull still diverges
within 90 steps. The bias factor at z_h ≠ 0 is highly nonlinear
because the Containership SDF has a *sharp* top edge at body z = 0
— as the hull moves, the BDIM kernel jaggily picks up and releases
volume.

## Fix (T1): add a deck to the Containership SDF

Extended `containership_sdf` with a `deck_h` parameter (committed in
ShipShapes.jl alongside the L1 Wigley deck). The deck interior uses
the true minimum-distance SDF (`−min(d_side, d_top, d_kel)`), same
form as the corrected Wigley deck. ShipShapes test count: 42 (no
regressions).

**Result with deck_h = T/2 + strong damping (β=0.5/dt):**

| Quantity                         | Value     |
|----------------------------------|-----------|
| ⟨z_h⟩ (last 25 % avg)            | **−0.21 cells** |
| ż_rms (last 25 %)                | **0.002** (settled) |
| F_hyz at equilibrium             | **4805**  |
| M·g                              | 4800      |
| **Bias factor with deck**        | **~1.00** (was 0.37 without deck) |
| Residual F_net                   | 0.1 % of weight |

**The deck SDF restores the analytical-Archimedes match.** The
Containership now floats at the expected equilibrium z_h ≈ 0 with
the correct sinkage signature. Strong damping is needed because
the natural heave frequency is high enough that explicit Euler
overshoots; a Newmark-β integrator would replace this with
unconditional stability.

The session's takeaway: **for any hull with corners or sharp top
edges, the BDIM bias is large enough to break 6-DOF integration
unless the SDF includes an above-waterline body**. L1 fixed the
Wigley; T1 fixed the Containership.

## See also

- `RESULTS-heave-1dof.md` (J1) — Wigley counterpart, stable.
- `RESULTS-hydrostatic-check.md` (L3) — Wigley BDIM bias (12 % off).
- `RESULTS-wigley-deck-sdf.md` (L1) — Wigley deck fix that this
  Containership case also needs.
- `runs/heave_1dof_containership/heave.png` — divergent z_h trajectory.
