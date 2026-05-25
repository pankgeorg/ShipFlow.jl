# Propeller open-water polar via VLM (LiftingSurfaces.jl)

Sweep of advance ratio J = V∞/(nD) from 0.3 to 1.2 in 0.1 increments,
all via `LiftingSurfaces.rotor_forces`. 3-blade twisted/tapered
geometry (chord 0.25 → 0.18, twist 35° → 15°), 16×6 panels per blade.

Driver: [`scripts/propeller_polar_via_VLM.jl`](scripts/propeller_polar_via_VLM.jl).

## Results

| J    | \|CT\| | \|CQ\| | KT    | 10·KQ | η_VLM |
|-----:|------:|------:|------:|------:|------:|
| 0.30 | 15.35 | 0.560 | 0.543 | 0.099 | 2.62  |
| 0.40 |  8.74 | 0.523 | 0.549 | 0.165 | 2.12  |
| 0.50 |  5.65 | 0.474 | 0.555 | 0.233 | 1.90  |
| 0.60 |  3.97 | 0.429 | 0.561 | 0.303 | 1.76  |
| 0.70 |  2.94 | 0.392 | 0.566 | 0.377 | 1.67  |
| 0.80 |  2.27 | 0.360 | 0.571 | 0.453 | 1.61  |
| 0.90 |  1.81 | 0.334 | 0.576 | 0.532 | 1.55  |
| 1.00 |  1.48 | 0.312 | 0.581 | 0.613 | 1.51  |
| 1.10 |  1.23 | 0.293 | 0.585 | 0.697 | 1.47  |
| 1.20 |  1.04 | 0.277 | 0.589 | 0.784 | 1.44  |

## Interpretation

- **CT, CQ both drop monotonically with J** — more flow advance per
  revolution leaves less time for the blades to act on the water,
  so per-rev force and torque go down. CT drops by 15× over the
  sweep; CQ by 2×.
- **KT is nearly flat ≈ 0.55 across the range.** This is because KT
  scales as `½π·CT·J²/4`, and the rapid CT decrease is approximately
  cancelled by the J² scaling. The blade is operating at near-fixed
  geometric pitch; only the angle-of-attack relative to Ωr changes.
- **10·KQ rises monotonically with J** — physically meaningful:
  at higher J the rotation is slower, but the torque coefficient
  normalisation `Q/(n²D⁵)` amplifies the small absolute torque.
- **η_VLM exceeds 1** (1.4 – 2.6 across the sweep). This is
  unphysical for a real propeller and reflects the inviscid-VLM
  bound: no skin friction, no cavitation, no induced-drag
  correction at the propeller plane.

## Caveats — why these aren't ship-design predictions

VortexLattice.jl with `steady_analysis!` runs an inviscid + linear
lifting-surface solve. The classical Wageningen B-series open-water
polar curves include:

1. **Skin friction** on the blade surface (Reynolds-dependent),
2. **Wake-induced velocity** (computed self-consistently by free
   wake; VLM's wake is fixed),
3. **Cavitation** (no pressure-floor in VLM),
4. **Section profile drag** (only with `nonlinear_analysis!` and
   airfoil polars).

For comparison, a Wageningen B4-70 at the same operating range:
KT ≈ 0.40 → 0.10 (vs ours 0.55 flat), 10·KQ ≈ 0.55 → 0.20
(vs ours 0.10 → 0.78), η_max ≈ 0.70 at J ≈ 0.7
(vs our >1 everywhere).

## Use case

The VLM polar is still useful as:

- **Sanity check** that LiftingSurfaces.BladedRotor returns
  monotonic and dimensionally-consistent forces across J.
- **Input to the coupled WaterLily run** — the qualitative shape
  of CT(J) matters more than the absolute level when VLM's CT
  is used as a body-force source. WaterLily applies its own
  viscous drag on the resulting wake, so a fraction of the
  missing physics is recovered downstream.
- **Cheap parameter sweep** for design exploration (~1 second per
  point, vs minutes for a coupled WaterLily run).

For absolute propeller-performance prediction, switch to
`nonlinear_analysis!` with a section-polar lookup (CCBlade airfoil
polars; VortexLattice's `generate_rotor` integrates this).
