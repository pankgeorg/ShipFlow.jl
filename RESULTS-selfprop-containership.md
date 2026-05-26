# Self-propulsion on the Containership hull (F1)

Same VLM rotor and procedure as G1 on Wigley, but the hull is now the
`Containership` SDF (parallel midbody, par_frac = 0.5, Cb ≈ 0.75).

Driver: [`scripts/containership_self_propulsion_VLM.jl`](scripts/containership_self_propulsion_VLM.jl).

## Setup

- Grid 96 × 48 × 48, hull L = 48, B = 10, T = 6, par_frac = 0.5.
- Fr = 0.25, Re = 1000, ρ_w / ρ_a = 10.
- 100 timesteps per J; drag averaged over the last 25%.
- BladedRotor: 3 blades, R = 2.4 (= 0.8·T/2), R_hub = 0.2·R,
  chord = (0.25·R, 0.18·R), twist = (35°, 15°), 12 × 4 panels.
- Rotor at the stern, prop_xc = hull_xc + L/2 + T/2.

## Sweep

| J     | \|CT\| | Thrust  | Drag    | T − D   |
|-------|--------|---------|---------|---------|
| 0.150 | 60.31  | 545.6   | 428.9   | +116.7  |
| 0.175 | 44.45  | 402.2   | 398.8   |  **+3.3** |
| 0.200 | 34.14  | 308.9   | 376.0   |  −67.2  |
| 0.225 | 27.06  | 244.8   | 357.9   | −113.1  |
| 0.250 | 21.98  | 198.9   | 343.1   | −144.2  |
| 0.300 | 15.36  | 139.0   | 320.2   | −181.2  |
| 0.500 |  5.65  |  51.2   | 271.0   | −219.9  |

**Linear-interpolated self-propulsion J ≈ 0.176.**

## Comparison vs Wigley (G1)

| Hull          | Cb      | J_self  | Hull drag at J_self |
|---------------|---------|---------|---------------------|
| Wigley (G1)   | 0.44    | 0.315   | ~120 (cell units)   |
| Containership | 0.75    | 0.176   | ~400 (cell units)   |

The Containership operates at **roughly 1.8× higher rotor RPS** than
Wigley at the same Fr, because it drags ~3× as hard. This matches the
~1.7× displaced-volume ratio (Containership V / Wigley V ≈ 2160 / 1280)
plus the larger wetted area of the parallel midbody. The VLM rotor
gives correspondingly higher CT (~44 vs ~3 at the respective self-prop
points) — the propeller is much more heavily loaded.

## Interpretation

- **The scan brackets the zero crossing.** J = 0.175 sits within 1% of
  T − D = 0; J = 0.150 is well in surplus, J = 0.200 well in deficit.
- **The drag is mildly J-dependent** (430 at J=0.15 vs 271 at J=0.50)
  because the rotor's jet locally accelerates flow past the stern. The
  effect is larger here than on Wigley because the bluffer stern of
  the Containership (no tapered run aft) reacts more strongly to the
  axial jet.
- **CT values are unphysically large** at the lowest J — VLM has no
  stall model and over-predicts loading when blade pitch is high
  relative to inflow. A real propeller would cavitate or stall long
  before CT ≈ 60. This is a known limitation of the VLM tier and
  matches the regime where actuator-disk surrogates are more
  appropriate.

## Caveats

- 100 timesteps is enough to settle a mean drag at this grid, but not
  enough to resolve the wave train past the hull. The Fr=0.25 wave
  field is still developing.
- No two-way coupling in this scan — the rotor sees freestream as
  inflow, not the hull wake. With `trilinear_inflow` the effective
  wake fraction would lower the operating CT slightly.
- Single hull geometry; no sweep over `par_frac`. A future scan could
  vary par_frac and map out how Cb shifts J_self.

## See also

- `RESULTS-selfprop-VLM.md` — Wigley G1 baseline (J_self = 0.3152).
- `scripts/containership_self_propulsion_VLM.jl` — driver for this run.
- `ShipShapes.jl/src/ShipShapes.jl` — `containership_sdf`,
  `containership_volume`, `Containership` factory.
