# Headline demo (F3) — integrated VLM stack, all features on

Single-purpose script that runs the full stack with every recent
upgrade enabled and dumps a 2-panel figure per frame. This is the
"show me what this thing does" deliverable.

Driver: [`scripts/headline_demo.jl`](scripts/headline_demo.jl).

## What's enabled

| Feature                          | Source                                                  |
|----------------------------------|---------------------------------------------------------|
| Containership hull               | `ShipShapes.Containership` (par_frac=0.5, Cb≈0.75)      |
| WALE LES sub-grid model          | `Turbulence.WALE`                                       |
| Sectional blade smear            | `LiftingSurfaces.smear_blades!` (3 blades × 4 sections) |
| Two-way coupling to rudder       | `LiftingSurfaces.trilinear_inflow(flow.u)`              |
| BladedRotor VLM thrust + torque  | `LiftingSurfaces.BladedRotor` at J=0.32                 |
| Rudder VLM CL/CD                 | `LiftingSurfaces.Rudder` at δ=10° fixed                 |

The rudder sits at `rud_xc = prop_xc + 0.5·R` — inside the rotor race,
the F2 amplification regime.

## Run

```bash
julia ShipFlow.jl/scripts/headline_demo.jl
# Default: 128×64×32, 80 frames + 15 burn-in, ≈ 100s wall
```

Override knobs (env): `WL_NX`, `WL_NY`, `WL_NZ`, `WL_NFRAMES`,
`WL_BURNIN`, `WL_FR`, `WL_RE`, `WL_DELTA`, `WL_PAR_FRAC`.

## Output

80 PNG frames in `ShipFlow.jl/runs/headline_demo/frames/frame_*.png`.
Each frame has:

- **Top panel**: η elevation heatmap (free surface) in plan view,
  with hull silhouette, rotor (×), rudder (▼).
- **Bottom panel**: centreline u_x velocity slice (side view), with
  nominal waterline (cyan dashed) and hull cross-section outline.

Assemble into a GIF with:

```bash
ffmpeg -framerate 24 -i runs/headline_demo/frames/frame_%05d.png \
  -filter_complex "[0:v]palettegen=stats_mode=full[p];[0:v][p]paletteuse=dither=bayer" \
  runs/headline_demo/headline.gif
```

## Numbers from the run

- Burn-in: 15 steps in 23.4s (≈ 1.6 s/step including WALE update).
- Rendering: 80 frames in 81.0s (≈ 1.0 s/frame including PNG I/O).
- Max |u| settles around 3.4–3.6 (the rotor race peak), well-bounded.
- Total wall: 104 s (1.7 min) on the default 0.26M-cell grid.

## See also

- `RESULTS-twoway-amplification.md` — F2 quantitative basis for placing
  the rudder inside the race.
- `RESULTS-selfprop-containership.md` — F1 quantitative basis for the
  Containership choice.
- `RESULTS-bladed-vs-swirl.md` — G4 quantitative basis for sectional
  smear (column 3).
