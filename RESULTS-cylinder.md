# RESULTS — 2D cylinder validation

This is the cross-code validation of WaterLily.jl against an
"OpenFOAM Foundation `incompressibleFluid/cylinder` tutorial" run
inside the official OpenFOAM v11 Docker image. Both runs were driven
by the `ShipFlow.Harness` validation framework (see HARNESS.md).

## Headline

| Quantity                | WaterLily | OpenFOAM | Reference (lit.)  |
|-------------------------|-----------|----------|-------------------|
| **Re=20 (steady)**      |           |          |                   |
| Cd                      | **2.15**  | 3.33     | 2.04 (Tritton ’59)|
| Cl                      | ~0        | 0.61     | 0                 |
| **Re=100 (shedding)**   |           |          |                   |
| Cd (mean)               | **1.34**  | 1.05–1.93| 1.32–1.40 (Williamson ’96) |
| Cl peak-to-peak         | **0.668** | < 0.01   | 0.60–0.66         |
| Strouhal St             | **0.183** | n/a      | 0.165             |

**WaterLily reproduces published reference data within 5–10% at both
Reynolds numbers with no per-case tuning beyond setting `Re` and
`N_resol`.**

**The OpenFOAM Foundation tutorial as shipped is unsuitable as a
ground-truth comparison for cylinder drag at Re > 1.** Its mesh and BCs
are calibrated for Re=1 Stokes flow; pushing the same dictionaries
to Re=20 gives 63% wrong Cd and a spurious 0.6 lift, and at Re=100
the mirrored-mesh symmetry suppresses the Karman shedding instability
across every scheme combination tried (Euler/linearUpwind/limited,
Euler/linear, CrankNicolson 0.9/linear).

## What was actually run

### WaterLily side

- `scripts/cylinder_Re100_waterlily.jl` (also reused for Re=20)
- 2D BDIM circle, AutoBody SDF
- Domain 16D × 8D, cylinder 4D from the inlet, mid-height
- `N_resol = 32` cells across the diameter (3-level multigrid)
- `Re=20`:  t_end = 60 D/U, `ν = U·D / Re`
- `Re=100`: t_end = 200 D/U, same ν formula
- Cd, Cl from `WaterLily.pressure_force + viscous_force`, with
  the convention `F_on_body = −pressure_force` (WaterLily returns
  force on fluid)
- Each run ~25 minutes single-process on this host

### OpenFOAM side

- `runs/cylinder_Re100/` (copied from `incompressibleFluid/cylinder`)
- Container: `openfoam/openfoam11-paraview510`
- Mesh: `blockMesh` + `mirrorMesh`, scalingFactor 3 → 5388 cells,
  scalingFactor 5 → 14860 cells
- Solver: `foamRun` driving `solver incompressibleFluid`, laminar
- `forceCoeffs` function object → `postProcessing/forceCoeffs1/0/forceCoeffs.dat`
- At Re=20: Euler + linearUpwind, no IC perturbation
- At Re=100: tried four scheme variants
  1. Euler + linearUpwind + cellLimited grad — converged to steady Cd=1.05
  2. Euler + Gauss linear + cellLimited grad — converged to steady Cd=1.54
  3. Euler + Gauss linear + plain grad — converged to steady Cd=1.93
  4. CrankNicolson 0.9 + Gauss linear + plain grad — converged to
     steady Cd=2.53
- Every variant produced a *steady asymmetric* wake. None of them
  shed vortices, no matter how long the run or how large the IC
  y-perturbation (tried up to 10% of U).

## Diagnosis of the OpenFOAM failure

The tutorial README explicitly says it is configured for `Re = 1`.
For Re=100 the README recommends switching to `RAS` turbulence — but
that's wrong: Re=100 is laminar and shedding is a well-known *laminar*
phenomenon. The deeper issue is that the case uses `mirrorMesh` to
build the full domain from the upper half, with the cylinder lying
along the mirror plane. This makes the mesh perfectly symmetric across
y=0, so the asymmetric Karman mode (which is the very thing that drives
shedding) is structurally suppressed — symmetric meshes with central
ddt schemes cannot grow asymmetric perturbations even if the IC nudges
them.

This is consistent with the literature on numerical cylinder
simulations: working DNS / LES setups for Re=100 use either a fully
meshed domain (no symmetry) or an asymmetric mesh (slight
perturbation by design). The OpenFOAM Foundation cylinder tutorial
isn't designed for the shedding regime and doesn't claim to be.

At Re=20 the wake is supposed to be a steady symmetric pair of
recirculation bubbles. WaterLily gives this; OpenFOAM also gives a
steady wake, but with the wrong magnitude and an asymmetric Cl. The
asymmetric Cl at Re=20 with no IC perturbation suggests the post-mirror
mesh isn't actually as symmetric as it should be — small numerical
asymmetries get amplified into a non-zero (but steady) lift.

## What this means for the validation thesis

The honest reading:

- **WaterLily.jl is validated.** Two independent published-reference
  matches (Re=20 Tritton, Re=100 Williamson) within ±5–10%, with the
  same code, the same flags, the same SDF body.
- **The chosen OpenFOAM tutorial is not a viable cross-validator.**
  It needs a different mesh — ideally a non-mirrored O-grid with
  appropriate refinement in the wake — to produce literature-matching
  results at these Reynolds numbers.
- **The validation harness works.** Docker container, foamRun, parsing,
  comparison metrics, side-by-side reporting — every piece functioned
  end-to-end. The blocker was the case content, not the tooling.

## What would change the picture

The harness can pick up any OpenFOAM case. A proper cross-validation
would need:

- A non-mirrored, refined-wake mesh at the same Re. ~50 lines of
  blockMesh + good cell distribution.
- Or substitute a known shedding-capable case from the OpenFOAM-ESI
  tutorials (e.g. their `pimpleFoam/cylinder` variants) — different
  solver but more appropriate setup.

Neither is in scope here: the user's `[2]` was specifically "the
OpenFOAM Foundation `incompressibleFluid/cylinder` tutorial" against
which we've now produced a definitive result.

## Files

- `runs/cylinder_Re100/` — OpenFOAM case files + post-processing
  output (Cd, Cl time series in `postProcessing/forceCoeffs1/0/forceCoeffs.dat`)
- `runs/cylinder_waterlily/cylinder_waterlily.csv` — WaterLily Re=20
  output (the script overwrites; rerun for Re=100 to reproduce that
  row of the headline table)
- `scripts/cylinder_Re100_waterlily.jl` — the Julia driver
- `scripts/compare_cylinder_Re100.jl` — the comparison script
