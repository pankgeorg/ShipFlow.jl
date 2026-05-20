# HARNESS — validation strategy across all five workloads

ShipFlow.jl is the one place that knows how to run OpenFOAM, parse its
output, and compare to WaterLily. Each downstream package's
`test/openfoam/` directory calls into ShipFlow.jl's harness functions —
*not* directly into Docker or `foamRun`. That keeps OpenFOAM
orchestration in exactly one repo.

## Two test layers

| Layer       | Frequency      | Where                              | Budget   | Truth source                 |
|-------------|---------------|------------------------------------|----------|------------------------------|
| **L1 unit** | every PR      | each package's `test/runtests.jl`  | < 5 min  | analytic / DNS / canonical   |
| **L2/3 OpenFOAM** | nightly + manual dispatch | `ShipFlow.jl/test/integration/`   | < 4 h    | OpenFOAM tutorial outputs    |

Reasoning: OpenFOAM cases take minutes to hours; running them on every PR
burns CI time. L1 tests catch 90% of regressions in seconds. L2/3 catches
the rest overnight.

## Harness API

```julia
using ShipFlow.Harness

# Run an OpenFOAM tutorial in a pinned Docker container.
# Caches result by (case path, OpenFOAM version, input hash).
of = Harness.run_tutorial(
    "incompressibleFluid/channel395";
    config_overrides = Dict("turbulenceProperties/RAS/RASModel" => "kOmegaSST"),
)
of.U                          # volVectorField, indexed by cell
of.sample_line(p1, p2; n=64)  # linear sample, returns (s, u, p, ...)

# Run a WaterLily simulation with matched setup
wl = Harness.run_waterlily(
    sim_factory = (; kw...) -> ShipFlow.channel395_sim(; kw...),
)

# Compare on a sampling line
report = Harness.compare(of, wl;
    metric = :u_plus_vs_y_plus,
    tolerance = 0.05,
)
@assert report.passed
```

## OpenFOAM-in-Docker

- **Image:** `opencfd/openfoam-default:2406` (free, public, official).
  Pin by SHA in `Harness.OF_IMAGE_SHA` and update via PR.
- **Fallback:** build from the cloned `OpenFOAM/` tree using its
  `Allwmake`. Slow (~30 min) but lets us track `OpenFOAM-dev` if upstream
  diverges from the OpenCFD release.
- **Mounting:** the harness creates a tmpdir, copies the tutorial in,
  applies `config_overrides`, mounts read-write, runs `Allrun` (or the
  case's `system/controlDict` directly), then reads outputs back into
  Julia.
- **Output extraction:** convert with `foamToVTK -ascii` then read with
  `ReadVTK.jl` (already a WaterLily extension dep). For probe/sample data,
  parse OpenFOAM's columnar ASCII dumps directly.

## Reference-data repo

`pankgeorg/cerulean-reference-data` (created lazily — only when something
actually needs to go in it). Git LFS for fields, plain text for
1D/2D sampled data.

Contents:

- Moser–Kim–Mansour DNS data (channel395, channel590, channel950).
- Lee–Moser high-Re channel DNS data.
- Hysing rising-bubble reference fields.
- Martin–Moyce 1952 dam-break experimental data (CSV).
- KCS, DTC, DTMB 5415 published-resistance datasets.
- Cached OpenFOAM tutorial outputs at pinned image SHA (regenerated
  yearly).

Schema: one top-level directory per source. A `MANIFEST.toml` lists hash,
source, citation, redistribution-rights status. **License check before
adding anything.**

## Comparison metrics — the canonical set

Implemented once in `ShipFlow.Harness.Metrics`. Tests pick what they need.

| Metric                              | What it catches                                |
|-------------------------------------|------------------------------------------------|
| `u_plus_vs_y_plus(profile)`         | wall layer; turbulence model bug               |
| `centerline_velocity(line)`         | bulk transport; convection scheme bug          |
| `reattachment_length(box)`          | separation/reattachment; pressure-velocity coupling |
| `front_position(α, t)`              | VoF advection accuracy                         |
| `mass_drift(α)`                     | VoF conservation                               |
| `wave_elevation_along_hull(z, x)`   | free-surface fidelity                          |
| `thrust_coefficient(K_T, J)`        | propeller body-force calibration               |
| `axial_induction(disk_plane)`       | actuator-disk wake                             |
| `resistance_coefficient(F_x, U, S)` | hull drag                                      |
| `wake_fraction(U_disk, U_∞)`        | hull-propeller interaction                     |

Each takes two field/profile pairs (OF and WL) plus a tolerance, returns
a `ComparisonReport` with `passed::Bool`, `error_norm`, and a plot.

## Performance tracking

Per workload, a `bench/` directory with `BenchmarkTools` cases. CI posts
a comparison comment per PR: baseline = `main` of the same repo. Hard
fail if any case regresses by >20%.

Cross-package bench: ShipFlow.jl owns a `bench/end_to_end.jl` script that
runs the headline DTC simulation and reports total wall time + per-step
breakdown. Tracked nightly.

GPU benchmarks: opt-in via a self-hosted runner. Not in the
required-pass set until we have a runner.

## CI infrastructure

- GitHub Actions for L1 (each repo has its own workflow).
- A single `nightly.yml` workflow in ShipFlow.jl that:
  1. Pulls every sibling repo at `main`.
  2. Boots the OpenFOAM container.
  3. Runs L2 + L3 tests for each package in dependency order.
  4. Posts a Markdown summary to a tracking issue in ShipFlow.jl.
- A `dispatch.yml` workflow that lets any sibling repo trigger the
  ShipFlow.jl nightly suite from a PR comment (`/integration-test`).

Total expected cost on free GitHub minutes: < 50 GB-h/month at the
nightly cadence. If we exceed that, switch to a self-hosted runner.

## When to skip the harness

Three cases where running the full harness is wasteful:

1. **Pure documentation changes.** Skip everything except `cargo doc`
   equivalent (`julia --project=. -e 'using Documenter; ...'`).
2. **Geometry-only changes in ShipShapes.jl.** L1 suffices; no fluid solver.
3. **Tagging a release.** Tag from the last green nightly; don't re-run.
