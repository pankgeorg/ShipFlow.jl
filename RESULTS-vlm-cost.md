# VLM solve cost per step

Benchmark of `LiftingSurfaces.rudder_forces` and
`LiftingSurfaces.rotor_forces` per call. Median over 200 calls (50 for
the high-res rotor). All on a single Julia 1.12 thread.

Driver:
[`scratch/vortex_lattice/benchmark_solve_cost.jl`](../scratch/vortex_lattice/benchmark_solve_cost.jl).

## Numbers

| Configuration                  | Panels | Time / call (ms) |
|--------------------------------|-------:|-----------------:|
| Rudder, 16 spanwise × 8 chordwise   |    128 |              5.6 |
| BladedRotor, 3 blades, 12 × 4       |    144 |              9.5 |
| BladedRotor, 3 blades, 24 × 8       |    576 |              149 |

VLM cost scales as ~N²-ish (dense AIC matrix solve), so the
high-res rotor jumps disproportionately when panel count grows from
144 → 576 (4×).

## Context

A typical WaterLily `mom_step!` at 192 × 96 × 48 with VoF/MULES
takes ~1500 ms. The low-res VLM solves (rudder + 12 × 4 rotor =
15 ms) are **about 1 %** of the step cost. The high-res
24 × 8 rotor is ~10 % — still tractable, especially since the
WaterLily step dominates and is independent of panel count.

## Conclusion

For ship-CFD coupling, **the cheap LiftingSurfaces tier costs
essentially nothing on top of the Eulerian solve.** The case for
using `BladedRotor` over `SwirlingDisk` is *almost free* —
the only question is whether the higher physics fidelity is
warranted for your application.

When to use what:

- **`ActuatorDisk`** — design-point sanity checks, drag estimation
  where blade detail doesn't matter. Microseconds per step.
- **`SwirlingDisk`** — same plus you want a swirling wake
  (cavitation-rough comparison, time-averaged wake studies).
- **`BladedRotor`** (this work) — radial loading distribution
  matters, tip vortices, or unsteady blade-passage frequency
  needed. Cost: 5-150 ms depending on panel resolution.
- **`Rudder`** (this work) — any time the rudder side-force matters
  (maneuvering, course-keeping, free-running simulations). Cost:
  ~5 ms.
