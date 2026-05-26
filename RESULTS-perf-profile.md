# Per-step allocation profile (J5)

Profiles ShipFlow's actual mom_step!-per-step usage across five
configurations, from bare scalar-ν 3D to the full integrated stack
(VoF + Wigley body + VLM rotor + VLM rudder). Reports allocations and
heap bytes per step on a small (64 × 32 × 32) grid where wall time is
irrelevant compared to allocation counts.

Driver: [`scripts/profile_mom_step.jl`](scripts/profile_mom_step.jl).

## Result

| Case                                              | KiB/step | Allocs/step |
|---------------------------------------------------|----------|-------------|
| A — scalar-ν 3D, no udf, no body                 | 257      | 8,170       |
| B — + VoF (array-ν via `vof.ν`), no udf          | 256      | 8,070       |
| C — + `smear_force!` udf (rotor proxy)           | 257      | 8,074       |
| D — + `smear_force!` + `smear_torque!` udf       | 257      | 8,074       |
| **E — + Wigley body + VLM rotor + VLM rudder**   | **1,008**| **16,150**  |
| `step_vof_mules!` (separate, not in mom_step!)   | **28,656**| **1,192,991** |

## Observations

- **Base 3D `mom_step!` allocates ~256 KiB / 8 k allocs.** Upstream
  WaterLily's published allocation gate (`test/alloctest.jl`) is for
  2D cylinder and reads ~2 KiB. The 3D path costs more because each
  loop direction adds boundary kernels and corner conditions.
- **VoF's array-ν path is essentially free** — the `_ν(ν, I)` and
  `_νf(ν, j, I)` dispatches compile to the same machine code as
  scalar ν, so cases A and B agree to 1 KiB.
- **Smear udfs add ~4 allocs per step** (cases C, D) — essentially
  free; the `smear_force!` and `smear_torque!` kernels are
  alloc-free per the unit tests.
- **VLM is the major mom_step! allocator** — case E hits 1 MiB /
  16 k allocs, **4× case D**. Each call into
  `rudder_forces(rudder, δ, V∞)` builds a fresh `VortexLattice.System`
  from the panel grid (it's not cached). Same for `rotor_forces`
  if called per step.
- **`step_vof_mules!` is the runaway allocator**: 28.6 MiB / 1.2 M
  allocs per call. Eight per-direction loops with `CartesianIndices`
  closures + Zalesak limiter passes — each allocates the inner
  index range eagerly.

## Recommendations (not implemented)

1. **Cache the VortexLattice `System` and `grid`** in the `Rudder`
   and `BladedRotor` structs. The geometry doesn't change between
   steps, only the inflow. Saves ~750 KiB per step in case E.
2. **Hoist `CartesianIndices` allocations** in `step_vof_mules!`
   loops to module-level pre-computed iterators, or rewrite the
   per-direction loops as `@loop ... over ...` macros (which fuse
   well and avoid the inner iterator).
3. **`step_vof_mules!` workspace is already eager** (`_mules_*`
   buffers in `VoFFlow`) — the 28 MiB is *transient* allocation,
   not persistent. A profile run with `@code_warntype` would
   identify the specific hot lines.

## What this is NOT

- Not a wall-time profile. At 0.5 s/step on the 96-grid, allocation
  pressure is modest (~150 MB/s GC pressure for the full stack).
  Wall time is dominated by Poisson solves, not allocation.
- Not a GPU readiness check. The `@loop ... over ...` macro is
  KA-ready; the explicit `for I in CartesianIndices(...)` loops in
  `step_vof_mules!` and in the VLM code are not.

## See also

- `WaterLily/test/alloctest.jl` — upstream baseline.
- `scripts/profile_mom_step.jl` — driver.
