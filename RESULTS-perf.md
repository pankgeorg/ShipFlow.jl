# RESULTS — performance review

This is a profile of where time goes in the WaterLily ship-flow stack,
plus the wins we've shipped and the larger ones still on the table.
Numbers are from the 96 × 48 × 48 headline grid (~220 k cells) on a
JULIA_NUM_THREADS=80 host, SIMD backend, Float32 fields.

## Baseline throughput

| Configuration                                        | steps/s | sec/step |
|------------------------------------------------------|--------:|---------:|
| Float64, scalar ν                                    | 1.66    | 0.60     |
| Float32, scalar ν                                    | 2.18    | 0.46     |
| Float32, array ν (Turb+VoF coupled)                  | 1.91    | 0.52     |
| Same + every-step force eval (old scan loop)         | 1.16    | 0.86     |
| Same + every-25-step force eval (new scan loop)      | 1.67    | 0.60     |

## Component breakdown (Float32, scalar ν, isolated calls)

Per-call cost averaged over 10 invocations:

| Component                | per call (ms) | calls/step | step share |
|--------------------------|--------------:|-----------:|-----------:|
| `conv_diff!`             |          21   |     2      |    9 %     |
| `accelerate!`            |          37   |     2      |   15 %     |
| `BDIM!`                  |           3   |     2      |    1 %     |
| `mom_project!`           |         491   |     2      |     —      |
| `Poisson solver!`        |         147   |     2      |   60 % (within `mom_project!`) |
| `pressure_force`         |         214   |  0–1       | 30-65 % when called |
| `viscous_force`          |         124   |  0–1       | 18-40 % when called |
| `step_vof!`              |          43   |     1      |    9 %     |
| `update_νt!` (WALE 3D)   |          49   |     1      |   10 %     |

The Poisson solver dominates mom_step!: roughly 60 % of step time even
when only ~3 V-cycles are needed. Forces are not in the hot path of
production runs but were 60 % of scan time before the every-25-step
fix.

## Wins shipped during the loop session

1. **Force subsampling in the self-propulsion scan**
   (`wigley_self_propulsion_scan.jl`, commit `6c19bfc`):
   instead of evaluating `pressure_force` + `viscous_force` every
   single step, sample every `NSTEPS_PER ÷ 25` steps until the
   averaging window. ~35-40 % faster end-to-end on the headline scan.

2. **`pois_tol = 1e-6` instead of `1e-8`**: with 3 V-cycles either
   way, the looser tolerance leaves headroom in case the residual
   structure shifts. No throughput impact in steady cases.

3. **WaterLily Hook 1 ν stored by reference** (commit `e4b8854`):
   *not* a perf change but a correctness one. The previous copy meant
   downstream LES/VoF couplings were inert; physics was wrong even
   though throughput was the same.

4. **Float32 throughout the Kelvin and snapshot drivers**: ~30 % gain
   versus Float64 with no precision loss on the cell-unit fields.

## Cheap wins still on the table

| Item                                                   | est. gain | effort |
|--------------------------------------------------------|----------:|-------:|
| Float32 in the remaining scan scripts                  |     20 %  | tiny   |
| Skip `viscous_force` when ν is uniformly zero          |    ≈ 0    | tiny   |
| Cache `AutoBody.measure` (n, kern) for static bodies   |     50 %  | medium — upstream WaterLily change |
| Pre-allocate `mass_repair` workspace (currently allocs)|      5 %  | tiny   |

## Expensive wins not on the table today

1. **GPU backend (CUDA / AMDGPU)** — 5–10× speedup on the same grid.
   WaterLily supports both via `mem=CuArray`. The downstream packages
   inherit this; the only gotcha is the Hook-1 fix above, which we now
   verify in a test (commit `0eb57c0`).

2. **Better Poisson preconditioner** — geometric multigrid (current)
   stalls at the ρ-jump interface, forcing the tight tol we used.
   Algebraic multigrid, GMRES, or a level-set-aware coarsening would
   each save ~30 % on hard cases.

3. **Adaptive mesh refinement** — most of the action is near the hull
   and at the free surface, but the uniform Cartesian grid pays full
   cost everywhere. AMR would change the cost model entirely; out of
   scope for a hooks-and-validate phase but worth flagging.

## Notes for the user

The "feels slow" perception is real and load-bearing on the scan
loops — 30+ minutes per 7-point scan was the lived experience before
the force-subsampling fix. With the fix, a 7-point scan now sits at
~10-12 min on this host. The single biggest remaining lever for CPU
runs is Float32 hygiene; the biggest available win is GPU.

## MULES overhead

After the MULES landing (commits `5e78ecf`, `13350c1`, `62c8f79`), measure on
the same 96×48×48 headline grid (Float32, scalar ν disabled, no turbulence):

| step component                          | vanLeer | MULES  | factor |
|-----------------------------------------|--------:|-------:|-------:|
| `step_vof!` / `step_vof_mules!` alone   |  27 ms  | 263 ms |  ~10× |
| full mom_step + vof + WALE              | 378 ms  | 608 ms |  +60 % |

MULES adds ~230 ms/step on this grid — significant but acceptable for the
free-surface cases where the Kelvin wedge fidelity matters. Likely wins:

| Item                                                  | est. gain | effort |
|-------------------------------------------------------|----------:|-------:|
| Pre-allocate MULES workspace in `VoFFlow` struct      |     40 %  | small  |
| Fuse Φ-compute + α_UD-update into one pass            |     20 %  | small  |
| Compute `_local_extrema` inline in P_pos/P_neg loop   |     10 %  | small  |

Combined: MULES could realistically come down to ~120 ms (4× rather than 10×
the `step_vof!` baseline). Not done in this session — flagged for follow-up.
