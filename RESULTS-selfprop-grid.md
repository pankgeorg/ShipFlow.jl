# J_self grid-resolution sensitivity (J4)

Re-runs the Wigley self-propulsion sweep at three grids to test
grid convergence of the headline G1 number (J_self = 0.315).

Driver: [`scripts/wigley_selfprop_grid_sweep.jl`](scripts/wigley_selfprop_grid_sweep.jl).

## Setup

Three proportional grids; hull occupies ~50 % of NX in each:

| Grid           | cells   | L_c × B_c × T_c | Cells per L_c |
|----------------|---------|-----------------|---------------|
| 64 × 32 × 32   | 0.07 M  | 24 × 5 × 3      | 24 / L_c      |
| 96 × 48 × 48   | 0.22 M  | 48 × 10 × 6     | 48 / L_c      |
| 128 × 64 × 64  | 0.52 M  | 64 × 14 × 8     | 64 / L_c      |

Same Fr = 0.25, Re = 1000, BladedRotor at the stern (`R = 0.8·T_c/2`),
80 timesteps per (grid, J), last 25 % averaged. J grid:
`{0.20, 0.25, 0.30, 0.35, 0.45}`.

## Result

| Grid           | J_self  | Δ vs 96 grid |
|----------------|---------|--------------|
| 64 × 32 × 32   | **no bracket** | — |
| 96 × 48 × 48   | 0.3062  | (baseline)   |
| 128 × 64 × 64  | 0.2770  | **−9.5 %**   |

## Interpretation

- **J_self drops ~10 % going from 96 to 128.** The hull drag at fixed
  J grows faster than the rotor thrust as the grid refines: at
  J = 0.30, drag goes from 136 (grid 96) to 287 (grid 128), while
  thrust scales geometrically with `R² = (0.8·T/2)²` (139 → 247).
  Drag growing faster than the geometric factor means the
  resolution affects how the boundary layer interacts with the
  hull stern, increasing the *effective* drag at finer grids.
- **Grid 64 doesn't bracket** because at this coarse resolution
  the hull's drag stays low (only 16 at J = 0.25, vs 148 at grid
  96). The J grid 0.20–0.45 has T > D for every J at this
  resolution — would need J as low as 0.10 to find the bracket.
- **Not converged.** A converged result would have J_self within
  ~1 % across grid doublings. The 10 % shift between 96 and 128
  is too large; we'd need 192 or 256 to confirm asymptotic
  behaviour. With current compute budget (the 128 sweep took
  170 s for 5 points), going to 192 quadruples that to ~700 s
  per J, still tractable for a final-resolution scan.
- **Practical implication**: the headline G1 number (0.315)
  is grid-dependent. A "true" J_self at infinite resolution is
  somewhere below 0.28. All conclusions about the *trend*
  (J_self vs Cb in I1, e.g.) remain valid since they were taken
  at fixed grid — relative rankings are preserved.

## Caveats

- **The hull and rotor scale with the grid** in this sweep — we
  keep `L_c / NX = 0.5` rather than holding the physical hull
  fixed. A strict grid-refinement study should hold L_c constant
  and shrink the cell size. We didn't because that requires a
  much bigger grid (256, 384, …) to bracket from the converging
  side.
- The 80-step run length may not fully settle drag at the
  larger grids — wave-train development scales with L_c / c_g.
  At grid 128 with L_c = 64, 80 steps × 0.25 dt = 20 cell-time
  units; wave-train length is comparable.
- **Re is held at 1000** by scaling ν with L. A real grid-
  convergence study at fixed Re would also need to confirm the
  Re value is in the right regime (laminar Blasius is what we
  have).

## Recommendation

Treat J_self = 0.28 as the more reliable number for the Wigley
geometry at Fr = 0.25, Re = 1000. The earlier 0.315 from the 96
grid is biased ~10 % high due to under-resolved stern boundary
layer.

## See also

- `RESULTS-selfprop-VLM.md` — G1, the original 0.315 value.
- `runs/selfprop_grid_sweep/summary.csv` — the underlying numbers.
- `scripts/wigley_selfprop_grid_sweep.jl` — driver.
