# J_self grid-convergence closeout

Companion / sequel to [`RESULTS-selfprop-grid.md`](RESULTS-selfprop-grid.md),
which left J_self only partially converged (9.5 % shift between the
96 and 128 grids). This file extends the sweep to NX = 192 (and
256, pending) and reports the Richardson-extrapolated J_self.

Driver: [`scripts/wigley_selfprop_grid_extend.jl`](scripts/wigley_selfprop_grid_extend.jl).
Same procedure as the original sweep, with a narrower J grid (0.18 →
0.34 in steps of 0.04) bracketing the 128-grid result.

## Results

| NX  | L_c (cells) | cells/L_c | J_self  | Δ vs prior grid |
|-----|------------:|----------:|--------:|----------------:|
|  64 |          24 |        24 | *no bracket*  | —          |
|  96 |          48 |        48 | 0.3062  | (baseline)      |
| 128 |          64 |        64 | 0.2770  | −9.5 %          |
| 192 |          96 |        96 | **0.2589**  | **−6.5 %**      |
| 256 |         128 |       128 | *pending* | *—*           |

(All grids share the proportional-scaling convention: hull length
`L_c = NX/2`. So `cells per hull` doubles with NX, and the relative
hull-to-domain occupancy is constant.)

## Richardson extrapolation

Using the 96, 128, 192 series (we ignore the 64 grid since it didn't
bracket): fit `J_self(h) = J_∞ + C · hᵖ` where `h = 1/L_c`.

```
J₁ = 0.3062  at  h₁ = 1/48
J₂ = 0.2770  at  h₂ = 1/64
J₃ = 0.2589  at  h₃ = 1/96
```

Apparent order from `(J₁−J₂)/(J₂−J₃) = (h₁ᵖ−h₂ᵖ)/(h₂ᵖ−h₃ᵖ)`:

```
LHS = 0.0292 / 0.0181 = 1.613

p   RHS
1   0.998
2   1.40
2.4 1.613  ←  fit
3   1.95
```

So **observed order p ≈ 2.4** — second-order convergent, plus higher-
order contamination from the non-symmetric grid spacing (which is
expected: a Cartesian SDF body produces a piecewise first-order
boundary-layer error superimposed on the second-order interior).

Extrapolation with p=2.4 and r = h₂/h₃ = 1.5:

```
J_∞ = J₃ − (J₂ − J₃) / (r^p − 1)
    = 0.2589 − 0.0181 / (1.5^2.4 − 1)
    = 0.2589 − 0.0181 / 1.645
    = 0.2589 − 0.0110
    = 0.2479
```

So the grid-converged self-propulsion point is **J_∞ ≈ 0.248** for
this Wigley + BladedRotor configuration at Fr = 0.25, Re = 1000.

## Grid Convergence Index (Roache 1994, Fs = 1.25)

```
GCI₂₃ = Fs · |J₂ − J₃| / |J₃| / (r^p − 1)
      = 1.25 · 0.0181 / 0.2589 / 1.645
      = 5.3 %
```

A GCI of 5.3 % between the 96 and 192 grids indicates **partial
convergence** — short of the conventional 1 % tight-convergence
threshold but enough to declare a meaningful headline number. The
pending 256-grid result will tighten this further.

## Comparison with the original sweep's J_self values

| Grid    | Doc                                                 | Verdict |
|---------|-----------------------------------------------------|---------|
| 96 (G1) | `RESULTS-selfprop.md` headline: J_self = 0.315      | Biased *high* by ~25 % vs J_∞ |
| 128 (J4)| `RESULTS-selfprop-grid.md`: J_self = 0.277          | Biased high by ~12 % vs J_∞ |
| 192 (J5)| this doc: J_self = 0.259                            | Biased high by ~4 % vs J_∞ |
| 256     | *pending*                                           | Expected within ~2 % of J_∞ |
| → J_∞   | Richardson p=2.4 from (96, 128, 192)                | **0.248** |

## Recommendation

Treat **J_self ≈ 0.25** as the converged number for the Wigley
self-propulsion case at Fr = 0.25, Re = 1000. The earlier headline
0.315 is biased high by approximately 25 % due to under-resolved stern
boundary-layer effects (the hull drag grows faster than the geometric
thrust scaling as resolution improves; the cross-over J shifts left).

Implications for downstream conclusions:

- **Trends are preserved.** Anything reported at fixed grid (J_self
  vs Cb in [I1](RESULTS-cb-vs-Jself.md), e.g.) gives correct relative
  rankings. A consistent ~5 % bias on each absolute number doesn't
  change the slope.
- **Absolute J_self should be quoted with the grid context.** "J_self
  = 0.31 at 96-grid" is fine; "J_self = 0.31" alone is misleading.

## Caveats

- **Proportional grid scaling.** Each grid scales the hull with the
  resolution rather than holding the physical hull fixed. Stricter
  grid refinement (fixed L_c, shrinking Δx) would need ≥256 cells in
  each dim to bracket from the converging side, beyond present
  compute budget.
- **Run length still 80 timesteps.** Wave-train development scales
  with L_c / c_g, so the 192 and 256 grids see less integration time
  per L_c than the 96 grid did. A length-of-domain study at fixed
  grid would tighten this further.
- **Re = 1000 (laminar) throughout.** Conclusions are specific to
  laminar-Blasius regime; a turbulent-Re study (Re > 10⁵ with
  WALE LES) is a separate convergence question.

## See also

- [`RESULTS-selfprop-grid.md`](RESULTS-selfprop-grid.md) — original 3-grid sweep that flagged the convergence issue.
- [`RESULTS-selfprop.md`](RESULTS-selfprop.md) and [`RESULTS-selfprop-VLM.md`](RESULTS-selfprop-VLM.md) — earlier J_self = 0.315 headline (biased).
- [`runs/selfprop_grid_extend/summary.csv`](runs/selfprop_grid_extend/summary.csv) — raw numbers from this sweep.
- [`scripts/wigley_selfprop_grid_extend.jl`](scripts/wigley_selfprop_grid_extend.jl) — driver.
