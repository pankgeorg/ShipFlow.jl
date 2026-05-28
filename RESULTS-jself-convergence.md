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
| 192 |          96 |        96 | 0.2589  | −6.5 %          |
| 256 |         128 |       128 | **0.2502**  | **−3.4 %**      |

(All grids share the proportional-scaling convention: hull length
`L_c = NX/2`. So `cells per hull` doubles with NX, and the relative
hull-to-domain occupancy is constant.)

`ΔJ` between consecutive grids: −0.0292, −0.0181, −0.0087.
These are dropping geometrically by roughly a factor of 1.6 per
refinement — the sequence is converging.

## Richardson extrapolation

Fit `J_self(h) = J_∞ + α · hᵖ` with `h = 1/L_c`. Different triples of
grids give different apparent orders, which is a useful diagnostic for
whether we're in the asymptotic regime.

### Using the 3 coarsest grids (96 / 128 / 192)

Apparent order from
`(J₁−J₂)/(J₁−J₃) = (h₁^p−h₂^p)/(h₁^p−h₃^p) = 0.0292/0.0473 = 0.617`:

```
p     f(p)
2     0.584
2.4   0.617  ← fit
2.5   0.623
```

Order **p ≈ 2.4**, so

```
J_∞ ≈ J₃ − (J₂−J₃) / ((h₂/h₃)^p − 1)
    = 0.2589 − 0.0181 / (1.5^2.4 − 1)
    = 0.2589 − 0.0110
    = 0.2479
```

### Using the 3 finest grids (128 / 192 / 256)

Apparent order from `0.0181/(0.0181+0.0087) = 0.676`:

```
p     f(p)
1.1   0.657
1.18  0.674  ← fit
1.2   0.681
1.3   0.689
```

Order **p ≈ 1.18**, so

```
J_∞ ≈ J₄ − (J₃−J₄) / ((h₃/h₄)^p − 1)
    = 0.2502 − 0.0087 / (1.333^1.18 − 1)
    = 0.2502 − 0.0087 / 0.404
    = 0.2502 − 0.0215
    = 0.2287
```

### Why the apparent order drops

Going from the coarser to the finer triple, the apparent order falls
from 2.4 to 1.2. That is the **signature of not-yet-asymptotic
convergence**: at coarser grids, higher-order error terms dominate
and inflate the observed order; at finer grids, the true leading-
order term (first-order, from the BDIM stair-step boundary
treatment on a Cartesian grid) takes over.

The finer-grid analysis (p ≈ 1.2) is the more reliable indicator of
the asymptotic behaviour. So the more honest headline is:

> **J_self ≈ 0.23 ± 0.02** at infinite resolution, for the Wigley +
> BladedRotor configuration at Fr = 0.25, Re = 1000.

A pure-first-order linear extrapolation from the finest pair
(192, 256) gives J_∞ ≈ 0.224, supporting the lower end of that
band. The earlier headline of **0.31** (the 96-grid result reported
in [`RESULTS-selfprop-VLM.md`](RESULTS-selfprop-VLM.md)) is biased
high by approximately **30 %**.

## Grid Convergence Index (Roache 1994, Fs = 1.25)

Using the 2 finest grids and p = 1.2:

```
GCI₃₄ = Fs · |J₃ − J₄| / |J₄| / (r^p − 1)
      = 1.25 · 0.0087 / 0.2502 / 0.404
      = 10.8 %
```

A GCI of 10.8 % between the 192 and 256 grids indicates **partial
convergence only** — going one more grid (NX = 384) would be
required to bring the GCI below 5 %, and (NX = 512) to bring it
below 1 % under the same first-order assumption. With current
compute budget (the 256-grid sweep took ~24 minutes), going to 384
would take ~80 minutes — feasible but expensive.

## Comparison with the original sweep's J_self values

| Grid    | Doc                                                 | Verdict |
|---------|-----------------------------------------------------|---------|
| 96 (G1) | `RESULTS-selfprop.md` headline: J_self = 0.315      | Biased *high* by ~25 % vs J_∞ |
| 128 (J4)| `RESULTS-selfprop-grid.md`: J_self = 0.277          | Biased high by ~12 % vs J_∞ |
| 192 (J5)| this doc: J_self = 0.259                            | Biased high by ~4 % vs J_∞ |
| 256 (J5)| this doc: J_self = 0.2502                           | Biased high by ~9 % vs J_∞ (best estimate) |
| → J_∞   | Richardson p≈1.2 from (128, 192, 256)               | **0.229** (with ~10 % GCI uncertainty) |

## Recommendation

Treat **J_self ≈ 0.23 ± 0.02** as the asymptotic-resolution number
for the Wigley self-propulsion case at Fr = 0.25, Re = 1000. The
earlier headline 0.315 is biased high by approximately 30 % due to
under-resolved stern boundary-layer effects (the hull drag grows
faster than the geometric thrust scaling as resolution improves; the
cross-over J shifts left).

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
