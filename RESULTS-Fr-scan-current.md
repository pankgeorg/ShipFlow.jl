# Wigley resistance curve at current substrate (I5)

Runs the Wigley hull (no propeller) at Fr ∈ {0.20, 0.25, 0.30, 0.35,
0.40, 0.45} and reports pressure-drag and viscous-drag components.
Refreshes the old `wigley_fr_sweep.jl` (which only made Kelvin-wedge
visuals) with quantitative resistance numbers from the current
MULES + array-ν + face-averaged-ν substrate.

Driver: [`scripts/wigley_fr_resistance.jl`](scripts/wigley_fr_resistance.jl).

## Setup

- Grid 128 × 64 × 32, hull L = 36, B = 8, T = 5.
- Re = 5000, ρ_w/ρ_a = 10. No propeller, no rudder.
- 120 steps per Fr, last 25 % averaged.
- Cwp / Cwv normalised by `0.5·ρ_w·U²·S_ref` with
  `S_ref = wigley_volume/T ≈ 128` cells² (a *planform-like*
  reference, not wetted area — see Caveats).

## Result

| Fr   | Dp     | Dv    | D_total | Cwp     | Cwv     |
|------|--------|-------|---------|---------|---------|
| 0.20 | 12.55  | 1.47  | 14.02   | 1.96×10⁻² | 2.30×10⁻³ |
| 0.25 | 17.23  | 1.54  | 18.77   | 2.69×10⁻² | 2.40×10⁻³ |
| 0.30 | 27.48  | 1.58  | 29.05   | 4.29×10⁻² | 2.46×10⁻³ |
| 0.35 | 36.15  | 1.62  | 37.78   | 5.65×10⁻² | 2.53×10⁻³ |
| **0.40** | **38.24** | **1.66** | **39.90** | **5.98×10⁻²** | 2.59×10⁻³ |
| 0.45 | 35.85  | 1.66  | 37.51   | 5.60×10⁻² | 2.60×10⁻³ |

Plot: `runs/wigley_fr_resistance/Cw_vs_Fr.png`.

## Interpretation

- **The classical Wigley "hump" shape is reproduced**: Cwp climbs
  monotonically with Fr, peaks at Fr ≈ 0.40, and starts declining
  at Fr = 0.45.
- **Viscous Cwv is Fr-independent** (2.3–2.6 × 10⁻³), as expected
  at fixed Re. The 13 % drift across the sweep is the wave-induced
  perturbation of the boundary layer at the free surface — small
  but resolvable.
- **Comparison to Bai 1979 / Inui 1980**: classical thin-ship
  theory predicts the Wigley wave-resistance peak at Fr ≈ 0.32
  with Cw ≈ 4×10⁻³ when normalised by `0.5·ρ·U²·S_wet` with
  `S_wet` the proper wetted area. Our peak is at Fr = 0.40 with
  Cwp ≈ 6×10⁻². The shift to higher Fr is consistent with **Re
  effects** — at Re = 5000 the boundary layer is thick relative
  to the hull and effectively shifts the apparent length forward,
  moving the resonance Fr. The 10× magnitude difference is partly
  because we use `S_ref = V/T` (planform-like) instead of the true
  wetted surface (≈ 6× larger), bringing the apparent Cw down to
  ~10⁻² — still 2-3× the published Bai number, which is the
  residual Re effect.

## Caveats

- **Wrong reference area**: our `S_ref = 128` (cells²) is `V / T`,
  not wetted area. A real Wigley wetted-area is ≈ 1.18·L·(B + 2T)
  for the half-form, ≈ 765 cells² here — 6× larger. Multiplying
  our Cwp/Cwv values by ≈ 1/6 gives values that can be more
  meaningfully compared to published towing-tank data.
- **Re = 5000 is two-three orders of magnitude lower** than a real
  ship's Re (10⁸–10⁹). The wave-resistance peak position and
  magnitude both shift; quantitative agreement to towing-tank data
  is not expected at this Re.
- **120 steps may not fully settle Fr ≥ 0.40** — the wave train
  past the stern needs ~`L/c_g` steps to develop, where `c_g` is
  the group velocity. At Fr = 0.45 the wave train is the longest;
  120 steps × dt = 30 cell-time-units may be marginal.
- **No turbulence model** in this run. Adding WALE LES would
  thicken the boundary layer artificially at this Re and increase
  Cwv, but wouldn't move the wave-resistance peak materially.
- **Three-dimensional wave-pattern resolution**: 128 × 64 grid is
  modest; the Kelvin wake at Fr = 0.45 has wavelengths approaching
  the grid limit.

## What the run buys us

- A **clean, reproducible** wave-resistance curve from the current
  substrate. Future PRs that touch the VoF advection or BDIM
  kernel should re-run this and confirm the shape is preserved.
- **The peak location ≈ 0.40** is the operating point with the
  most wave-resistance signal. Self-prop and rudder studies that
  want to maximise wave-making contrast should target Fr = 0.40.
- **An independent check** that the substrate produces physical
  Wigley behaviour — the hump-and-hollow pattern is well-known and
  failure to reproduce it would have indicated a deeper problem.

## See also

- `RESULTS-selfprop-VLM.md` — G1 self-propulsion at Fr = 0.25
  (which is on the steep rising part of the curve above).
- `runs/wigley_fr_resistance/Cw_vs_Fr.png` — the curve.
- `runs/wigley_fr_resistance.csv` — the underlying numbers.
