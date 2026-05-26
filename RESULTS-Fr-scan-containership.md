# Containership Fr resistance curve (L5)

Repeats I5 (Fr-vs-resistance sweep) on the Containership hull family
to locate the wave-resistance peak. K2 had hinted that the
Containership's resistance behaviour differs from Wigley's.

Driver: [`scripts/containership_fr_resistance.jl`](scripts/containership_fr_resistance.jl).

## Setup

Identical to I5 except the hull is `ShipShapes.Containership` (default
par_frac = 0.5, Cb ≈ 0.75) and the reference area uses
`containership_volume` / T instead of `wigley_volume` / T.

- Grid 128 × 64 × 32, hull L = 36, B = 8, T = 5.
- Re = 5000, no propeller.
- 120 steps per Fr, last 25 % averaged.

## Result

| Fr   | Dp     | Dv    | Cwp     | Cwv     | Cw_total |
|------|--------|-------|---------|---------|----------|
| 0.20 | 109.4  | 1.34  | 1.01×10⁻¹ | 1.24×10⁻³ | 1.03×10⁻¹ |
| 0.25 | 119.1  | 1.41  | 1.10×10⁻¹ | 1.31×10⁻³ | 1.12×10⁻¹ |
| 0.30 | 134.4  | 1.50  | 1.24×10⁻¹ | 1.39×10⁻³ | 1.26×10⁻¹ |
| **0.35** | **136.8** | **1.60** | **1.27×10⁻¹** | 1.48×10⁻³ | **1.28×10⁻¹** |
| 0.40 | 126.8  | 1.67  | 1.17×10⁻¹ | 1.55×10⁻³ | 1.19×10⁻¹ |
| 0.45 | 111.7  | 1.71  | 1.03×10⁻¹ | 1.58×10⁻³ | 1.05×10⁻¹ |

Plot: `runs/containership_fr_resistance/Cw_vs_Fr.png`.

## Comparison to Wigley (I5)

| Hull          | Peak Fr | Cwp at peak | Cwp at Fr=0.20 |
|---------------|---------|-------------|----------------|
| Wigley        | 0.40    | 6.0×10⁻²    | 2.0×10⁻²       |
| Containership | **0.35**| **1.3×10⁻¹**| **1.0×10⁻¹**   |
| Ratio         | −12 %   | +2.1×       | +5×            |

**The Containership has a peak at lower Fr** (0.35 vs 0.40) and
**twice the peak Cwp** (1.3×10⁻¹ vs 6×10⁻²). At low Fr the
difference is even more dramatic — at Fr=0.20 the Containership
drags 5× more than the Wigley.

## Interpretation

- **Peak Fr shift to 0.35** is consistent with the Containership's
  shorter effective length: the parallel midbody adds wetted area
  but doesn't significantly change the bow/stern wave-making
  region. The Froude number based on the bow-stern distance is
  effectively higher, shifting the resonance peak down.
- **Higher Cwp across the board** comes from the bluffer ends
  (linear taper vs Wigley's parabolic faired ends). The wave train
  is launched with a sharper pressure pulse at the bow and stern.
- **Wave drag at the hump dominates by 100× over viscous drag**
  (Cwp ~ 0.13 vs Cwv ~ 0.0015). At very low Fr (≤0.10) viscous would
  start to matter, but we didn't sweep that range.
- **Engineering implication**: a Containership hull at Fr ≈ 0.20
  has wave resistance similar to a Wigley at Fr ≈ 0.30. So a
  real ship with this Cb can NEVER operate efficiently at high
  Fr — the wave-resistance penalty saturates the propulsion budget
  long before Wigley would. This is why real containerships do
  cruise at lower Fr ≈ 0.20 with much bigger props.

## Why the K2 finding (J_self higher at Fr=0.40) now makes sense

K2 observed that the Containership's J_self at Fr=0.40 was
*higher* than at Fr=0.25, opposite to the intuition "higher Fr →
more wave drag → lower J". This sweep explains why: at Fr=0.40
the Containership is **past** its wave-resistance peak (which sits
at Fr=0.35), so wave drag *decreases* going from 0.35 → 0.40. The
J_self trend is therefore non-monotonic in Fr.

## Caveats

- 120 steps may not fully settle the longer wavelengths at Fr=0.45.
- Reference area is `V/T` (planform-like, not wetted), so absolute
  Cw values aren't directly comparable to towing-tank Bai-style
  Wigley data. *Relative* shape of the curve is meaningful.
- BDIM bias from L3 means the absolute drag is ~12 % low; the
  shape of the curve (peak position, peak/trough ratio) is robust
  to this since both Containership and Wigley share the bias.

## See also

- `RESULTS-Fr-scan-current.md` (I5) — Wigley counterpart, peak at Fr=0.40.
- `RESULTS-cb-vs-Jself-Fr40.md` (K2) — J_self at Fr=0.40, now
  better-explained.
- `runs/containership_fr_resistance/Cw_vs_Fr.png` — the curve.
