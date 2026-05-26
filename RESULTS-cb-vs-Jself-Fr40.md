# Cb sweep at Fr = 0.40 (K2)

Repeats I1 (Cb sweep) at Fr = 0.40 where I5 placed the Wigley
wave-resistance peak. Tests whether the J_self(Cb) relationship is
sensitive to the operating Fr.

Driver: [`scripts/containership_cb_sweep_Fr40.jl`](scripts/containership_cb_sweep_Fr40.jl).

## Setup

Identical to I1 except `Fr = 0.40` (was 0.25). Grid 96 × 48 × 48,
Containership hull, 100 steps per (par_frac, J), drag averaged over
last 25 %.

## Result

| par_frac | Cb    | J_self @ Fr=0.40 | J_self @ Fr=0.25 (I1) | Δ        |
|----------|-------|------------------|----------------------|----------|
| 0.30     | 0.65  | **0.2592**       | 0.2106               | **+23 %** |
| 0.50     | 0.75  | **0.2264**       | 0.1762               | **+28 %** |
| 0.70     | 0.85  | **0.1808**       | 0.1448               | **+25 %** |
| 0.85     | 0.93  | **0.1375**       | 0.1106               | **+24 %** |

**Empirical fit at Fr = 0.40:** `J_self ≈ 0.534 − 0.430 · Cb`
(slightly steeper than the Fr=0.25 slope of −0.367).

## Interpretation

- **J_self is higher at Fr = 0.40** by ~25 % across the entire Cb
  family. Lower drag → less rotor thrust required → higher J.
- **This is counter-intuitive relative to Wigley I5**, which found
  the wave-resistance peak at Fr = 0.40. Two possible explanations:
  1. **Containership ≠ Wigley.** I5 was run on a Wigley hull
     (L = 36, B = 8, T = 5). K2 uses the Containership family. The
     Containership's wave-resistance peak likely sits at a different
     Fr — parallel-midbody hulls have broader and lower wave-resistance
     curves than slender Wigleys.
  2. **Including the rotor changes the picture.** The rotor jet
     induces a low-pressure region behind the hull that *reduces*
     the hull's pressure drag. At higher Fr the jet's contribution
     scales with `ρU²` faster than the wave drag grows, so the net
     drag drops.
- **The Cb scaling is preserved**: relative ranking of par_frac
  is the same as Fr=0.25. The slope `dJ/dCb` is slightly steeper
  at Fr=0.40 (−0.43 vs −0.37).
- **Practical implication**: at design Fr = 0.40, the rotor needs
  to be ~30 % smaller (or running at higher RPS) than at Fr = 0.25
  for the same Cb. This kind of cross-Fr scaling rule is exactly
  what a preliminary-sizing tool would want.

## Caveats

- 100 steps may be marginal at Fr = 0.40 — the wavelength scales
  with Fr², so the wave train past the stern is longer than at
  Fr = 0.25, and 100 steps × 0.25 dt = 25 cell-time-units may not
  fully resolve the steady wake.
- VLM CT values for par_frac=0.93 reach 134; these are well in the
  unphysical regime (see K3 / `RESULTS-rotor-cavitation.md`). The
  *trend* across Cb is still meaningful but the absolute thrust
  numbers should not be quoted in isolation.
- 96 grid only; per J4 the J_self values shift downward at finer
  grids. Absolute J_self ~0.13–0.26 should be treated as biased
  high by ~10 %.

## Cross-Fr comparison at par_frac=0.5 (Cb=0.75)

Adding the M3 data point at Fr=0.35 (the actual Containership
wave-resistance peak per L5) completes the J_self vs Fr picture
for this hull:

| Fr       | J_self  | D at J_self | Notes |
|----------|---------|-------------|-------|
| 0.25 (I1) | 0.176   | ~400        | "off-peak" Fr, fully on rising flank |
| **0.35** (M3) | **0.210** | ~310    | wave-res peak per L5 |
| 0.40 (K2) | 0.226   | ~280        | past peak |

Non-monotonic: J_self rises monotonically with Fr in this hull
range, but D_hull at J_self is lowest at Fr=0.40 and highest at
Fr=0.25. The rotor + free-surface coupling produces a different
trend than the bare-hull L5 resistance curve, because the rotor's
race modifies the hull's pressure field.

## See also

- `RESULTS-cb-vs-Jself.md` — I1, the Fr=0.25 baseline.
- `RESULTS-Fr-scan-current.md` — I5, the Wigley resistance curve.
- `RESULTS-Fr-scan-containership.md` — L5, the Containership
  wave-resistance curve (peak at Fr=0.35).
- `runs/cb_sweep_Fr40/scans.csv` — detail data.
- `runs/containership_selfprop_Fr35/scan.csv` — M3 single-point.
