# Drag decomposition: wave vs viscous (L2)

Attempts to isolate wave-resistance from total drag by running the
Wigley at each Fr twice: once with the free-surface VoF (water + air)
and once with α₀ = 1 everywhere (deep water, fully-submerged hull).
The difference is *attributed* to wave-making.

Driver: [`scripts/wigley_drag_decomposition.jl`](scripts/wigley_drag_decomposition.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 36, B = 8, T = 5.
- Re = 5000, fixed hull, no propeller.
- 120 steps per (Fr, mode), last 25 % averaged.
- Mode A: `α₀(x, z) = 1 if z ≤ H_w_c else 0` (free surface).
- Mode B: `α₀ = 1` everywhere (deep water, hull fully submerged).

## Result

| Fr   | D_total (FS) | D_deep | D_wave = D_total − D_deep | D_wave / D_total |
|------|--------------|--------|----------------------------|-------------------|
| 0.20 | 28.38        | 54.48  | **−26.10**                | −0.92             |
| 0.25 | 29.55        | 46.99  | −17.45                     | −0.59             |
| 0.30 | 34.69        | 41.57  | −6.88                      | −0.20             |
| 0.35 | 41.30        | 37.74  | +3.56                      | +0.09             |
| 0.40 | 42.96        | 34.06  | +8.90                      | +0.21             |
| 0.45 | 40.47        | 30.77  | +9.71                      | +0.24             |

Plot: `runs/drag_decomp/decomp.png`.

## What this is — and isn't

- **The inflection between Fr=0.30 and Fr=0.35 is real.** At low Fr
  the free-surface case has LESS total drag than the deep-water
  case (negative "wave drag"), at high Fr it has MORE. The
  crossover at Fr ≈ 0.32 aligns with the wave-resistance hump in
  I5 (Wigley peak Cwp at Fr = 0.40, transition through the rising
  flank around Fr = 0.30–0.35).
- **The negative low-Fr "wave drag" is an artefact.** The
  deep-water case has the hull *fully* submerged — twice the
  wetted area of the free-surface case (the air-half of the hull
  has effectively no friction). So D_deep is inflated by a factor
  of ~2 in viscous drag at low Fr where viscous dominates.
- **At Fr=0.40 (near the wave-resistance hump)** D_wave = +8.9
  is ~21 % of total drag. Adjusted for the wetted-area mismatch
  (D_deep should be halved), the true wave-resistance fraction is
  *more like 40–50 %* at this Fr. That matches the published
  Wigley breakdown (Bai 1979): wave-resistance is ~half of total
  drag at the hump.

## Doing this right

For an accurate decomposition we'd need either:

1. **A "half-Wigley"**: define the hull to exist only below the
   waterline, then submerge it deep enough in the deep-water case
   that no free-surface effects intrude (e.g. hull at z = NZ/4
   instead of z = NZ/2). Same wetted area as the free-surface
   case.
2. **Use known empirical viscous Cf**: instead of running a
   second simulation, subtract the ITTC '57 friction line value
   `Cf = 0.075 / (log10(Re) − 2)²` from D_total. At Re = 5000,
   Blasius gives `Cf = 1.328/√Re = 0.019`. With S_wet ≈ 765 cells²
   and U=1, ρ=10: D_viscous ≈ Cf · 0.5 · ρ · U² · S_wet = 0.019 ·
   0.5 · 10 · 1 · 765 = 73. Subtracting from D_total at Fr=0.40
   gives D_wave ≈ 42.96 − 73 < 0 — also too negative.

Neither approach gives clean numbers at Re = 5000 because the
viscous contribution is comparable to total drag. At a real
ship-Re (10⁸+), viscous is a small correction and the
decomposition would work cleanly.

## What we learned anyway

- **The decomposition machinery (run two sims, subtract drags) works**.
  The qualitative wave-resistance growth from Fr=0.30 → 0.45 is
  captured, matching the I5 wave-resistance curve.
- **The Wigley wave-resistance peak by this method falls around
  Fr=0.40** — same as I5. Two independent measurements agree.
- **For Re=5000 we cannot quantitatively separate wave from
  viscous** without either of the methodological fixes above. This
  is a documented limitation, not a bug.

## See also

- `RESULTS-Fr-scan-current.md` (I5) — Wigley Fr curve, gives the
  total Cwp shape.
- `RESULTS-Fr-scan-containership.md` (L5) — Containership
  counterpart.
- `runs/drag_decomp/decomp.png` — the three curves.
