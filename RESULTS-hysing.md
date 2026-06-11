# RESULTS — Hysing rising bubble (CSF surface-tension validation)

VoF.jl PLAN milestone 3 Layer-1 benchmark: Hysing et al. 2009 test
case 1 (ρ-ratio 10, μ-ratio 10, Re=35, Eo=10, D=0.5 bubble in
[0,1]×[0,2], t→3 s). Driver: `scripts/hysing_bubble.jl` —
MULES + c_α=1 advection, Brackbill CSF via the `udf` hook
(`σκ∇α·(1/ρ_face)` with the face density shared with the projection),
4 smoothing passes for κ.

**Verdict: rise dynamics VALIDATED, shape fidelity NOT in envelope.**

| metric | reference (Hysing groups) | N=64 | N=128 |
|---|---|---|---|
| max rise velocity | 0.2417–0.2421 @ t≈0.92 | **0.2417** @ 0.85 | **0.2405** @ 0.85 |
| centroid y_c(t=3) | 1.0813–1.0817 | 1.0965 (+1.4 %) | 1.0946 (+1.3 %) |
| min circularity | 0.9011–0.9013 | 0.636 | 0.605 |
| mass drift | — | 1e-6 | 1e-6 |

- **Rise velocity magnitude is spot-on** (0.0–0.6 % at both grids);
  the peak arrives ~8 % early. Centroid lands within 1.4 %. The
  buoyancy–drag balance through the variable-density projection and
  the capillary force scale are right.
- **The bubble deforms too much** and refinement makes it slightly
  worse (0.636 → 0.605) — systematic, not resolution. The circularity
  estimator itself carries only ≈5 % staircase bias (a perfect initial
  circle reads 0.947–0.955), so most of the deficit is real
  deformation.

## Why the shape is off (analysis — updated after the HF experiment)

Two hypotheses were tested:

1. ~~**Smoothed-CSF curvature under-restores the rim.**~~ **Ruled out
   2026-06-11.** VoF.jl gained Popinet-style height-function curvature
   (`curvature!(vof, Val(:height))` — on a fractional disc its
   face-sampled κ is within 1.9 % with ~30 % lower scatter than the
   smoothed estimate). Re-running N=128 with `WL_KAPPA=height` gives
   **c_min = 0.6052 vs 0.6047 smoothed** — unchanged. The κ estimator
   is not the binding error.
2. **Tangential viscous stress jump is unmodeled — now the isolated
   cause.** Our momentum equation is kinematic with ν = μ(α)/ρ(α); for
   this case ν is then *uniform* (10/1000 = 1/100), so the 10×
   dynamic-viscosity contrast never enters — and the missing
   `∇·(μ∇uᵀ)` transpose term is exactly the piece that carries the
   tangential stress balance at an interface with a μ jump. The
   bubble's internal circulation is over-driven, and it over-flattens
   regardless of how accurately the capillary force is evaluated.

Fixing (2) means a conservative variable-μ stress formulation — the
mass-momentum-consistent territory of InterfaceAdvection.jl, out of
VoF.jl's algebraic scope by design. The benchmark quantifies the cost
of that design choice: ~0.3 of circularity at Eo=10, with rise
dynamics unaffected.

## Implication for ship scales

Hull wave-making operates at Weber numbers where σ is irrelevant
(VoF PLAN: surface tension matters for propeller cavitation, not
resistance). The validated parts — buoyancy/drag through the
projection, capillary force scale, exact mass — cover the ship use
case. The shape envelope matters if we ever do bubble/cavity dynamics
seriously; then the path is height-function κ (+ transpose viscous
term), filed as follow-ups.

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
WL_N=64  WL_TEND=3.0 WL_TAG=N64  julia +1.12 --project=. scripts/hysing_bubble.jl
WL_N=128 WL_TEND=3.0 WL_TAG=N128 julia +1.12 --project=. scripts/hysing_bubble.jl
```

CSVs: `runs/hysing_bubble/bubble_N{64,128}.csv`
(t, y_c, v_c, circularity, area_rel, mass_rel).
