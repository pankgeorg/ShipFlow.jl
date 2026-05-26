# Rotor cavitation regime analysis (K3)

Analytical sanity check: at every operating point we've used, is the
`BladedRotor` in the cavitation regime, the stall regime, or both?

Driver: [`scripts/rotor_cavitation_check.jl`](scripts/rotor_cavitation_check.jl).

## Setup

- Rotor geometry: same as G1 / F1 (3 blades, R = 2.4, R_hub = 0.2·R,
  tapered + twisted, 12 × 4 panels).
- Depth: 3 cells below waterline (≈ T/2 for our hull).
- Cell-units pressure: hydrostatic `p_∞ = ρ_w·g·h_depth ≈ 13.3`.
- Vapor pressure proxy: `p_v = −1.5` (a real-water analogue at
  ~25 °C scaled to our cell units).

Cavitation number:
```
σ = (p_∞ − p_v) / (0.5·ρ·V_local²)
```
Cavitation onset criterion:
```
Cp_min < −σ   ⇒   σ + Cp_min < 0   (cavitates)
```

`Cp_min` for a thin airfoil at lift CL ≈ 2π·α_eff is approximately
`−4·CL − 4·CL²` (rule of thumb for the upper-surface minimum).

## Result

| J    | CT      | CL_section_est | Cp_min_est | σ_local | σ + Cp_min | Verdict     |
|------|---------|----------------|------------|---------|------------|-------------|
| 0.10 | 134.8   | 821            | −2.7×10⁶   | 0.07    | very neg.  | CAVITATES   |
| 0.15 | 60.3    | 367            | −5.4×10⁵   | 0.15    | very neg.  | CAVITATES   |
| 0.20 | 34.1    | 208            | −1.7×10⁵   | 0.25    | very neg.  | CAVITATES   |
| 0.25 | 22.0    | 134            | −7.2×10⁴   | 0.35    | very neg.  | CAVITATES   |
| **0.32 (G1)** | 13.5 | 82       | −2.7×10⁴   | 0.51    | very neg.  | **CAVITATES** |
| 0.40 | 8.7     | 53             | −1.2×10⁴   | 0.70    | very neg.  | CAVITATES   |
| 0.60 | 4.0     | 24             | −2400      | 1.14    | very neg.  | CAVITATES   |
| 1.00 | 1.5     | 9              | −360       | 1.79    | very neg.  | CAVITATES   |

## Interpretation

- **σ values are realistic** (~0.5 at G1, ~1.8 at light load J=1.0) —
  these are *low* compared to towing-tank ship-prop conditions where
  σ is typically 0.5–2 for design-loaded screws. Our σ values
  bracket the realistic range from "heavily cavitating" to "lightly
  cavitating".
- **CL_section is absurd** at all J values we've used. Realistic
  airfoil CL_max is ~1.4–1.7 (NACA-style, with stall). Our VLM
  reports CL_section ~ 80–820. This confirms the rotor is operating
  in a fantasy regime where VortexLattice's no-stall assumption is
  badly violated.
- **Cp_min is dominated by the CL² term** at these loads, giving
  Cp_min ~ −4·CL² ~ −10000 at G1. A realistic Cp_min for a stalled
  blade is around −10 to −15.
- **All operating points cavitate** under the analytical model. The
  margin (σ + Cp_min) is so negative that no plausible refinement
  of the cavitation model would change the verdict.

## What this means for the project

- **VLM at our loading is qualitatively useful, quantitatively
  fantasy.** The right sign of CT and the right trends with J are
  correct; the absolute numbers are 10–100× too large.
- **SwirlingDisk is the right model for quantitative open-water
  studies** at our loading. We've already calibrated SwirlingDisk
  thrust/torque from VLM (RESULTS-bladed-vs-swirl.md), and the
  bladed-vs-swirl comparison was for *the same prescribed thrust*
  — that calibration step is the right place to absorb the VLM's
  quantitative bias.
- **BladedRotor's value is in coupling**: the two-way coupling path
  with `trilinear_inflow` (F2, I2, J2) responds to the local flow.
  As long as we interpret the resulting CT as a *relative shift*
  rather than an absolute prediction, the bladed model gives us
  manoeuvring physics that SwirlingDisk can't.
- **For real ship-design work**, we'd need to either:
  1. Add a stall + cavitation model to LiftingSurfaces.jl
     (e.g., post-stall airfoil polars from XFOIL data)
  2. Use VLM only for off-design polar generation at light loads
     and switch to SwirlingDisk for heavy-load self-propulsion
  3. Couple to a BEM (Blade Element Momentum) code that includes
     stall and cavitation natively (CCBlade.jl etc.)

## Caveats

- The CL_section formula assumes the loading is uniformly
  distributed across blade panels, which it isn't (VLM has a tip-
  loaded bias). A panel-by-panel cavitation check would need the
  VortexLattice circulation γ(panel), not implemented.
- The thin-airfoil Cp_min rule is conservative — actual airfoil
  pressure distributions are smoother, but the CL² term still
  dominates at high lift.
- p_v = −1.5 is a placeholder; the verdict doesn't change for any
  reasonable vapor-pressure choice.

## See also

- `RESULTS-Jself-analytical.md` (I10) — companion analytical
  caveat about the heavily-loaded regime.
- `RESULTS-bladed-vs-swirl.md` (G4) — empirical confirmation that
  SwirlingDisk reproduces BladedRotor's hull-drag effect when both
  use the same prescribed thrust.
- `project_session_findings.md` (memory) — "VLM regime caveats"
  section.
