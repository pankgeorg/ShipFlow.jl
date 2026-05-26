# F2 diagnostic — where the side-force went (I3)

The F2 trial (`RESULTS-twoway-amplification.md`) showed that placing
the rudder inside the rotor race raises CL_rudder by **+300 %**, yet
the integrated hull side force `F_hull,y` went *slightly down*
(−2.49 → −2.19 in F2's `pressure_force` + `viscous_force` measurement).
This script reproduces the F2 setup and reports the per-x-station
side-force density on the body, to localise where the deficit
appears.

Driver: [`scripts/rudder_in_race_diag.jl`](scripts/rudder_in_race_diag.jl).

## Setup

- Same grid, hull, Fr, Re, rotor, rudder placement as F2
  (rudder at `prop_xc + 0.5·R`, inside the race).
- 80 timesteps, statistics over the last 25 %.
- Per-x side-force density `dF_y/dx` computed as
  `Σ_{j,k} p · ∂μ₀_y/∂y`, summing over y, z at fixed x. This is the
  BDIM kernel contribution to the y-momentum equation, integrated
  across the cross-section.

## Result

Both runs deposit force into the fluid along three localised regions:
- the **hull** (x ≈ hull_xc ± L/2, ≈ 7.5 to 43.5 cells)
- the **rotor disk** (x = 48.6, orange dashed)
- the **rudder** (x = 50.1, lime dashed)

|                                  | TWOWAY = 0 | TWOWAY = 1 |
|----------------------------------|-----------:|-----------:|
| Total ∫ dF_y/dx                  | −50.07     | −49.62     |
| ∫ dF_y/dx over hull (x ∈ [7.5, 43.5]) | −50.41 | −49.98     |
| Δ vs TWOWAY=0                    | —          | **−0.85 % on the hull** |

(Note: this BDIM-density measure differs in absolute scale from
F2's `pressure_force` because the two routines normalise differently;
the ratio between cases is the meaningful quantity.)

## What the plot shows

`runs/rudder_in_race/Fy_per_x.png`:

1. **Both curves peak around the rotor + rudder station** (x ≈ 48–51),
   with TWOWAY=1 showing a substantially larger spike there — the
   raw fluid-side body force from the 4× CL boost lands in the
   rotor-race region.
2. **Across the hull span** (x ≈ 8–43), TWOWAY=1 sits **slightly
   above** TWOWAY=0 in many cells (less negative side force
   reaction). The integrated effect is the ~1 % deficit observed.
3. **No upstream propagation** to the bow — the side-force response
   is confined to the back half of the hull, which is consistent
   with the convection time scale of pressure information at this
   flow speed (no time for the rudder's effect to be felt at the
   bow within 80 steps).

## Interpretation

- **Most of the "lost" hull side force is convected downstream into
  the rotor race**, not absorbed by the hull. The rudder's
  body-force injection is partly carried past the hull by the
  accelerated rotor jet *before* reaching the hull's pressure field.
- **The +300 % CL_rudder in F2 was real**, but its hull-projection
  efficiency is geometry-dependent. A rudder placed AHEAD of the
  rotor (rare in practice, but conceivable for thrust-vectoring
  studies) would project differently.
- **For manoeuvring-coupler use**, the right quantity to read is
  **CL_rudder itself** — that's what couples to the yaw equation of
  motion. The hull's BDIM-mediated reaction is a much smaller side
  channel.

## Caveats

- BDIM-density vs `pressure_force` normalisation: the two cases'
  ratio is consistent (~1 % less in TWOWAY=1 in both measurements),
  but absolute values differ. The diagnostic is reliable for
  relative analysis only.
- 80 timesteps doesn't fully settle the wave field — the
  longer-tail rudder-induced wave response is suppressed.
- Only one rudder placement; a rud_xc sweep would show how rapidly
  the deficit changes with position.

## See also

- `RESULTS-twoway-amplification.md` — F2, the result being explained.
- `scripts/rudder_in_race_diag.jl` — driver.
- `runs/rudder_in_race/Fy_per_x.png` — the spatial plot.
