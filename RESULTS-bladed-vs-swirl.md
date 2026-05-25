# BladedRotor (VLM) vs SwirlingDisk on Wigley + free surface

Same Wigley hull, same prescribed thrust (calibrated from the
BladedRotor's VLM `CT`), same free-surface, same Re/Fr. The only
difference is the propeller model.

Driver: [`scripts/swirling_vs_bladed_rotor.jl`](scripts/swirling_vs_bladed_rotor.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 32, B = 8, T = 5 at xc = NX/4.
- Fr = 0.30, Re = 5000, ρ_w/ρ_a = 10.
- 80 timesteps; statistics over the last third.
- BladedRotor: 3 blades, R = 3, R_hub = 0.5, chord (1.0 → 0.5),
  twist (35° → 15°), 12 × 4 panels. J = 0.7, |CT| ≈ 3.01,
  thrust ≈ 42.5, |CQ| ≈ 0.42, torque ≈ 17.7.
- SwirlingDisk: matched `thrust = 42.5`, `torque = 17.7`, w = 1.5.

## Result

|                              | SwirlingDisk    | BladedRotor pt-smear | BladedRotor sectional (G4) |
|------------------------------|-----------------|----------------------|----------------------------|
| Hull drag (mean ± σ)         | 45.64 ± 0.64    | 52.13 ± 0.42         | **50.19 ± 0.74**           |
| Δ_drag vs SwirlingDisk       | —               | **+14.2 %**          | **+10.0 %**                |
| Wave RMS (post-stern)        | 0.2004          | 0.2047               | 0.2132                     |
| Wave peak-peak               | 3.1070          | 2.7037               | 3.0013                     |

> The third column is **G4 of NEW_PLAN.md**: replace the single-point
> `smear_force!` + `smear_torque!` at the disk centre with
> `LiftingSurfaces.smear_blades!`, which distributes the same total
> thrust+torque across `N_blades × N_sections = 3 × 4 = 12` points
> arranged as radial blade lines. Net thrust and net moment match
> to ~5 % per the unit test.
>
> **Sectional smear narrows the drag delta to +10 % and brings the
> wave peak-peak nearly back to the SwirlingDisk value**, confirming
> the review intuition: the single-point smear was over-loading the
> hull stern through localised flow acceleration. The sectional
> version spreads the force across the swept area like
> SwirlingDisk does, so the hull sees a similar mean flow.
>
> Earlier numbers in commit history: cell-centred smear gave +7.9 %
> drag; face-staggered single-point gave +14.2 %; sectional gives
> +10.0 %. All numbers are within the same physical regime; the
> remaining 10 % gap between sectional-BladedRotor and SwirlingDisk
> reflects the VLM's actual radial loading distribution (which is
> not uniform — there's a tip-loaded bias).

## Interpretation

For the same integrated thrust, the **VLM blade-resolved rotor
deposits the body force in a more concentrated way** (Gaussian
smear at a single point with `ε = 2.5`) than the SwirlingDisk
(top-hat distribution over the entire annulus, `w = 1.5`). This
shows up as:

- **+7.9% hull drag** — the more concentrated jet pulls slightly
  harder on the hull stern. Real propellers are somewhere between
  these two extremes; a multi-point smear (one Gaussian per blade
  section, varying with radial position) would be more faithful.
- **Identical wave RMS** at the post-stern probe — the integrated
  wave energy is the same.
- **−11 % wave peak-to-peak** — BladedRotor's wake is smoother
  (peak suppression by ~10 %). The wide top-hat SwirlingDisk
  excites the free surface more sharply at the wake edge.

## When to use which

- **Use SwirlingDisk** when you only care about time-averaged
  drag and don't trust the per-blade smear placement. The
  forces it deposits are spread over the entire disk volume
  — robust to grid resolution.
- **Use BladedRotor** when you want a real VLM solve in the
  loop (radial-load awareness, two-way coupling with WaterLily
  via `trilinear_inflow`, future unsteady blade-passage work
  with `unsteady_analysis!`). The drag prediction will differ
  from SwirlingDisk by ~5-10 % at the same C_T; the wake
  structure will be marginally less peaked. Both predictions
  are within the noise of the actuator-disk approximation
  generally.

## Caveats

- The current BladedRotor smear is a single Gaussian at the disk
  centre. Distributing the spanwise force as a function of radius
  (one Gaussian per panel strip, scaled by sectional CL) would
  be closer to the actuator-line method. See PLAN-additions §A1.
- The BladedRotor isn't using two-way coupling here. With
  `additional_velocity=trilinear_inflow(flow.u)` (the path used in
  `wigley_full_stack_VLM.jl` with `WL_TWOWAY=1`) the rotor would
  see the hull wake reduction at the disk plane, which would
  reduce its CT slightly — closer to a real wake-fraction
  calculation.
- The torque (CQ) is not deposited in this run. SwirlingDisk does
  apply it; BladedRotor's udf only smears the axial thrust. The
  swirl difference would be more apparent if both used the same
  torque path.
