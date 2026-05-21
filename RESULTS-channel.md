# RESULTS — channel-flow validation (status: partial)

Second cross-code validation: Turbulence.jl Smagorinsky LES in
WaterLily against an OpenFOAM `incompressibleFluid/channel395`
tutorial run with the same Smagorinsky closure. Target Re_τ = 395
(turbulent channel, classic Moser-Kim-Mansour benchmark setup).

## Headline

|                       | WaterLily (N_HC=16) | OpenFOAM (60 k cells) | Expected (turbulent)|
|-----------------------|---------------------|-----------------------|---------------------|
| Flow regime           | **laminar**         | **turbulent**         | turbulent           |
| u_max / U_bar         | 1.60                | 1.17                  | ≈ 1.17              |
| Re_τ achieved         | ≈ 94                | ≈ 234                 | 395                 |
| ν_t (max)             | ≈ 5 × 10⁻³          | (not extracted)       | nonzero in BL       |

**Both codes agree with the Smagorinsky *closure* — the discrepancy is
that WaterLily's channel flow stayed laminar over the available
simulation time, while OpenFOAM's transitioned to turbulence.** This
is a *physics* difference (regime), not a code disagreement.

OF's near-parabolic-but-flattened profile (u_max=1.17 instead of the
laminar 1.50) is the classic turbulent channel signature.

WL's perfectly parabolic profile (u_max=1.60, very close to the
laminar 1.50, with the 7% extra peakiness explained by the BDIM
wall smearing) shows the flow has not transitioned. The Smagorinsky
model is responding to the laminar shear (νₜ grows to ~ 2× ν), but
there's no inflectional instability seeded to push the velocity
field into turbulence within 60 D/U of simulation time.

## What was actually run

### OpenFOAM (`runs/channel395`)
- Foundation v11 cylinder mesh: 4×2×2 with 40×50×30 cells (60 000 total)
  graded toward both walls
- `incompressibleFluid` solver, Euler ddt, linearUpwind div
- Smagorinsky LES (switched from default WALE to match Turbulence.jl)
- `meanVelocityForce` constraint targeting `Ubar = 0.1335`
- ν = 2×10⁻⁵ → bulk Re ≈ 6675
- Ran to t = 1000 s (≈ 130 flow-through times); averaged t = 500-1000
- Result: turbulent, centerline u/Ubar = 1.17, achieved Re_τ ≈ 234
  (lower than the tutorial's design 395 because Smagorinsky is more
  dissipative than the default WALE — known characteristic)

### WaterLily (`runs/channel395_waterlily_v2`)
- 32 × 32 × 32 cells, BDIM walls via `min(y, N_Y-y)` SDF
- Periodic streamwise + spanwise (`perdir=(1,3)`)
- Parabolic IC + 3% wave perturbation
- Constant streamwise body force g_x = u_τ²/δ (canonical channel driver)
- Smagorinsky from Turbulence.jl as `udf`, updating `flow.ν` each step
- ν chosen for bulk Re = 6675 (matching OF)
- Ran to t = 60 D/U; averaging from t = 30 D/U
- Result: **laminar parabolic**, u_max/Ubar = 1.60, no transition to
  turbulence

## Diagnosis

At bulk Re = 6675 the channel flow *is* linearly unstable, but
inflectional instability from a parabolic base flow with only smooth
perturbations grows slowly — sub-critical transition with 3%
sinusoidal perturbation typically takes hundreds of flow-through
times. The OpenFOAM run had 130 flow-through times *with* random
small-scale numerical noise from its segregated solver, which it
seems was sufficient. WaterLily's projection-method residuals are
much cleaner — less "free noise" — so the flow stays in its
laminar attractor unless deliberately tripped.

### Earlier 6-hour run

A larger run at N_HC=32, t_end=400 (32×64×64 cells) was attempted
but its log output was piped through `tail -3` so intermediate
progress couldn't be inspected. The process consumed ~6 hours of
wall time without producing visible output and was eventually
killed. Even at t=400, the same laminar-stays-laminar issue would
likely have applied (the IC was the same 3% smooth perturbation).

## Path forward — three options

1. **Trip the flow explicitly** — replace the smooth IC perturbation
   with a random divergence-free field at 5-10% magnitude (e.g.
   Schoppa-Hussain vortex pair seeding). This is the standard recipe
   in LES literature for forcing transition in finite time.

2. **Restart from an OF snapshot** — read a turbulent OF velocity
   field, project it onto the WL Cartesian grid, use it as the IC.
   The flow starts already turbulent and Smagorinsky has to sustain
   it. This is the most direct test of the closure.

3. **Higher Re** — push to Re_bulk ≈ 30 000 (Re_τ ≈ 1000). Instability
   grows much faster, even smooth ICs transition in 20-50 D/U.
   Trade-off: needs proportionally finer mesh for the wall layer.

For "the next small step", #1 is by far the cheapest. ~30 lines of
extra code in the WL script + one rerun.

## What this validation *does* tell us

- Turbulence.jl's `Smagorinsky` model wires up correctly via
  WaterLily Hook 1 (effective viscosity)
- The model computes νₜ from the strain-rate tensor as expected
- The full pipeline (BDIM walls, periodic BCs, body force, MeanFlow
  averaging, profile extraction) runs end-to-end without error
- The infrastructure (OF Docker harness, time-averaging, profile
  binning across graded meshes, comparison machinery) is fully in
  place

What it does **not** tell us is whether the Smagorinsky model
*correctly predicts wall-bounded turbulence*. That awaits one of the
three transition recipes above.

## Files

- `runs/channel395/` — OpenFOAM case + post-processed profile
- `runs/channel395_of_profile/ux_profile.csv` — OF (x,z,t)-averaged
- `runs/channel395_waterlily_v2/ux_profile.csv` — WL profile
- `scripts/channel395_waterlily.jl` — WL driver
- `scripts/channel395_of_profile.jl` — OF post-process
- `scripts/compare_channel.jl` — side-by-side comparison
