# Handoff — `julia-boat-cfd` session

Goal: build a Julia CFD stack on top of WaterLily.jl that can ultimately
simulate a ship with a rotating propeller in water. Project background,
package layout, and decision gates live in [`ShipFlow.jl/MASTER_PLAN.md`](./ShipFlow.jl/MASTER_PLAN.md);
this file is the resume-from-here snapshot.

## Where we are

Six private GitHub repos under `pankgeorg`:

| Repo | Branch | State |
|---|---|---|
| [`WaterLily.jl`](https://github.com/pankgeorg/WaterLily.jl/tree/plan/upstream-hooks) | `plan/upstream-hooks` | Hook 1 (per-cell ν) + Hook 2 (`transport!`) — both implemented and used by downstream packages |
| [`Turbulence.jl`](https://github.com/pankgeorg/Turbulence.jl) | `main` | Smagorinsky LES — **validated** against OF channel395, bulk RMS 0.028, u_max within 0.5% |
| [`VoF.jl`](https://github.com/pankgeorg/VoF.jl) | `main` | α-advection + new `VoFFlow` (variable-density) — **in progress** (see below) |
| [`Propellers.jl`](https://github.com/pankgeorg/Propellers.jl) | `main` | Actuator-disk unit-tested, no OF cross-validation yet |
| [`ShipShapes.jl`](https://github.com/pankgeorg/ShipShapes.jl) | `main` | Wigley hull SDF unit-tested |
| [`ShipFlow.jl`](https://github.com/pankgeorg/ShipFlow.jl) | `main` | OpenFOAM-Docker harness + RESULTS docs |

Local working tree: `/home/pgeorgakopoulos/foam/` (= `$ROOT`)

## Closed validations

- **Cylinder Re=100** — WaterLily vs OpenFOAM vs Williamson 1996. Cd, Cl_pp, St all within ~5–10%. Forces matched; wake bubble length differs (documented BDIM kernel-vs-D effect, not a bug). See [`ShipFlow.jl/RESULTS-cylinder.md`](./ShipFlow.jl/RESULTS-cylinder.md) and [`ShipFlow.jl/RESULTS-cylinder-profiles.md`](./ShipFlow.jl/RESULTS-cylinder-profiles.md).
- **Channel Re_τ=395** — WaterLily + Turbulence.jl Smagorinsky vs OpenFOAM Smagorinsky. Bulk u_x(y) RMS 0.028, u_max/Ubar 1.180 vs 1.174 (0.5%). See [`ShipFlow.jl/RESULTS-channel.md`](./ShipFlow.jl/RESULTS-channel.md). Key insight: needed 20% broadband 3D IC perturbation to trip transition from laminar.

## Where VoF.jl crashed

The session was attempting to validate `VoF.jl` against OpenFOAM `damBreak`. Three blockers in order:

### 1. ✅ Done — OpenFOAM reference

OF damBreak runs successfully (`$ROOT/ShipFlow.jl/runs/damBreak/`). Front-position extractor at [`scripts/damBreak_of_front.jl`](./ShipFlow.jl/scripts/damBreak_of_front.jl). Required v12→v11 syntax patches: `alpha.water` IC zonal→`setFields`, fvSchemes from parent, `fvSolution` MULES key renames (see commit `9d10f8b`).

### 2. ✅ Done — Variable-density coupling in VoF.jl

`VoFFlow` struct + `step_vof!` (commit `e1dbf54` on `VoF.jl`):
- α + face-staggered `L = 1/ρ_face` (the Poisson coefficients)
- After each `WaterLily.sim_step!`, `step_vof!` advects α with the just-projected velocity, refreshes ν=μ/ρ and L from the new α, propagates L into the MultiLevelPoisson levels via `WaterLily.update!`
- 18 unit tests pass (10 new for VoFFlow)
- Critical fix: `BC!(L, ntuple(_->0, D), false, perdir)` — zero L on the ghost layer is what `WaterLily.Flow.μ₀` does and is required for no-flux walls. Without it the solver blows up.
- Equal-density case (ρ_w=ρ_a) reproduces baseline WaterLily exactly

### 3. ⛔ Where it crashed — unit-scaling bug in `damBreak_waterlily.jl`

WaterLily internally assumes ΔX=1 cell with U_ref=1 nominal. The damBreak driver I wrote was passing **SI units mixed with cell-units**: physical `g=9.81 m/s²`, physical `ν=1e-6 m²/s`, but feeding them into a WaterLily Flow whose internal time-step is in cell-units. Result: WaterLily's CFL formula `inv(max(σ)+5ν)` returns a meaningless number, and after one step `|u|` jumps to `1e14 m/s`.

#### Diagnostic showing the fix

After rescaling to **WL cell-units** (`U_ref = sqrt(g·H)`, `g_cell = g·ΔX/U_ref²`, `ν_cell = ν/(U_ref·ΔX)`, ρ as dimensionless ratios), the simulation runs without blowup. Diagnostic at end of last turn shows after step 1 with ρ_w/ρ_a=10:1:

```
g_cell = 0.062, Re_water = 30967

Pressure at column centre (x_cell=4):
  j= 2 (water bottom)  p = +4.09   u_y = 0
  j=17 (water top)     p = −0.35   u_y = −0.0043
  j=18 (air top)       p = −0.39
Hydrostatic prediction: ρ_w · g_cell · H_water = 10 · 0.062 · 15 ≈ 9.3
```

So:
- Pressure profile is **qualitatively right** (high in water, low above)
- Magnitude is **~50% of hydrostatic** — Poisson under-converges across the ρ-jump
- Velocity is bounded but flow is still falling (~½ free-fall rate instead of zero in water)
- `pois.n[end] = 1` ← **only ONE Vcycle was used by the solver** — that's the root cause

## Resume here

### Immediate next step

The MultiLevelPoisson `solver!` defaults to `tol=1e-4, itmx=32`. It exited after 1 iteration because the L₂ norm of the residual already dropped below `tol` — but the L∞ norm at the ρ-jump interface is still large. Two fixes to try in order:

1. **Tighter tolerance + more iterations**: pass `tol=1e-8, itmx=200` to the solver. Easy first test — just modify `damBreak_waterlily.jl` to call `WaterLily.solver!(sim.pois; tol=1e-8, itmx=200)` manually, but `mom_step!` calls solver! internally. Cleanest: temporarily lower `tol` inside WaterLily for this experiment, or override via a custom `mom_project!`.

2. **L∞-based convergence**: WaterLily's solver currently exits when `L₂(p) < tol`. For variable-coefficient problems the L₂ norm hides interface error. Either patch `MultiLevelPoisson.solver!` to use `L∞(p)` or accept the L₂ exit and just iterate more.

3. **Face-density interpolation**: I'm using `1/ρ_face = 0.5(1/ρ_L + 1/ρ_R)` (arithmetic mean of 1/ρ). OpenFOAM `interFoam` uses linear interpolation of ρ then inverts. Try `1/ρ_face = 1 / (0.5(ρ_L + ρ_R))` and compare.

### The damBreak driver

Located at [`$ROOT/ShipFlow.jl/scripts/damBreak_waterlily.jl`](./ShipFlow.jl/scripts/damBreak_waterlily.jl). **Already has the SI-unit bug** — it computes `dt_phys = dt * ΔX` which is wrong; needs to use `dt_phys = dt_cell * ΔX / U_ref` once the unit-scaling is sorted.

#### Rewrite outline

```julia
const ΔX    = 0.585 / NX
const G_p   = 9.81
const H_w   = 16                # water column in cells
const U_ref = sqrt(G_p * H_w * ΔX)
const G_c   = G_p * ΔX / U_ref^2     # gravity in WL cell-units
# Phase properties as dimensionless ratios:
const ρ_w = 1000.0; const ρ_a = 1.0     # same ratio as physical
# Kinematic viscosities in WL cell-units:
const ν_w = (1e-3/1000) / (U_ref * ΔX)   # μ_w/ρ_w in m²/s, divided by U·L
const ν_a = (1.8e-5/1.0) / (U_ref * ΔX)
# μ_w, μ_a passed to VoFFlow must satisfy μ/ρ matching the above:
# VoFFlow computes ν[I] = μ[I]/ρ[I] internally, so feed those.
const μ_w_cell = ρ_w * ν_w
const μ_a_cell = ρ_a * ν_a
α₀(i,x) = (x[1] < W_w && x[2] < H_w) ? 1.0 : 0.0
vof  = VoFFlow((NX,NY); α₀, ρ_w, ρ_a, μ_w=μ_w_cell, μ_a=μ_a_cell, T=Float64)
flow = WaterLily.Flow((NX,NY), (0.0, 0.0); T=Float64, ν=vof.ν,
       g=(i,x,t)-> i==2 ? -G_c : 0.0, Δt=0.1)
pois = WaterLily.MultiLevelPoisson(flow.p, copy(vof.L), flow.σ)
# Time integration: dt_phys = sim.flow.Δt[end-1] * ΔX / U_ref
```

The clean cell-unit code lived in a one-off `julia -e` diagnostic at the end of the session (see the conversation history). Migrate that into `damBreak_waterlily.jl`.

### Beyond that — toward an actual damBreak comparison

Once the Poisson converges to give hydrostatic balance in water (u_y → 0 there), the next steps:

1. **Run damBreak end-to-end** with ρ ratio 10:1 first (easier solver), then push to 1000:1 (real water/air)
2. **Extract front-position vs time** from WL (similar to the OF extractor)
3. **Side-by-side comparison** vs `runs/damBreak_front/front_vs_t.csv` (OF reference)
4. **Cross-check vs Martin-Moyce 1952** experimental data (X-vs-τ curve)

## Reproducing the validations we already have

```sh
cd $ROOT/ShipFlow.jl
julia --project=. scripts/cylinder_Re100_waterlily.jl   # WL cylinder
julia --project=. scripts/compare_cylinder_Re100.jl     # vs OF
julia --project=. scripts/channel395_waterlily.jl       # WL channel  (set WL_N_HC=16 to keep it short)
julia --project=. scripts/compare_channel.jl            # vs OF
```

OpenFOAM cases are already run; their output is committed under `runs/`. To re-run, the harness needs the `openfoam/openfoam11-paraview510` Docker image.

## Environment setup (server resume)

1. Clone `pankgeorg` repos:
   ```sh
   for r in WaterLily.jl Turbulence.jl VoF.jl Propellers.jl ShipShapes.jl ShipFlow.jl; do
       git clone git@github.com:pankgeorg/$r.git
   done
   cd WaterLily.jl && git checkout plan/upstream-hooks && cd ..
   ```
2. Julia 1.10+ (project was developed on 1.12.5). Set `JULIA_NUM_THREADS=auto`.
3. Docker for OpenFOAM: `docker pull openfoam/openfoam11-paraview510`.
4. From `ShipFlow.jl/`, run `julia --project=. -e 'using Pkg; Pkg.instantiate()'`. The `[sources]` block in `Project.toml` points at sibling repos via relative paths — make sure all six repos sit alongside each other.

## Memory files

This session has persistent auto-memory at
`/home/pgeorgakopoulos/.claude/projects/-home-pgeorgakopoulos-foam/memory/`
with entries for user role, project goal, repo layout, validation
cases, and two pieces of feedback (don't-lick-the-cookie, OpenFOAM
license). These will load automatically if you continue in this
working directory; copy them over if you move directories.

## Files modified at crash time

- `VoF.jl/src/VoF.jl` — committed (`e1dbf54`). VoFFlow + step_vof!.
- `VoF.jl/test/runtests.jl` — committed. 10 new tests.
- `ShipFlow.jl/scripts/damBreak_waterlily.jl` — committed (`9d10f8b`) **but contains the unit-scaling bug** described above. Needs the cell-unit rewrite shown in §3.
- `ShipFlow.jl/runs/damBreak/` — OpenFOAM case + reference output, committed.
- `ShipFlow.jl/runs/damBreak_front/front_vs_t.csv` — committed.

No uncommitted work in flight at crash.
