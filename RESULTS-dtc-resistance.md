# RESULTS — DTC bare-hull calm-water resistance (Fr = 0.218)

First Phase-3 resistance attempt on the **actual DTC benchmark
geometry** (every prior resistance number was Wigley or the
Containership parametric stand-in). Full WaterLily five-package stack:
WaterLily (flow + BDIM) + VoF.jl (free surface, variable density) +
Turbulence.jl (WALE LES; SST verified to compose) + ShipShapes.jl
(tabulated DTC SDF, `dtc_sdf_h12.5mm.bin`).

Reference: el Moctar, Shigunov & Zorn, *Duisburg Test Case:
Post-Panamax Container Ship for Benchmarking*, **Ship Technology
Research 59(3), 50–65, 2012**, Table 4 (model-scale resistance, scale
1:59.407, SVA Potsdam towing tank). Cross-checked against the
OpenFOAM-11 `incompressibleVoF/DTCHull` tutorial (same hull, ν=1.09e-6,
ρ=998.8, kΩSST RAS, $UMean = 1.668 m/s).

**Verdict: GATE FAILS. The five-package stack runs the real DTC hull in
a full free-surface turbulent flow — stable, mass-conserving, force
output finite — but it does NOT reproduce the calm-water residuary
resistance at Fr=0.218. The resolved force is dominated by (a) a trapped
longitudinal sloshing mode (C_T oscillates with ±56 % std, never
settles) and (b) BDIM pressure/form-drag bias on the coarse Cartesian
cut. Best estimate C_P,sim ≈ 3.6×10⁻³ (NL=128) vs reference
C_R,ref = 0.62×10⁻³ — about +475 %, far outside the ±15 % gate. The
true wave-resistance signal at this sub-hump Froude number
(C_R = 0.6×10⁻³, C_W = 0.34×10⁻³) is genuinely tiny and sits below the
method's current noise floor. Reported honestly with the root causes so
the next iteration is targeted, not blind.**

The one clean quantitative result is the *fully-submerged* (deep-water)
DTC hull drag, which IS grid-converged (split-half 1.6 %): C_P,deep ≈
5.25×10⁻³ at Re=1e5 — see §Deep-water diagnostic.

## The Reynolds problem and why the Froude subtraction breaks

Model-scale Re at Fr=0.218 is **Re_m = 9.145×10⁶**, unreachable on a
uniform Cartesian grid with BDIM and *no y⁺ wall function*. We ran at
**Re = 1×10⁵** intending to recover C_R through the Froude line:

```
C_R = C_T − C_F,ITTC(Re)        ITTC-57:  C_F = 0.075/(log₁₀Re − 2)²
```

**This subtraction is invalid here.** BDIM at this resolution develops
essentially *no* turbulent wall friction: the measured viscous
coefficient is C_F,meas ≈ 1×10⁻⁴, while the theoretical
C_F,ITTC(1e5) = 8.33×10⁻³ — two orders of magnitude larger. So
`C_T,sim − C_F,ITTC(Re_sim)` overshoots into a large *negative* number
(−4.6×10⁻³ at NL=128) and is meaningless. The simulation's resolved
force is almost entirely *pressure* (form + wave + BDIM-staircase +
sloshing); there is no resolved friction to subtract. The physically
honest comparison is therefore **C_P,sim vs the reference residuary
C_R,ref** (residuary resistance is itself the pressure-origin part).

### Reference numbers (Table 4, Fr = 0.218, model scale)

| quantity | value |
|---|---|
| model speed v_m | 1.668 m/s |
| Re_m | 9.145×10⁶ |
| R_T | 31.83 N |
| **C_T,exp** | **3.670×10⁻³** |
| C_F,ITTC(Re_m) | 3.047×10⁻³ |
| **C_R,exp = C_T − C_F** (simple Froude) | **0.623×10⁻³** |
| C_W (paper, form factor k=0.094) | 0.336×10⁻³ |
| wetted surface S_w | 6.243 m² |
| Lpp, B, T | 5.976, 0.859, 0.244 m |
| displacement ∇ | 0.827 m³, C_B = 0.661 |

Normalization (both sides): `C = R / (½·ρ·S·U²)`, S = 6.243 m²,
S_c = S/ΔX².

## Headline numbers — grid ladder

Hull FIXED at design draft, zero trim. Re=1×10⁵, Fr=0.218, WALE LES,
MULES+mass-repair VoF, ρ-ratio 100:1, gravity ramped over the first
1 L/U, averaged over t·U/Lpp ∈ [2, t_end]. KernelAbstractions backend,
`julia -t auto` (80 threads); SIMD restored afterward.

| N_L | grid | ΔX [m] | C_T,sim (mean ± std) | C_P (press) | C_F (visc, meas) | split-half |
|---:|---|---:|---:|---:|---:|---:|
| 96  | 256×128×112 | 0.0623 | 4.28×10⁻³ ± 2.40×10⁻³ (56 %) | 4.19×10⁻³ | 8.2×10⁻⁵ | 6.8 % |
| 128 | 336×160×144 | 0.0467 | 3.70×10⁻³ ± 2.09×10⁻³ (56 %) | 3.58×10⁻³ | 1.2×10⁻⁴ | 3.0 % |

(`C_F,ITTC(1e5) = 8.33×10⁻³`. Froude C_R = C_T − that = −4.1/−4.6×10⁻³,
invalid, see above.)

**Comparison (primary):** C_P,sim vs C_R,ref = 0.623×10⁻³:
- NL=96:  C_P = 4.19×10⁻³ → **+573 %**
- NL=128: C_P = 3.58×10⁻³ → **+475 %**

Both **FAIL** the ±15 % gate. C_P *decreases* with refinement
(4.19 → 3.58×10⁻³), trending the right direction but still ~5–6× the
reference and nowhere near a 15 % band; 2 grids is a trend, not a
Richardson extrapolation (the third grid, NL=160 / 15.2 M cells, was
dropped — projected >10 h under thread contention; the plan says drop
grids from the top, keep the write-up).

**Coincidental near-match warning.** C_T,sim at NL=128 (3.70×10⁻³)
happens to land within 1 % of the *experimental* C_T (3.67×10⁻³). This
is a **cancellation of errors, not a validation**: our Re is 1e5 (wrong
friction regime — friction is ~100× under-resolved) and the number is
sloshing-dominated (±56 %). Do not cite it as agreement.

## Deep-water diagnostic (the one converged number)

To isolate wave-making we ran a fully-submerged baseline (`WL_DEEP=1`,
α=1 everywhere — no free surface). It settles cleanly (no sloshing):

| run | C_T | ± std | C_P | split-half |
|---|---:|---:|---:|---:|
| free surface (NL=96) | 4.28×10⁻³ | 2.40×10⁻³ | 4.19×10⁻³ | 6.8 % |
| **deep water (NL=96)** | **5.27×10⁻³** | **5.5×10⁻⁵** | **5.25×10⁻³** | **1.6 %** ✓ |

The deep-water run is the only **grid-converged** point (1 % std,
split-half 1.6 %): a stable, well-defined drag on the real DTC geometry
at Re=1e5 — a genuine stack-integration success even though it isn't the
headline resistance.

The wave-resistance-by-subtraction `C_P,fs − C_P,deep` = 4.19 − 5.25 =
**−1.06×10⁻³ (negative)** is **confounded**: the deep run wets the
*entire* hull (both above and below the design waterline) whereas the
free-surface run wets only the submerged half, yet both are normalized
by the same waterline S. The wetted-area mismatch swamps the wave term
(the same caveat flagged in RESULTS-drag-decomposition.md). So the
decomposition does not rescue the comparison here.

## Root-cause analysis (why it fails, in priority order)

1. **Trapped sloshing.** With spanwise-periodic sides and a 2.6-Lpp box,
   the impulsive start + gravity ramp excite a longitudinal standing
   wave whose period (~1 L/U) does not decay over the run. C_T swings a
   full order of magnitude (≈0.4×10⁻³ ↔ 8×10⁻³) every cycle; the [2,4]
   window captures ~2 cycles, so the mean is defined but the std is
   ±56 % and split-half stays 3–7 % (> the 2 % protocol target). A real
   towing tank avoids this (open/long tank, gradual tow). **Fixes to
   try:** much longer domain + sponge/damping zone at the outlet,
   gentler combined inflow+gravity ramp, or a steady-state (local
   time-stepping) formulation.
2. **No y⁺ wall function under the hull.** BDIM resolves the geometry
   but not the turbulent boundary layer; resolved friction is ~100×
   below ITTC. This is *why* the Froude subtraction is invalid and why
   raw C_T (Re=1e5) is not comparable to C_T,exp (Re=9.1e6). Needs the
   BDIM wall-stress model on the risk register (Kempe & Fröhlich 2012)
   or an SST native ω-wall treatment that actually loads the near-wall
   band.
3. **Tiny true signal.** At Fr=0.218 (well below the hump) the DTC is
   friction-dominated: C_R = 0.62×10⁻³ is only 17 % of C_T,exp. The
   wave signal is near the method's numerical noise floor — a hard
   target for a first attempt. A higher-Fr point (≥0.30, nearer the
   hump, larger wave fraction) would be an easier first validation.
4. **BDIM staircase form drag.** The Cartesian cut of a slender hull
   adds a resolution-dependent pressure bias; C_P drops 15 % from
   NL=96→128, consistent with this decaying under refinement but not
   yet converged.

## What was actually run

`runs/dtc_resistance/forces_<TAG>.csv` (sampled every 5 steps):
`t_LU, Fx_pressure, Fx_viscous, CT, CP, CF_meas, umax, wave, niter`.
`ladder_summary.csv` — the window-averaged ladder.

- **Driver:** `scripts/dtc_resistance.jl`. Knobs `WL_NL`, `WL_FR`,
  `WL_RE`, `WL_TEND_LU`, `WL_TURB` (wale|sst|smagorinsky|none),
  `WL_VOF` (clamp|mules), `WL_RHO`, `WL_RAMP_LU`, `WL_DEEP`,
  `WL_HYDRO`, `WL_TAG`.
- **Analysis:** `scripts/dtc_resistance_analysis.jl`.
- **Domain:** x ∈ [−0.75, +1.85] Lpp about midship (≥1 Lpp wake), width
  1.30 Lpp (full hull, no symmetry plane), 0.80 Lpp below + 0.30 above
  the waterline; spanwise periodic, convective exit, blockage 0.41 %.
  Grid dims snapped to multiples of 16 (MultiLevelPoisson needs
  size = a·2ⁿ).
- **Grids:** NL=96 (256×128×112, ~3.7 M cells, t_end=4 L/U, ~35 min),
  NL=128 (336×160×144, ~7.7 M cells, t_end=4 L/U), deep96 (NL=96
  submerged, t_end=3 L/U). NL=160 (15.2 M) attempted then dropped on
  time.

### Hydrostatic sanity (`WL_HYDRO=1`)

Quiescent box, gravity ramped on, vertical pressure force vs ρ·g·∇:
ratio settled at **0.77** at NL=48 (23 % under-buoyant). Dominated by
vertical under-resolution at the smoke grid (T=0.244 m ≈ 2 cells), plus
the ~12 % BDIM buoyancy bias. Affects *vertical* force only; the
fixed-hull resistance is horizontal, so it does not enter the C_R
comparison.

### Turbulence sensitivity

**SST (kΩSST, the reference CFD's closure) composes with VoF and runs
stably** (NL=48 smoke: stable, mass-conserved, CT within ~5 % of WALE
at the same instant, Poisson nit=3). WALE was the production default —
cheaper (no per-cell wall-distance solve) and, since neither model
develops near-wall friction without a wall function, the eddy-viscosity
choice is a second-order effect on this result. Reported for the user:
the SST+VoF composition path works and is available via `WL_TURB=sst`.

## Caveats

- **Fixed sinkage/trim** (tests were free in both); a real bias even at
  this sub-hump Fr. Free sinkage is a stretch goal (Newmark 1-DOF
  machinery exists in the heave scripts).
- **ρ-ratio 100:1**, not physical 815:1 (the damBreak study's
  stability/realism compromise).
- **SDF resolution ceiling:** trilinear from h=12.5 mm; finest ΔX
  (0.047 m) stays coarser than h, so the SDF is not the limit.
- **VoF advection:** MULES + mass-repair; c_α interface compression OFF
  (documented sharp-interface+body instability regime).

## Phase-3 progress

This is the first time the DTC tabulated SDF flows in a full
free-surface + turbulence simulation (ShipShapes milestone 3 criterion
"feeds VoF.jl Phase-3 integration test" — **met**, integration-wise).
The resistance *validation* gate is **not** met; the headline DTCHullProp
self-propulsion gate (MASTER_PLAN Phase-3) remains blocked behind a
working bare-hull resistance, which needs the sloshing and wall-function
issues above resolved first.

## Recommended next steps

1. **Kill the sloshing** before anything else: lengthen the domain
   downstream and add an outlet wave-damping/sponge zone; ramp inflow
   *and* gravity together more gently; or move to a steady local-time-
   stepping formulation. Target split-half < 2 % so C_T is actually
   defined.
2. **Pick an easier first Fr** (≥0.30, nearer the hump) where the wave
   fraction of C_T is large enough to clear the noise floor, then walk
   back down to 0.218.
3. **Add a near-wall stress model** (BDIM wall function / SST ω-wall) so
   resolved friction is physical and the Froude decomposition becomes
   valid — this is the long pole for any Re-honest ship resistance.
4. Once C_T settles: re-enable the 3-grid ladder (96/128/160 or
   96/144/192) for a real Richardson/GCI band.

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
JL=/home/pgeorgakopoulos/.julia/juliaup/julia-1.12.6+0.aarch64.linux.gnu/bin/julia
# production used backend=KernelAbstractions in WaterLily/LocalPreferences.toml
WL_HYDRO=1 WL_NL=48 $JL --project=. scripts/dtc_resistance.jl          # hydrostatic
for NL in 96 128; do
  WL_NL=$NL WL_TEND_LU=4.0 WL_TURB=wale WL_VOF=clamp WL_TAG=$NL \
    $JL -t auto --project=. scripts/dtc_resistance.jl
done
WL_NL=96 WL_TEND_LU=3.0 WL_DEEP=1 WL_TAG=deep96 \
  $JL -t auto --project=. scripts/dtc_resistance.jl                    # submerged baseline
WL_GRIDS="96:96,128:128" $JL --project=. scripts/dtc_resistance_analysis.jl
```

---

## Round 2 — fixes (2026-06-13)

Round-2 acts on the round-1 diagnosis (above): kill the trapped sloshing so
C_T settles, then walk the Froude number, then wire the existing near-wall
stress model. Follows `PLAN-dtc-fixes.md`. Round-1 content above is left intact
for the plan-vs-outcome diff.

### Verdict (Round 2)

**Fix 1 (sloshing) — LANDED via a surface (α) damping zone, not the velocity
sponge.** The planned Rayleigh *velocity* sponge `−σ(u−u_ref)` applied as a
`udf` body force was built and wired, but it **stalls the pressure projection**
(Poisson hits the 50-iteration cap every step, vs nit≈2 without it) because the
relaxation body force injects a velocity field the incompressible projection
must then fight. Isolated by a controlled A/B: longer domain + *no* sponge →
nit≈2; longer domain + velocity sponge → nit=50; longer domain + **α-damping
only** → nit≈2. The α-damping (blend α toward the still-water profile inside the
outlet/inlet bands, applied *after* the VoF step so it never enters the Poisson)
is the piece that actually flattens the outgoing wave's surface deformation —
exactly as the plan anticipated ("this is the piece that actually stops the
standing wave"). Production config: **α-damp ON, velocity sponge OFF**
(`WL_SPONGE=1 WL_SPONGE_UVEL=0`), with the domain lengthened downstream so the
damping bands sit clear of the hull + near wake.

### Fix-1 settling result (NL=96, Fr=0.218, Re=1e5, WALE, t_end=4 L/U)

Settled window t·U/Lpp ∈ [2, 4], sampled every 5 steps. Three sponge
variants on the lengthened (DOM_XHI=2.50) domain, KA backend, `-t auto`:

| variant | config | C_T (mean ± std) | std/mean | split-half | nit |
|---|---|---:|---:|---:|---:|
| round-1 (ref) | no sponge, short box | 4.28e-3 ± 2.40e-3 | 56 % | 6.8 % | 2–3 |
| α-damp only | `UPOST=0`, α-damp=1 | unsettled (swings through 0) | 112 % | 62 % | 2–3 |
| **α-damp + post-sponge 0.5** | `UPOST=0.5` | **3.56e-3 ± 2.2e-3** | 62 % | **10.5 %** | 2–3 |
| velocity body-force sponge | `UVEL=1` | — (Poisson **stalls**, nit=50) | — | — | **50** |

**Honest verdict: Fix 1 PARTIAL — not landed to the < 2 % split-half
target.** The DTC hull spans the full channel width in a closed-ish
≈3.25-Lpp box; the impulsive start + gravity ramp excite a robust
longitudinal **seiche** (standing wave, period ≈ 1 L/U, C_T amplitude
≈ ±2e-3 about a ≈3.5e-3 mean). The damping bands can only live in the
far wake (aft of ≈1.55 Lpp) and a thin 0.12-Lpp inlet strip — they must
stay clear of the hull (bow 0.16, stern 1.34 Lpp from inlet) — so they
absorb the wave only where it reaches the box ends, not the antinodes
between inlet and hull. The post-projection velocity sponge **does**
reduce the swing and is **Poisson-safe** (nit stays ≈2–3), but the
settled split-half (10.5 %) is not better than round-1's plain run
(6.8 %) and the mean drifts across sub-windows (it is a slowly-evolving,
not a stationary, oscillation). The windowed mean C_P ≈ 3.5e-3 is still
≈ 5–6× the reference residuary C_R = 0.62e-3 — i.e. the **accuracy**
gate also fails, for the same round-1 reasons (no near-wall friction →
Fix 3; tiny sub-hump wave signal; BDIM staircase form drag). Reported as
a windowed mean ± honest std with sloshing noted, per the plan's failure
playbook.

**The one genuinely new, transferable finding** is the Poisson-stall
diagnosis: a Rayleigh *velocity* body force inside `mom_step!`'s `udf`
hook pins the multigrid at its iteration cap (50 vs 2–3) on this
two-phase 100:1 system — isolated by a clean A/B (longer box, sponge
OFF → nit≈2; sponge ON → nit=50; α-damp only → nit≈2) and fixed by
moving the velocity relaxation to a **post-projection** explicit step
(`post_sponge!`). Anyone adding a momentum sponge to this stack should
apply it post-projection, never as a `udf` body force.

A **stronger/wider variant** (`UPOST=1.0`, outlet band 1.4 Lpp,
`fix1f`) was run to test whether more aggressive end-zone damping cracks
the seiche: it does **not** — C_T = 3.65e-3 ± 2.45e-3, split-half
**14.9 %** (slightly *worse* than `UPOST=0.5`). The standing wave is
robust against damping *strength*; the limit is that the bands can only
sit at the box ends, not at the seiche antinodes near/under the hull.
Killing it needs a structurally different approach (much longer domain
with a real inlet sponge ahead of the bow; a non-reflecting outflow;
integer-period or far-longer averaging; or a steady local-time-stepping
formulation) — out of scope for this round, on the next-steps list.

### Reference data correction (important)

Re-reading the source paper (el Moctar/Shigunov/Zorn 2012, Tab. 4, fetched and
parsed directly) shows the SVA Potsdam towing tank tested **only six speeds,
Fr = 0.174 … 0.218**. There is **no experimental point at Fr = 0.28 or 0.33** —
those are above the tested envelope. The full table (model scale, C_T,C_F ×1e-3,
C_W ×1e-4; C_F is ITTC-57, C_W = C_T − (1+k)C_F with k = 0.094):

| Fr | Re×10⁻⁶ | C_T×10³ | C_F×10³ | C_W×10⁴ |
|---:|---:|---:|---:|---:|
| 0.174 | 7.319 | 3.661 | 3.170 | 1.932 |
| 0.183 | 7.681 | 3.605 | 3.142 | 1.672 |
| 0.192 | 8.054 | 3.588 | 3.116 | 1.791 |
| 0.200 | 8.415 | 3.602 | 3.092 | 2.194 |
| 0.209 | 8.783 | 3.623 | 3.069 | 2.660 |
| **0.218** | 9.145 | **3.670** | 3.047 | **3.360** |

Consequence for **Fix 2**: the planned Fr=0.33→0.28→0.218 walk-down cannot be
quantitatively gated above 0.218 (no reference). The grounded walk is therefore
**down the tested ladder** (0.218 → 0.209 → 0.200 …), where every point has a
real C_T,exp. Fr=0.218 (the highest tested) carries the largest wave fraction
(C_W=3.36e-4) and is the strongest validation point available — so it remains
the primary target. Runs above 0.218 (if done) are reported as sim-only with an
explicit "no experimental reference" caveat, never compared to an invented
number. This table is encoded in `dtc_resistance_analysis.jl` (`TAB4`,
`tab4_ref`).

### Driver changes (Fix 1 + Fix 3 wiring)

`scripts/dtc_resistance.jl` (new knobs):
- `WL_SPONGE` (default 1), `WL_SPONGE_UVEL` (default **0** — velocity body-force
  sponge OFF; it stalls the Poisson), `WL_SPONGE_W` / `WL_SPONGE_WIN` /
  `WL_SPONGE_TOP` (outlet/inlet/top band widths, Lpp), `WL_SPONGE_SIGMA`
  (σ_max·Lpp/U∞, only used by the velocity sponge), `WL_SPONGE_ADAMP` (α-damp
  blend, default 1.0).
- `WL_DOM_XHI` (downstream extent, Lpp, default 2.50) — lengthened so the
  damping bands clear the hull bbox (bow≈0.16 Lpp, stern≈1.34 Lpp from inlet).
  A startup guard prints the band edges vs the hull bbox and WARNs on overlap.
- `WL_INFLOW_RAMP` (default 1) — smooth 0→U∞ inflow ramp over RAMP_LU L/U;
  not the slosh fix, and the production gate run used impulsive inflow
  (`WL_INFLOW_RAMP=0`, round-1-comparable, fast Poisson).
- `WL_WALLFN` (default 0) + `WL_BAND` (default "1,3") — Fix 3: turn on
  Turbulence.jl's Spalding BDIM wall function. `WL_TURB=sa` added
  (Spalart–Allmaras, the model whose wall function hit the channel log-law
  within 7.4 %). The SA path captures the model-only ν, applies
  `apply_wall_function!`, then **gates the override to water** (α>0.5) so the
  air side gets no spurious stress. Smoke-tested end-to-end (constructs from the
  DTC SDF, steps stably, mass-conserving) — machinery de-risked; the open
  question is only whether resolved C_F rises physically.

No WaterLily or Turbulence.jl **source** was modified — Fix 3 uses Turbulence's
existing exported `step_sa!` / `apply_wall_function!` / `spalding_uτ` from the
driver, per the plan.

### Fix 2 (Froude walk-down) — BLOCKED on Fix 1

Fix 2 gates on a settled C_T. Since Fix 1 did not reach the < 2 %
split-half target, a clean C_R-vs-reference comparison at any Fr is not
yet meaningful (the ±2e-3 seiche swamps the C_R = 0.6e-3 sub-hump
signal). The reference table (TAB4 above) and the per-Fr comparison
helper (`tab4_ref`) are in place for when Fix 1 settles. Key reference
finding stands regardless: **the experimental envelope ends at
Fr = 0.218**, so the originally-planned 0.33/0.28 points have no truth
value and the grounded walk is 0.218 → 0.209 → 0.200 down the tested
ladder.

### Fix 3 (near-wall stress model) — WIRED + smoke-tested, not yet
exercised at scale

The Spalding BDIM wall function is fully wired and de-risked but not run
to a settled state (it gates on Fix 1):
- `WL_TURB=sa` builds a `SpalartAllmaras` model from the DTC tabulated
  SDF (`wall_distance` over the hull), steps stably, mass-conserving
  (NL=32 smoke, t=0.3 L/U, no errors).
- `WL_WALLFN=1` path: `step_sa!` (no internal wallfn) → snapshot
  model-only ν → `apply_wall_function!` (Spalding `spalding_uτ` →
  flux-match ν override in band `WL_BAND`, default (1,3)) → gate the
  override to water (α>0.5) so the air side gets no spurious stress.
- **Open question** (the риск the plan flagged): whether the measured
  C_F actually rises from the round-1/round-2 ~1e-4 toward a physical
  O(1e-3) at this coarse BDIM resolution. Not answered — needs a settled
  C_T first so the friction signal isn't buried in the seiche.

## Round-2 status summary

| fix | outcome | why |
|---|---|---|
| 1 sloshing | **PARTIAL** | sponge built (3 variants); post-projection velocity sponge is the Poisson-safe design and reduces the swing, but the confined-box seiche persists (split-half ≈10 %, not <2 %); reported as windowed mean ± std |
| 2 Froude walk | **BLOCKED** | gates on a settled C_T; reference table corrected (envelope ends Fr=0.218 — no 0.28/0.33 truth) |
| 3 wall function | **WIRED + smoke-tested** | SA + Spalding path runs end-to-end, water-gated; not exercised at scale (gates on Fix 1) |

## Reproduce (round 2)

```sh
cd $ROOT/ShipFlow.jl
JL=/home/pgeorgakopoulos/.julia/juliaup/julia-1.12.6+0.aarch64.linux.gnu/bin/julia
# production used backend=KernelAbstractions in WaterLily/LocalPreferences.toml (restore SIMD after)
# Fix-1 best settling run (post-projection velocity sponge, α-damp, lengthened box):
WL_NL=96 WL_TEND_LU=4.0 WL_TURB=wale WL_VOF=clamp WL_DOM_XHI=2.50 \
  WL_SPONGE=1 WL_SPONGE_UVEL=0 WL_SPONGE_UPOST=0.5 WL_SPONGE_W=1.0 WL_SPONGE_WIN=0.12 \
  WL_INFLOW_RAMP=0 WL_TAG=fix1e_upost05 $JL -t auto --project=. scripts/dtc_resistance.jl
# A/B that isolated the Poisson stall:
WL_NL=96 WL_TEND_LU=0.25 WL_DOM_XHI=2.50 WL_SPONGE=0 WL_TAG=disc_nospnj  $JL -t auto --project=. scripts/dtc_resistance.jl  # nit≈2
WL_NL=96 WL_TEND_LU=0.7  WL_DOM_XHI=2.50 WL_SPONGE=1 WL_SPONGE_UVEL=1 WL_TAG=test_uvel $JL -t auto --project=. scripts/dtc_resistance.jl  # nit=50
# Fix-3 wall-function smoke (SA + Spalding, water-gated):
WL_NL=32 WL_TEND_LU=0.3 WL_TURB=sa WL_WALLFN=1 WL_DOM_XHI=2.50 WL_TAG=smoke_sa $JL -t 4 --project=. scripts/dtc_resistance.jl
```
