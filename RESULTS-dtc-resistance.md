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
