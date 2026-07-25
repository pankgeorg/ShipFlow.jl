# RESULTS — blade-resolved rotating DTMB 4381, open water (Phase-3 Layer-3)

**First blade-resolved rotating propeller in the stack**, run in BOTH
WaterLily (moving-body BDIM, analytic SDF) and OpenFOAM v2512 (MRF,
snappyHexMesh) off **one geometry definition** —
`NavalArchitectToolbox.blade_sdf(dtmb4381)` + a capped-cylinder hub; the
OpenFOAM STL is marching cubes over that same SDF
(`NavalArchitectToolbox.jl/examples/propeller_stl.jl`, 389 k triangles,
h = 0.95 mm sampling).

OpenFOAM is the **native arm64 v2512** install (`/usr/lib/openfoam/openfoam2512`)
— the old "amd64-images-only" cross-validation blocker is gone.

## Verdict (design point J = 0.889)

| source | Re_D | KT | 10KQ | η |
|---|---|---|---|---|
| **Experiment** (Boswell 1971, via SMP'24) | 6.2×10⁵ | **0.208** | **0.445** | **0.661** |
| NAT VLM (`openwater_vlm(dtmb4381,…)`)¹ | ~inviscid | 0.174 | 0.388 | 0.634 |
| OpenFOAM k-ω SST, 1.6 M cells, no layers | 6.2×10⁵ | 0.154 (−26 %) | 0.363 (−18 %) | 0.601 (−9 %) |
| OpenFOAM laminar (matched-Re twin) | 5×10³ | 0.078 | 0.489 | 0.225 |
| WaterLily BDIM, D = 64 cells (matched twin) | 5×10³ | 0.112 ± 0.030 | 0.102 ± 0.084 | (1.56 — see below) |

¹ VLM leading-edge-suction constant `C_LE = 0.80` was calibrated on the
**skewed 4382** (where it hits the experiment to ~1 %); on the unskewed
4381 it under-predicts KT by ~16 %. A skew-dependence of the calibration,
worth its own note.

### The accidental first code-to-code agreement

The rotation handedness was initially wrong **in both codes in the same
way** (+Ω about the downstream shaft axis; the correct sense is −Ω:
the section apparent wind `U ẑ − Ωr θ̂` only aligns LE→TE for Ω < 0).
In that reversed, bluff, pressure-dominated braking state the two codes
agreed to ≈5 %:

|  | KT | 10KQ |
|---|---|---|
| WaterLily +Ω | −0.75 | 1.48 |
| OpenFOAM +Ω | −0.79 | 1.47 |

That regime is insensitive to thin-blade resolution — worth recording as
evidence the two force integrations, geometries, and kinematics are
mutually consistent.

## Honest read, working state (−Ω)

- **OpenFOAM SST at experimental Re: KT −26 %, KQ −18 %, η −9 %.**
  Uniform under-prediction with the classic under-resolved-LE signature:
  the STL was sampled at h = 0.95 mm while the root LE radius is
  ~0.35 mm, the snappy mesh has **no prism layers** (y⁺ ≈ 100, wall
  functions) and 1.6 M cells vs the ~5.4 M body-fitted reference RANS of
  Wu & Kinnas (SMP'24) which lands within a few %. Improvement path is
  mechanical: finer STL sampling, `addLayers`, more cells.
- **Matched-Re (5×10³) laminar pair — the real code-to-code test:**
  KT 0.112 (WaterLily) vs 0.078 (OpenFOAM): same ballpark, both far below
  the high-Re values because Re_chord ≈ 3×10³ decambers the sections
  (OpenFOAM's full laminar J-sweep shifts the zero-thrust J from ~1.25
  down to ~1.05 — see table below).
- **WaterLily torque is NOT credible at any affordable D** — the ladder
  settles it:

  | | KT | 10KQ |
  |---|---|---|
  | WL D = 64 | 0.112 ± 0.030 | +0.102 ± 0.084 |
  | WL D = 96 | 0.132 ± 0.019 | **−0.044 ± 0.039** (sign flips into noise) |
  | OF laminar (same Re) | 0.078 | 0.489 |
  | VLM (inviscid) | 0.174 | 0.388 |

  The blade is sub-cell-thick outboard of ~0.4R (t(0.7R) ≈ 0.6–0.9
  cells; immersed-volume ratio 0.14 at D = 64), so BDIM captures no
  blade-surface shear and the pressure-torque contribution cancels into
  noise. **KT meanwhile converges AWAY from the viscous low-Re answer
  and TOWARD the inviscid VLM** — coherent: BDIM with an unresolved
  boundary layer behaves quasi-inviscid. A credible KQ needs several
  cells across the section (D ≳ 256 — not a CPU option here, but a
  **single-A100 GPU job**: shrunk 2.5D×1.5D×1.5D domain at D = 256 is
  94 M cells ≈ 10 GB of fields; WaterLily's KA backend makes it
  `mem=CuArray` away. Uniform-grid resolution is a *hardware* limit, not
  a solver one — what WaterLily can't do is *concentrate* resolution).

  **D = 128 hits a hard stability wall — the CPU ladder ends here.**
  Two independent failure modes, same onset (rev ≈ 0.5, healthy
  KT ≈ 0.19 before it):
  - defaults (ϵ = 1): Δt collapses, KT → −30 (divergence);
  - stabilized (ϵ = 2, Poisson `tol=1e-5` via the solver-control kwargs
    chain): no divergence, but a **saturated spurious state** — KT climbs
    monotonically to ≈ 2.0 (10× physical) with η staying plausible
    (KT and KQ grow together). Signature of a pumping feedback: the
    ≥kernel-porous blade acts as a fan, through-flow accelerates,
    effective J at the blade drops, loading rises, repeat.
  - D = 64 is flat-to-decreasing across revs and D = 96 only mildly
    drifting — the runaway is specific to the ~1-cell-thickness regime.
  - **λ = vanLeer probe: fails identically.** Onset at the same
    rev ≈ 0.5, later full divergence (KT → O(10³), Δt → 1e-4). The
    limiter changes the failure trajectory, not the cause.
  - **The instability is robust to every affordable knob** — (ϵ 1→2,
    Poisson tol 1e-4→1e-5, λ quick→vanLeer) all fail with the same
    onset. Conversely, the *pre-onset* value is strikingly repeatable:
    KT ≈ 0.19 / 0.19 / 0.20 across the three D = 128 attempts (rev
    0.3–0.45 window) — a consistent early-window estimate right on the
    →inviscid trend (VLM 0.174).
  - **Final call:** blade-resolved WaterLily is closed pending a
    kernel-consistent thin-shell BDIM treatment (upstream-grade work —
    the ~1-cell-thickness pumping feedback is a publishable finding in
    itself). GPU D = 256 is on hold: added resolution sharpens, not
    cures, the feedback. The stack's propeller lane remains the
    validated actuator/VLM disk; blade-resolved verification remains
    OpenFOAM's.

  **Division-of-labor conclusion:** blade-resolved propellers belong to
  OpenFOAM (or an AMR/multi-resolution solver); WaterLily's propeller
  lane in this stack is the validated actuator/VLM-coupled disk
  (±0.3 % round-trip, RESULTS-propeller-layer2.md) behind the hull —
  with OpenFOAM supplying/verifying the radial loading.

## OpenFOAM laminar J-sweep (Re_D = 5×10³, 1.6 M cells)

| J | KT | 10KQ | η | VLM KT | VLM 10KQ |
|---|---|---|---|---|---|
| 0.600 | 0.220 | 0.638 | 0.33 | 0.287 | 0.546 |
| 0.889 | 0.078 | 0.489 | 0.23 | 0.174 | 0.388 |
| 1.100 | −0.044 | 0.337 | (windmilling) | 0.082 | 0.258 |

## Numerics that bit (worth remembering)

1. **Handedness**: fixed analytically (velocity triangle) + empirically
   (both codes reverse). A 20-step frozen-rotor force read is a
   *transient*, not a thrust sign.
2. **WaterLily time bases**: `sim_time()` is convective (t·U/L); Ω, the
   ramp and rev counting are flow-time. Mixing them made the loop 64×
   too long (and mimicked a "CFL collapse" — Δt was healthy all along).
3. **Spin-up ramp** (quadratic, 0.25 rev): removes the impulsive-start
   spike of full-Ω on sub-cell blades.
4. **OpenFOAM parallel snappy workflow**: `decomposePar` before snappy
   prunes not-yet-existing patches from the fields — keep fields in
   `0.orig`, `restore0Dir -processor` after meshing, and add
   `#includeEtc "caseDicts/setConstraintTypes"` for procBoundary*.
5. **Marching-cubes STL hygiene**: weld vertices (adjacent-voxel float
   scatter), keep the dominant connected component (~400 few-triangle
   specks from the razor-thin tip), leaving 1 part / 0 illegal /
   4 open edges (one ~1 mm pinhole snappy doesn't care about).

## Reproduce

```bash
# geometry (STL for OpenFOAM; WaterLily uses blade_sdf directly)
julia --threads=32 NavalArchitectToolbox.jl/examples/propeller_stl.jl

# WaterLily, matched-Re: D cells, J, revs via env
WL_D=64 WL_REVS=2 julia --project=ShipFlow.jl --threads=32 \
    ShipFlow.jl/scripts/prop_bladeresolved_openwater.jl

# OpenFOAM: laminar matched-Re twin, or SST at experimental Re
ShipFlow.jl/openfoam/dtmb4381_openwater/Allrun 0.889
TURB=sst ShipFlow.jl/openfoam/dtmb4381_openwater/Allrun 0.889 <target-dir>
```

Data: `runs/prop_bladeresolved/kt_kq_*.csv`, `runs/of_dtmb4381_J*/postProcessing/forces/`.
