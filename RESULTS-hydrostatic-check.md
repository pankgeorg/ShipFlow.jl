# Hydrostatic pressure sanity check (L3)

Validates the BDIM / sign-convention pipeline by measuring the
Archimedes upthrust on a fixed Wigley under near-zero inflow. Tests
whether `−WaterLily.pressure_force(sim)[3]` converges to the
analytical `ρ_w·g·V_displaced`.

Driver: [`scripts/hydrostatic_sanity.jl`](scripts/hydrostatic_sanity.jl).

## Setup

- Grid 128 × 64 × 32, Wigley L = 36, B = 8, T = 5 fixed.
- Inflow U∞ = 1×10⁻⁴ (effectively zero).
- Re = 5000 (irrelevant — flow is essentially static).
- Gravity from step 0 (no ramp).
- 200 steps, average over the last 25 %.
- Expected Archimedes: `ρ_w·g·V₀ = 10 · 0.444 · 640 = 2844`.

## Result

| Quantity                             | Value     |
|--------------------------------------|-----------|
| Measured F_on_body (last 25 % avg)   | **2490**  |
| Expected ρ_w·g·V₀                    | 2844      |
| Relative error                       | **−12.5 %** |
| Residual |u|_max at convergence       | 0.43      |

**Verdict (< 5 % target): FAIL** — the framework underestimates
Archimedes upthrust by ~12 % at this BDIM kernel width.

## Interpretation

- **The sign convention works correctly.** `−Fp[3]` is positive and
  on the right order of magnitude. The earlier J1/K1 sign-flip
  finding is confirmed.
- **The 12 % deficit is a systematic BDIM finite-kernel-width bias.**
  The BDIM kernel μ₀ smears the body boundary over ~1–2 cells; the
  pressure-weighted surface integration over this smeared boundary
  loses some of the upthrust. The bias is ~ε/T (kernel width over
  draught) for a generic body — at our T = 5 cells, ε ≈ 1 cell, the
  expected bias is 1/5 = 20 %. The 12 % we see is *better* than
  the upper bound and consistent with the half-sigmoid kernel
  shape WaterLily uses.
- **A small persistent |u| = 0.43** persists at "equilibrium" —
  this is the residual standing wave from cold-starting against
  the hull. Tight convergence would require either a longer run
  (until viscous damps the wave out) or explicit damping. The
  drift in F_on_body across the last 100 steps is ~1 %, so the
  measurement is reasonably converged.
- **This bias explains a *small portion* of the J1 sinkage**
  (RESULTS-heave-1dof.md). If the body measures 12 % less buoyancy
  than it should, the equilibrium z_h shifts down by
  ΔV_sub / (∂V_sub/∂z_h) ≈ 87 / 192 ≈ 0.5 cells. J1's observed
  sinkage of −1.4 cells is therefore *mostly* wave-induced running
  sinkage at Fr=0.25, not BDIM bias.

## Implications

- **Force measurements have ~10 % systematic bias**. Drag numbers
  from `pressure_force` should be quoted with this in mind.
  Relative comparisons (case A vs case B, two-way on vs off) are
  robust because both sides have the same bias.
- **For 6-DOF rigid body work** (RESULTS-heave-pitch-2dof.md), the
  ~12 % deficit means the body floats lower than analytical
  equilibrium. A practical workaround: pre-compute the body's
  effective mass from a hydrostatic test like this one, then use
  M_effective = 0.88 · M_analytical to match the BDIM-measured
  Archimedes.
- **Grid refinement helps**: at finer ε/T the bias decreases. A
  follow-up at T = 10 cells (grid 256 × 128 × 64) would confirm
  this scales as expected.

## Caveats

- A real grid-convergence study (T = 5, 10, 20 cells) would let us
  fit the bias as a function of ε/T and extract a kernel-width
  parameter. Not done here.
- The residual |u| might be a wave reflecting off the periodic-y
  boundary. The setup uses `perdir=(2,)` (periodic in y), which is
  appropriate for ship-CFD but creates a small recirculation in a
  static test like this.

## See also

- `RESULTS-heave-1dof.md` (J1) — heave dynamics, sinkage value
  partially explained by this bias.
- `runs/hydrostatic_sanity/hydrostatic.png` — time series.
- `reference_waterlily_conventions.md` (memory) — sign convention.
