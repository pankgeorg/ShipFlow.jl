# RESULTS — Spalart–Allmaras RANS channel validation

Validates the Turbulence.jl **Spalart–Allmaras** one-equation RANS model
(implemented from AIAA 92-0439) on a turbulent channel at Re_τ = 395,
against the law of the wall — the canonical SA/DNS target that OpenFOAM's
SA also reproduces.

Driver: [`scripts/channel_sa.jl`](scripts/channel_sa.jl). 2D channel
(periodic streamwise, BDIM walls in y), segregated loop
`step_sa!(model, u, dt)` → `sim_step!`. RANS needs no transition trick —
the model supplies all turbulent viscosity, so a steady profile develops
directly. 25 000 steps to convergence.

## Headline

| Quantity                  | WaterLily SA | Target            |
|---------------------------|-------------:|-------------------|
| Re_τ (force balance)      | 395.0        | 395 (exact)       |
| peak ν_t/ν                | 33.5         | ~33–41 (channel)  |
| **log-layer slope 1/κ**   | **≈ 2.5 (κ ≈ 0.40)** | 2.44 (κ=0.41) |
| log-layer additive B      | ≈ 2.4        | 5.2               |
| centreline u⁺             | 16.8         | ~20.5 (DNS)       |
| mean rel. error vs log law| 17 %         | (gate ±10 %)      |

## The key result: slope right, offset shifted by the BDIM wall

The u⁺(y⁺) profile in the log layer has the **correct slope** — i.e.
the von Kármán constant the SA closure reproduces is κ ≈ 0.40 (target
0.41):

| y⁺    | WL u⁺ | log law u⁺ | Δ      |
|------:|------:|-----------:|-------:|
| 3.1   | 1.63  | 7.95       | (sublayer) |
| 52    | 12.21 | 14.86      | −2.65  |
| 102   | 13.88 | 16.48      | −2.60  |
| 151   | 14.91 | 17.44      | −2.53  |
| 201   | 15.63 | 18.13      | −2.50  |
| 250   | 16.16 | 18.67      | −2.51  |
| 299   | 16.53 | 19.11      | −2.58  |
| 349   | 16.75 | 19.48      | −2.73  |

The deviation is a **near-constant downward shift of ≈ 2.6 u⁺ units**,
not a slope error. Reading off the slope between y⁺ = 100 and 200:
`Δu⁺/Δ ln y⁺ = 1.75 / ln 2 = 2.52`, i.e. `κ ≈ 0.40` — the SA model is
producing the right log-layer mixing. What's reduced is the additive
constant `B`: `B_eff ≈ 2.4` vs the smooth-wall 5.2.

A correct slope with a reduced additive constant is the textbook
signature of **wall roughness** — and here the "roughness" is the BDIM
smeared wall. Independent confirmation: the friction velocity from the
near-wall slope (`u_τ = 0.726`) is 27 % below the exact force-balance
value (`u_τ = √(g_x·δ) = 1.000`). The BDIM wall cannot reproduce the
sharp near-wall velocity gradient, so it absorbs momentum across a
smeared band — behaving like distributed roughness that downshifts B.

## Interpretation

- **The SA closure is correct.** The log-layer slope (κ ≈ 0.40), the
  ν_t magnitude (peak ≈ 33, right for Re_τ = 395), and the force-balance
  Re_τ (exact) all confirm the transport equation, the fv1/fv2/fw
  functions, and the production/destruction balance are implemented
  faithfully from the paper. Layer-1 unit tests (fv1 limits, monotonic
  ν_t vs shear, laminar quiescence) pass independently.
- **The ±10 % log-law gate is not met (17 %)** — but the miss is
  entirely the constant B-offset from the BDIM wall, *not* a model
  error. This is the wall-treatment open question flagged in
  `Turbulence.jl/PLAN.md` ("Wall treatment with BDIM … novel, there is
  no OpenFOAM equivalent to copy").
- Same structural cause as the **LES channel** near-wall deviation
  (`RESULTS-channel.md`, `RESULTS-channel-wale.md`), but it bites harder
  in RANS because the whole profile is anchored to the wall shear, which
  BDIM smears.

## What would close the gap

A **BDIM wall function**: rather than reading `ν_t` straight from the
SA field in the smeared band, impose the log-law wall shear as a
boundary condition on the first off-wall cell (Menter-style automatic
wall treatment, or a Spalding blend). This recovers the correct `B`
without resolving the sublayer. It is genuinely new work — there is no
body-fitted-mesh recipe to port — and is the recommended next step
before SA (or any RANS model) can be quoted to ±10 % on wall-bounded
flows in this stack.

## Verdict (closure)

SA is **correctly implemented and produces the right log-layer physics
(κ ≈ 0.40, physical ν_t)**. Without a wall treatment, quantitative
law-of-the-wall agreement is limited to ~17 % by the BDIM wall (an
effective-roughness B-offset). See the wall-function results below for
the fix.

---

## BDIM wall function (Spalding-law ν_t override)

The fix for the B-offset: a Spalding-law eddy-viscosity override
(`Turbulence.apply_wall_function!`, enabled via
`step_sa!(...; wallfn=true)`). In a band of cells off the wall it solves
Spalding's universal profile for `u_τ` from the local tangential
velocity and overrides `ν_t` so the diffusive flux carries that wall
shear. The Spalding solver itself recovers `u_τ` to machine precision
across y⁺ = 1–436 (unit-tested).

### Results (Re_τ = 395, N_HC = 64, 25 000 steps)

| Configuration            | centreline u⁺ | log-layer mean err | outer log layer (y⁺ > 70) |
|--------------------------|--------------:|-------------------:|---------------------------|
| no wall function         | 16.8          | 17.0 %             | ~15 % deficit             |
| wall fn, band (1,3)       | 17.3          | 12.2 %             | improving                 |
| **wall fn, band (4,12)**  | **18.7**      | **7.4 %**          | **< 2 % (excellent)**     |
| (DNS / log law)          | ~20.5 / —     | —                  | —                         |

With the band anchored in the log layer (cells d ∈ [4,12], i.e.
y⁺ ≈ 25–74), the **outer log layer (y⁺ > 70) matches the law of the
wall to within 2 %** and the mean log-layer error drops to **7.4 %** —
under the ±10 % gate. The centreline u⁺ rises from 16.8 to 18.7 (vs the
~19.8 log-law value at y⁺ = Re_τ — within 6 %). The B-offset is
substantially recovered.

### Residual artifact

The wide band imposes a near-linear profile *within* the band
(overriding ν_t across many cells forces `du/dy ≈ u_t/d`, not the log
gradient), producing a local velocity deficit at y⁺ ≈ 25–50 — the
max-error 35 % is entirely this buffer-band feature, not the outer
flow. A **tapered/blended override** (ramp the override weight to zero
at the band edges rather than a hard switch) is the next refinement to
remove it; a single-cell anchor under-propagates (band (5,7) gave 13 %
mean), so the trade is between lift and locality.

### Verdict

The BDIM wall function **works**: it recovers the law of the wall in the
outer log layer to < 2 % and brings the mean log-layer error under the
±10 % gate (17 % → 7.4 %), confirming the B-offset was a wall-treatment
artifact, not a closure error. A tapered override to remove the
buffer-band deficit is the tracked refinement. **Milestone 3 closed on
the closure; the BDIM-wall-function milestone is substantially met
(mean gate passed) pending the taper.**

## See also

- [`scripts/channel_sa.jl`](scripts/channel_sa.jl) — driver.
- `runs/channel_sa/uplus.csv` — full u⁺(y⁺) profile.
- [`RESULTS-channel.md`](RESULTS-channel.md), [`RESULTS-channel-wale.md`](RESULTS-channel-wale.md) — the LES-channel BDIM-wall counterpart.
- `Turbulence.jl/src/rans.jl` — SA implementation.
