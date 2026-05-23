# RESULTS — Propellers.jl actuator-disk validation

Fifth cross-validation: Propellers.jl uniform-thrust actuator disk vs
1D Froude–Rankine momentum theory in a uniform stream. This is the
Layer-1 gate of the [Propellers.jl PLAN](../Propellers.jl/PLAN.md).
OpenFOAM cross-validation (Layer 2) requires a working OF Docker
image which is currently blocked on arm64 (no qemu-user-static).

## Headline

C_T target = 0.5, 128 × 48 × 48 cells, periodic y/z, exit BC in x.
Run to steady state (300 mom_step!).

| Quantity      | Theory  | Measured | Error    |
|--------------:|--------:|---------:|---------:|
| U at disk     | 1.1124  | 1.0801   | −2.9 %   |
| U at 5R aft   | 1.2247  | 1.1712   | −4.4 %   |
| induction a   | 0.1124  | 0.0801   | −28.7 %  |

**PASS** — within the Layer-1 gate of U_disk ±5 %, U_wake ±10 %.

The induction-factor error (28 %) is large because a = U_disk/U∞ − 1 ≈ 0.1
amplifies small absolute errors in U_disk. The headline U values are
the more honest metric.

## Why theory and convention

Two competing conventions for actuator-disk theory:

| Convention | Body force direction | Wake | Induction | C_T (definition)         |
|-----------:|---------------------:|-----:|----------:|:--------------------------|
| Wind turbine | opposes flow         | slower | a = 1 − U_disk/U∞ | C_T = 4a(1−a) |
| Propeller   | same as flow         | faster | a = U_disk/U∞ − 1 | C_T = 4a(1+a) |

Propellers.jl ActuatorDisk applies a positive `thrust` value as a body
force along the propeller `axis`. Mathematically this matches the
**propeller convention** (Carrica et al. 2010 *J. Ship Res.*). My first
attempt used wind-turbine signs and gave a 27 % wrong-direction error
before fixing the theory line.

Result: U downstream of disk is FASTER than U∞ — the disk is doing
work on the fluid. Wake velocity ratio ≈ 1.17.

## What was actually run

- 3D box 128 × 48 × 48 cells
- Disk: R=6, R_hub=0, w=2, centred at i=43 (one-third along x)
- Uniform inflow U∞ = 1.0 (cell-units), Re_D = U·D/ν = 5000 (laminar)
- Periodic in y,z; convective exit BC in x
- Body force = ActuatorDisk(...; thrust=½·ρ·A·U²·C_T = 28.27 cell-units)
- 300 mom_step! suffice for steady state on this grid

Driver: [`scripts/actuator_disk_uniform.jl`](./scripts/actuator_disk_uniform.jl).

## Outstanding follow-ups

| # | Item                                                                  |
|---|------------------------------------------------------------------------|
| 1 | OF Docker cross-validation (Layer 2). Needs qemu-user-static on host. |
| 2 | Goldstein (radially-graded) actuator disk                              |
| 3 | Actuator line (Sørensen & Shen 2002)                                  |
| 4 | DTC hull + actuator disk = Layer 3 (Phase 3 milestone)                 |

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
JULIA_NUM_THREADS=auto WL_NX=128 WL_NY=48 WL_NZ=48 WL_R=6 WL_NSTEPS=300 \
    julia +1.12 --project=. scripts/actuator_disk_uniform.jl
```
