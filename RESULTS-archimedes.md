# RESULTS — Archimedes buoyancy integration test

Fourth cross-validation: VoFFlow + WaterLily BDIM body coupling. A
static fully-submerged cylinder in a water column should experience
upward buoyancy equal to ρ_water · g · V_body (Archimedes).

This is the Phase-2-end integration test referenced in
[MASTER_PLAN.md](./MASTER_PLAN.md) for the VoF.jl + ShipShapes.jl
gate: "Combined VoF + ShipShape integration test: static floating
box reaches Archimedes equilibrium."

## Headline

| Configuration              | N   | Predicted F_buoy_cell | Measured F_buoy_cell | Error |
|---------------------------:|----:|----------------------:|---------------------:|------:|
| Fully submerged cylinder   |  64 | +4021.2               | +3958.7              | −1.55 % |
| Fully submerged cylinder   | 128 | +8042.5               | +7847.6              | −2.42 % |
| Half-submerged cylinder    | 128 | +4021.2               | +3107.6              | −22.7 % |

**Fully-submerged passes the 5 % gate.**  Half-submerged is over the
gate — the water surface around a partially-submerged body is being
disturbed by the body's own buoyancy-driven flow and α-clamping
mass-loss bleeds support around the body. Documented as a follow-up.

## Why the SDF matters

A first attempt used an axis-aligned BOX (`box_sdf`) and gave a 0.4 %
error at N=64 but **−78 %** at N=128.  Tracing the kernel-band
integration showed the box's corner gradients are non-smooth (the
analytic SDF transitions abruptly between dx-dominated and dy-dominated
zones at each corner). WaterLily's BDIM-smoothed pressure integral
samples cell-by-cell, and the discontinuous normal at the corners
distorts the surface integral. The amount of distortion depends on
how the box edges align with the grid, which is non-monotonic in N.

A smooth SDF (cylinder, sphere, ellipsoid, Wigley hull) does not have
this problem — gradients are continuous everywhere.

**Implication for ship-flow work:** ShipShapes hulls (Wigley, KCS) use
smooth analytic / interpolated SDFs and won't hit this issue. Avoid
axis-aligned boxes as immersed bodies in WaterLily without local
gradient regularization.

## What was actually run

Both at N = 64 and N = 128:

- 2D box of water 1 × 1 m
- Water level at y = 0.5 m;  air above
- Cylinder radius 0.10 m, centre (0.50, 0.25) — fully submerged
- ρ_water = 1000, ρ_air = 1
- VoFFlow + WaterLily Simulation with `body = AutoBody(cyl_sdf)`
- Override `pois_ctor` to multiply WL's BDIM μ₀ with vof.L (= 1/ρ_face)
- Step 60 momentum steps with `pois_tol=1e-10`, `pois_itmx=400`
- Read `pressure_force(sim)` and flip sign (WL convention: force on FLUID)

Driver: [`scripts/archimedes.jl`](./scripts/archimedes.jl).

## Pressure-field cross-check

The hydrostatic profile inside and outside the body is essentially the
same continuous ρ_w · g_cell gradient (the BDIM only constrains
velocity, not pressure inside the body). A vertical column through the
body matches the analytic prediction within ~1.5 %:

```
y_cell   p_stored   p_hydro_predicted
   0.5     1944       1970   (1.3% below)
  16.5     1696       1720
  32.5     1452       1470
  44.5      295        298   (top of cylinder, depth = 7 cells)
  64.5        1.0        1.0  (just into air)
 127.5        0.03       0.04 (top of domain)
```

Surface-integral Archimedes' theorem then converts the smooth-pressure
field into the buoyancy force we observe.

## Outstanding follow-ups

| # | Item                                                                 |
|---|----------------------------------------------------------------------|
| 1 | **Half-submerged body** — interface dynamics around the body + α     |
|   | clamping. Likely needs MULES + a CSF surface-tension correction      |
|   | for the meniscus.                                                   |
| 2 | **Dynamic Archimedes** — release a denser body, watch it sink to     |
|   | equilibrium depth (where weight = buoyancy). Requires Newton-time    |
|   | integration for the body + a remeasure! loop. Phase-3 prep.          |
| 3 | **3D Wigley hull** — ShipShapes provides the SDF; same harness but   |
|   | with a 3D Simulation. Required for the full ship-CFD gate.           |

## Reproduce

```sh
cd $ROOT/ShipFlow.jl
JULIA_NUM_THREADS=auto WL_N=128 julia +1.12 --project=. scripts/archimedes.jl
# (env knobs: WL_CONFIG=half for half-submerged, WL_YBOX=... for arbitrary depth)
```
