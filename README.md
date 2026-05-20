# ShipFlow.jl

Integration repo for the ship-flow CFD stack on top of
[WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl). Private,
WIP.

This repo holds:

- The **application** — end-to-end ship simulation scripts (DTC,
  KCS resistance, DTC self-propulsion).
- The **validation harness** — OpenFOAM-in-Docker orchestration,
  reference-data plumbing, comparison metrics, performance tracking.
- The **MASTER_PLAN.md** that sequences the five sibling repos
  (Turbulence.jl, VoF.jl, Propellers.jl, ShipShapes.jl, and the
  upstream-hooks PR in `pankgeorg/WaterLily.jl`).

See [MASTER_PLAN.md](./MASTER_PLAN.md).
