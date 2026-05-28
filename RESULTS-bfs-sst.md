# RESULTS — backward-facing step (SST), Layer-3 — WIP

Status: **harness scaffolded, not yet stable.** The k–ω SST backward-
facing-step (Driver–Seegmiller) Layer-3 validation is set up in
[`scripts/bfs_sst.jl`](scripts/bfs_sst.jl) but the run diverges early
(|u|→NaN before ~1 flow-through) at both Re_H = 37500 and 5000, so no
reattachment length is reported yet.

## What's there

- 2D expansion-ratio-2 geometry: inlet channel of height H above a step
  of height H; floor drops from y=H to y=0 at `x_step`. BDIM walls.
- SST stepped with its native ω-wall treatment (the channel showed the
  Spalding override double-counts for SST, so it is *not* used here).
- Reattachment detection: sign change of the near-wall streamwise
  velocity on the bottom wall downstream of the step.
- Target: x_r/H ≈ 6.0–6.3 (Driver–Seegmiller experiment; kΩSST
  literature 6–7).

## Why it diverges (diagnosis, not yet fixed)

The instability appears within the first ~2000 steps, independent of Re,
so it is numerical/BC-driven, not a physical high-Re effect:

1. **Inlet discontinuity.** The uniform `uBC=(U,0)` injects `U` across
   the *whole* inlet plane, including into the solid step block, where
   BDIM forces it back to zero — a strong shear singularity at the step
   lip on the inflow boundary. Fix: pass `uBC` as a function that
   injects `U` only above the step (`y>H`), zero below.
2. **Sharp re-entrant corner.** The `min`-of-two-planes floor SDF has a
   kink at the step corner; BDIM + the SST gradient terms (vorticity,
   ω cross-diffusion) spike there. Fix: round the corner SDF slightly,
   and/or cap ν_t.
3. **CFL / explicit coupling.** The segregated explicit SST substep at
   the corner may need a tighter Δt than the momentum CFL allows.

## OF cross-check — unblocked (use the arm64-native image)

Earlier this looked blocked: the locally-cached
`openfoam/openfoam11-paraview510` image is linux/amd64 and this host is
arm64 (`exec format error`). That was the *wrong image* — it's the
Foundation `.org` build, amd64-only.

The harness-prescribed image (`HARNESS.md`) is
**`opencfd/openfoam-default:2406`**, which is **multi-arch with a native
arm64 build** (`linuxARM64GccDPInt32Opt`). Verified here: `simpleFoam` /
`pimpleFoam` run natively, and the `pitzDaily` backstep tutorial
(`$FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily`) runs end-to-end
(blockMesh + simpleFoam) with no emulation.

```
docker run --rm --platform linux/arm64 opencfd/openfoam-default:2406 ...
```

Image digest at time of writing:
`sha256:dd5aa20630a55722663bf83ba0cb74870cba130081303e32e3865007fa2aa35a`.

So the OF `pitzDaily` kΩSST cross-check **can run on this host** after
the WaterLily BFS harness is stabilised; no amd64 machine needed.

## Next steps (tracked)

1. Function inflow BC (inject only above the step).
2. Corner SDF smoothing + ν_t cap.
3. Re-run; measure x_r/H vs the experimental 6.0–6.3.
4. (If an amd64 host is available) OF pitzDaily kΩSST cross-check.

This is the one remaining open item from the RANS validation plan; the
channel-flow law-of-the-wall gates for both SA (`RESULTS-channel-sa.md`)
and SST (`RESULTS-channel-sst.md`) are passed.
