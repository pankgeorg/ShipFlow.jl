# Feature — Kato–Launder production for k–ω SST

Opt-in turbulent-production variant for the Turbulence.jl SST model
(`step_sst!(...; production=:kato_launder)`). **A new feature, not a fix**
for the BFS over-prediction — the characterization below shows why it is
*not* a drop-in cure for that case.

## What it is

Standard SST production is `P_k = ν_t·S²` (S = strain-rate magnitude).
The Kato–Launder (1993) modification uses

```
P_k = ν_t · S · Ω        (Ω = vorticity magnitude)
```

Since `S = Ω` in pure shear, the two are identical in a shear layer; they
diverge only where strain dominates vorticity (stagnation, impingement,
reattachment), where Kato–Launder suppresses the spurious production that
plagues `S²` near stagnation points. It is the classic remedy for the
turbulence-kinetic-energy over-production at leading-edge stagnation in
standard two-equation models.

## Validation

- **Pure-shear / channel: no-op (as designed).** The SST channel at
  Re_τ = 395 is *identical* with Kato–Launder: u⁺ mean error 1.7 %,
  centreline 19.26 — bit-for-bit the same as standard SST, confirming
  `S = Ω` in the attached log layer. A Layer-1 unit test checks the two
  forms coincide in pure shear and differ in irrotational strain.
- The implementation reuses the existing strain and vorticity helpers;
  `Ω` is added to the SST per-cell closure.

## Backward-facing step: over-suppresses → unsteady (the interesting part)

On the clean ER=2 BFS at Re_H = 5000 (standard SST gives a steady
x_r/H = 9.90; OpenFOAM kΩSST 8.35), Kato–Launder does **not** settle to a
cleaner steady reattachment. Instead it drives the solution mildly
**unsteady**:

| step  | ν_t/ν max | x_r/H |
|------:|----------:|------:|
| 2000  | 0.0       | 9.35  |
| 8000  | 9.2       | 8.45  |
| 10000 | 0.0       | 9.90  |
| 12000 | 8.0       | 8.55  |
| 16000 | 6.7       | 8.60  |
| 20000 | 0.0       | 8.75  |
| 24000 | 0.0       | 8.95  |

ν_t periodically **collapses to ~0** and the reattachment oscillates
between 8.45 and 9.90 (|u|max rises to ~1.5 from the steady 1.16). The
mechanism: large parts of the BFS downstream/recovery flow are
strain-dominated (S > Ω), so `S·Ω ≪ S²` kills production broadly; ν_t
decays, the separated shear layer loses turbulent damping and sheds, then
production briefly regenerates — a limit-cycle, not a steady RANS.

So the snapshot x_r/H ≈ 8.5 that appears at the troughs is *not* a genuine
improvement toward OF — it is one phase of an oscillation whose
time-mean sits around ~9.0–9.5, no better than standard SST, and the
solution is no longer a clean steady state.

## Verdict

Kato–Launder is implemented, unit-tested, and **safe on attached flows**
(channel unchanged), available as an opt-in `production=:kato_launder`.
For *this* separated BFS it is **not** a cure for the reattachment
over-prediction: it over-suppresses ν_t and destabilises the RANS into
shedding rather than shortening the steady bubble. Its proper use case is
**stagnation-region over-production** (leading edges, impinging jets,
blade noses), where `S ≫ Ω` and the standard form spuriously generates
turbulence — not separated-shear-layer reattachment.

This reinforces the BFS conclusion (`RESULTS-bfs-sst.md`): the +19 %
reattachment over-prediction is a property of the BDIM substrate's
separated shear layer, and is not closed by resolution, advection scheme,
or production form. A sharper immersed-corner/shear-layer treatment is the
real lever.

## See also

- `Turbulence.jl/src/rans.jl` — `step_sst!(...; production=:kato_launder)`.
- [`RESULTS-bfs-sst.md`](RESULTS-bfs-sst.md), [`RESULTS-channel-sst.md`](RESULTS-channel-sst.md).
