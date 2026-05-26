# Rudder frequency sweep — manoeuvring transfer function (I4)

Drives `δ(t) = δ_max · sin(2π · t / T)` for four periods `T ∈ {80, 40,
20, 12}` steps and reports the single-frequency Fourier projection of
`F_hull,y(t)` relative to `CL_rud(t)`. Mentioned as future work in
`RESULTS-rudder-trial.md`.

Driver: [`scripts/rudder_frequency_sweep.jl`](scripts/rudder_frequency_sweep.jl).

## Setup

- Same Wigley + rotor + rudder geometry as the G2 / F2 trials.
- Rudder placed at `prop_xc + 2·R` (G2 baseline, OUT of the race —
  separate study from F2's in-race placement).
- δ_max = 10°, dt = 0.25 cell-time-unit.
- 3 full cycles per frequency; first cycle discarded as transient.

## Result

| Period (steps) | f (1/cell-time) | \|Aδ\| (°) | \|A_CL\| | \|A_F\| | Gain \|A_F/A_CL\| | Phase F vs CL (°) |
|----------------|-----------------|------------|----------|---------|--------------------|----------------------|
| 80             | 0.050           | 10.00      | 0.303    | 0.109   | **0.359**          | −187                |
| 40             | 0.100           | 10.00      | 0.303    | 0.029   | 0.096              | −192                |
| 20             | 0.200           | 10.00      | 0.303    | 0.026   | 0.085              | −177                |
| 12             | 0.333           | 10.00      | 0.303    | 0.045   | 0.147              | −22                 |

Bode plot: `runs/rudder_freq_sweep/bode.png`.
Time series: `runs/rudder_freq_sweep/timeseries.png`.

## Interpretation

- **Quasi-steady rudder VLM**: The CL_rud amplitude is independent of
  frequency (≈ 0.303 at every f), confirming that `rudder_forces`
  treats each call as a static VLM solve. Dynamic lift coefficients
  (Theodorsen unsteady aero) would require time-history aware
  modelling — out of scope here.
- **Low-frequency manoeuvring gain ≈ 0.36** at f = 0.050
  (cell-time)⁻¹. This matches the order of the G2 ramp result
  (where `|F_hull|/CL ≈ 5.75` in different units — the
  normalisations differ).
- **Out-of-phase response at low frequencies** (phase ≈ −180°):
  F_hull,y is *opposite in sign* to CL_rud. The rudder lifts to
  starboard → the integrated hull pressure reaction is to port. This
  is the standard Newton's-third-law response with no significant
  lag at these low frequencies (180° in the convention `F·cos(ωt+φ)`
  is just a sign flip).
- **Phase wraps to ≈ 0° at f = 0.333**: this is the highest frequency
  tested with only 36 total steps × 0.25 = 9 cell-time-units total,
  ≈ 3 wavelengths. The Fourier projection is under-sampled and the
  result is **unreliable**. Treat the rightmost point in the Bode
  plot as a caveat marker, not a number.
- **Gain rolls off above f ≈ 0.1**: from 0.36 (f=0.05) to ≈ 0.09
  (f=0.1, 0.2) — about 12 dB attenuation per decade. A
  1st-order low-pass roll-off with corner ≈ 0.07 cell-time-units⁻¹
  is consistent with the data.

## Caveats

- 3 cycles per frequency is the minimum that supports Fourier
  projection. 5+ would tighten the gain numbers, especially at the
  highest frequencies.
- The 80×64×32 grid resolves wave numbers down to ≈ 2 cells; the
  rudder's CL amplitude (≈ 0.3) deposits pressure perturbations that
  may be at the grid limit. Doubling the grid would test the
  resolution sensitivity.
- The rudder is placed *outside* the rotor race (G2 geometry); a
  separate sweep with rudder INSIDE the race (F2 geometry) would
  show how the rotor-jet convection modifies the transfer function.

## See also

- `RESULTS-rudder-trial.md` — G2: ramp response (step-like input).
- `RESULTS-twoway-amplification.md` — F2: amplification with rudder
  inside the race.
- `scripts/rudder_frequency_sweep.jl` — driver.
- `runs/rudder_freq_sweep/bode.png` — gain + phase plot.
- `runs/rudder_freq_sweep/timeseries.png` — δ, CL, F per frequency.
