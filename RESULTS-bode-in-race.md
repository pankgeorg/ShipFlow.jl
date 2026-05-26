# Rudder Bode plot inside the rotor race (J2)

Cross of F2 (rudder placed at `prop_xc + 0.5·R`, inside the race,
two-way coupling on) and I4 (sinusoidal δ(t) frequency sweep).
Quantifies how the +300 % CL amplification (F2's static finding)
shows up in the dynamic transfer function.

Driver: [`scripts/rudder_in_race_freq_sweep.jl`](scripts/rudder_in_race_freq_sweep.jl).

## Setup

- Same Wigley + rotor + rudder geometry as F2.
- **rud_xc = prop_xc + 0.5·R** (inside the race, NOT the G2/I4
  placement at 2·R).
- `inflow=trilinear_inflow(flow.u)` on the rudder so it samples
  the local accelerated flow.
- δ(t) = 10° · sin(2π · t/T) for periods `T ∈ {80, 40, 20}` steps.
- 3 cycles per frequency; 1 cycle burn-in discarded.

## Result

| Period | f (1/cell-time) | \|A_CL\| | \|A_F\| | Gain \|A_F/A_CL\| | Phase F vs CL (°) |
|--------|-----------------|----------|---------|---------------------|----------------------|
| 80     | 0.050           | **1.218**| 0.214   | 0.176              | −5.0                |
| 40     | 0.100           | **1.215**| 0.329   | 0.271              | −5.4                |
| 20     | 0.200           | **1.214**| 0.347   | 0.286              | +3.2                |

## Comparison to I4 (out-of-race baseline)

| Quantity         | I4 (out of race) | J2 (in race) | Ratio |
|------------------|-----------------|--------------|-------|
| \|A_CL\| (f=0.05)| 0.303           | **1.218**    | **4.0×** |
| \|A_F\| (f=0.05) | 0.109           | 0.214        | 1.96×  |
| Gain (f=0.05)    | 0.359           | 0.176        | 0.49×  |
| Phase (f=0.05)   | −187°           | **−5°**      | flip   |

## Interpretation

- **CL amplification persists across frequencies.** The static F2
  finding (+300 % CL at δ=10°) is reproduced in the AC amplitude:
  |A_CL| jumps from 0.30 (out of race) to 1.21 (in race) — exactly
  4×, independent of frequency. The race acts as a constant
  multiplier on rudder CL, not a frequency-dependent filter.
- **Hull side-force amplification is smaller.** |A_F| grows only
  2×, not 4×. Same effect as in the I3 diagnostic: most of the
  amplified rudder force is convected into the race, never reaching
  the hull's pressure field. The gain ratio |A_F/A_CL| therefore
  *drops* to roughly half its out-of-race value.
- **Phase flips by ~180°.** This is the striking new finding.
  Out-of-race, the rudder side force lifts the hull pressure field
  in the opposite direction (180° out of phase — Newton's third
  law). In-race, the rudder lifts WITH the rotor jet to push the
  fluid sideways, and the hull's pressure response is in phase
  with the rudder CL.
- **Implication for manoeuvring**: a rudder inside the race
  produces a hull-side-force in the SAME sense as a rudder outside
  the race produces (port rudder ⇒ port side force on hull), but
  the *mechanism* differs. Out of race, the hull pressure field
  reacts directly to the rudder. In race, the rotor jet pushes the
  fluid sideways and the hull is dragged with it. Both produce
  yaw, but the lag characteristics differ — important for fast
  manoeuvres where the time history matters.

## Caveats

- 3 cycles per frequency is the minimum for Fourier projection;
  the highest f (period=20) has only 40 useful steps after burn-in,
  marginal for Fourier accuracy.
- Single δ_max = 10°. The system is quasi-linear at this amplitude;
  larger δ_max would test the nonlinear regime where rudder stall
  and rotor-race choking matter.
- No comparison run with both rudders simultaneously (e.g., a
  ducted-propeller geometry) — there's no clear physics motivation
  for this, but worth noting we sampled only one rudder position.

## See also

- `RESULTS-twoway-amplification.md` — F2: static in-race CL amplification.
- `RESULTS-rudder-bode.md` — I4: out-of-race frequency sweep.
- `RESULTS-twoway-F2-diagnostic.md` — I3: spatial diagnostic
  showing the rudder force is convected into the race.
- `runs/rudder_in_race_freq/bode_in_race.png` — the plot.
