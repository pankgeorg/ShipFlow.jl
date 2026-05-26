# Containership Cb sweep (I1)

Generalises F1 to a one-parameter family: vary `par_frac` over
`{0.30, 0.50, 0.70, 0.85}` and find the self-propulsion advance ratio
`J_self` at each. Same VLM rotor, same Fr, same Re, same procedure.

Driver: [`scripts/containership_cb_sweep.jl`](scripts/containership_cb_sweep.jl).

## Setup

- Grid 96 × 48 × 48, hull L = 48, B = 10, T = 6.
- Fr = 0.25, Re = 1000, ρ_w/ρ_a = 10.
- 100 steps per (par_frac, J), drag averaged over last 25 %.
- BladedRotor at the stern: 3 blades, R = 0.8·T/2, R_hub = 0.2·R,
  chord (0.25·R → 0.18·R), twist (35° → 15°), 12 × 4 panels.

## Result

| par_frac | Cb    | J_self  |
|----------|-------|---------|
| 0.30     | 0.65  | 0.2106  |
| 0.50     | 0.75  | 0.1762  |
| 0.70     | 0.85  | 0.1448  |
| 0.85     | 0.93  | 0.1106  |

**Monotonic, near-linear** — J_self decreases by roughly 0.36 per unit
increase in Cb. An empirical fit:

```
J_self ≈ 0.452 − 0.367 · Cb
```

A naive heuristic — fuller hull → more drag → propeller works harder
→ lower J — is borne out cleanly.

## Interpretation

- **Hull drag scales steeply with Cb.** Going from Cb 0.65 → 0.93
  (par_frac 0.30 → 0.85), the drag at fixed J = 0.20 goes from ~291
  to ~705 (cell units) — a 2.4× increase for a 1.43× increase in Cb.
  Most of that is wave-making (this is at Fr = 0.25, near the
  steepest part of the Wigley resistance curve); the rest is form
  drag from the bluffer ends.
- **CT at J_self grows accordingly.** At par_frac=0.85, CT ≈ 100+
  at the self-prop point — VLM is well outside its linear range.
  Real propellers cannot operate there (cavitation, stall);
  practical containerships are limited to Cb ≲ 0.85 and run with
  much larger propellers than the 0.8·T/2 we used.
- **The empirical fit is useful for preliminary sizing.** Given a
  target Cb and a fixed rotor geometry, this relation gives a
  first-guess J for the resistance/self-propulsion match before
  running a full sweep.

## Caveats

- Single rotor geometry (constant R, hub fraction, twist). A more
  realistic study would scale R with hull beam or draught.
- No two-way coupling on the rotor (see [I2] below — work-in-progress)
  so the rotor sees freestream not hull wake. Two-way would lower CT
  by the wake fraction and shift J_self slightly upward.
- The 4-point Cb grid is coarse; a finer sweep would let us bound
  the fit's quadratic term.
- Fr = 0.25 only — at higher Fr the wave-resistance hump moves and
  the slope of J_self(Cb) would change.

## See also

- `RESULTS-selfprop-containership.md` — F1, par_frac = 0.50 detail.
- `RESULTS-selfprop-VLM.md` — original Wigley G1 (different hull
  family, but for reference J_self = 0.315 there).
- `runs/cb_sweep/scans.csv` — the underlying per-J data.
