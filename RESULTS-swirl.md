# Swirl vs uniform actuator disk — comparison behind Wigley hull

**Question.** Does the new `SwirlingDisk` (thrust + torque) materially
change the converged hull drag and free-surface response compared to
the uniform-thrust `ActuatorDisk`, for the same prescribed thrust?

**Setup.** Wigley L=48, B=10, T=6 at NX=128 × NY=48 × NZ=48. Disk
behind stern: R=1.8, w=1.5 at half-draft. Fr=0.30, Re=5000, ρ_w/ρ_a=10.
C_T = 0.60 (thrust=3.054). For the SWIRL case, Q/(T·R)=0.5
(torque=2.748). 180 timesteps. Drag stats taken from the last third.

Driver: [`scripts/swirl_vs_uniform_disk.jl`](scripts/swirl_vs_uniform_disk.jl).

## Result

|                              | UNIFORM         | SWIRL           | Δ%       |
|------------------------------|-----------------|-----------------|----------|
| Hull drag (mean ± σ)         | +63.749 ± 1.565 | +63.915 ± 1.287 | **+0.3%** |
| Wave RMS (post-stern)        | 0.1753          | 0.1800          | +2.7%    |
| Wave peak-to-peak            | 1.7733          | 1.8710          | +5.5%    |
| Thrust (prescribed)          | 3.054           | 3.054           |          |
| Torque (prescribed)          | 0.000           | 2.748           |          |

## Interpretation

For a moderate torque ratio Q/(T·R) = 0.5 (typical of merchant
propellers), **swirl has effectively no impact on the converged hull
drag** — the 0.3 % delta is well inside the σ ≈ 1.5 of the drag
fluctuations. The swirl energy mostly shows up in the **wake** rather
than the **hull boundary layer**: wave RMS rises 2.7 % and peak-to-peak
5.5 %. Physically reasonable — the tangential momentum injected at the
disk is shed downstream as helical wake, not back into the hull.

**Caveats.**
- Single Q/(T·R) value. Real propellers vary 0.2–0.8 across operating
  conditions; the effect on drag may grow at the high end.
- Single Re=5000. Boundary-layer/wake interaction strength scales with
  Re; the difference may be larger at ship-scale Re.
- Disk is centred behind the stern, axis +x. Asymmetric placement
  (offset from hull centreplane) would change the result.

## Conclusion

For the canonical stern-screw geometry at moderate Q/T, the uniform-thrust
ActuatorDisk is a defensible cheap approximation for hull-drag
prediction. SwirlingDisk's added cost is justified mainly when the
**wake structure** matters (cavitation, hull-propeller interaction
studies, free-surface signature).
