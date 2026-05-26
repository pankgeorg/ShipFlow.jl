# Self-propulsion J vs analytical models (I10)

How does our measured `J_self` for Wigley (0.315) and Containership
(0.176) compare to what classical actuator-disk + thin-ship-wave
theory would predict?

Short answer: **it can't be cleanly compared**, and that's a useful
finding. Below the reasoning.

## Actuator-disk momentum theory at our operating point

Momentum theory relates disk thrust to inflow:

```
CT_disk_loading = T / (0.5 · ρ · V² · A_disk) = 4 · a · (1 + a)
```

where `a` is the axial induction factor. For a real propeller at
near-optimal advance ratio, CT_disk_loading is ~0.5–1.5
(a ~ 0.1–0.3).

Our measurements at the respective self-prop points:

| Hull          | J_self | VLM CT |
|---------------|--------|--------|
| Wigley        | 0.315  | ~3.0   |
| Containership | 0.176  | ~44    |

Solving `4a(1+a) = CT` gives `a ≈ 0.50` for Wigley (just at the
edge of momentum-theory validity), and `a ≈ 2.9` for Containership
(well outside any classical regime).

**Conclusion**: the VLM is operating in a heavily-loaded regime
where classical actuator-disk theory does not predict thrust at all
— momentum theory assumes the slipstream contraction is small, and
at CT > 4 the slipstream would invert (suction face). The VLM
itself has no stall model, so it confidently reports unphysically
large CT at low J. We use this as a **methodology test only**, not
a quantitative prediction.

## Thin-ship-wave theory for the resistance

Wigley wave resistance from Bai 1979 at Fr = 0.25:
`Rw / (0.5 ρ U² ∇^(2/3)) ≈ 5×10⁻³`.

In cell units with `V = 1`, `ρ_w = 10`, `∇ = 4·L·B·T/9 ≈ 1280`:
`Rw ≈ 5e-3 · 0.5 · 10 · 1 · 1280^(2/3) ≈ 3.0`.

Our measured drag at J_self (rotor providing the thrust) is `D ≈
144`. So **viscous drag dominates by ~50×** in our setup — Re = 1000
is laminar Blasius-flat-plate regime, not the Re ≈ 10⁹ of a real
ship. The wave-resistance signal is buried in viscous drag.

For a meaningful wave-resistance comparison we'd need:
1. Higher Re (10⁴–10⁵ with WALE LES), or
2. A no-wave deep-water reference run to subtract viscous.

I5 (`RESULTS-Fr-scan-current.md`) gives the Fr-dependent total
resistance curve from which the wave portion can be estimated
empirically.

## What J_self does test cleanly

- **Methodology**: the matching procedure (sweep J, find T − D = 0)
  is consistent. The Wigley and Containership values trend in the
  physically expected direction (higher Cb → lower J_self).
- **Empirical scaling**: J_self ≈ 0.453 − 0.367 · Cb across the
  Containership family (I1) is a usable preliminary-sizing rule
  *for our setup* — different Re/Fr would give different
  coefficients.

## What J_self does not test

- **Absolute thrust at a given J** — VLM overestimates CT in the
  heavily-loaded regime by ~10× compared to towing-tank data.
- **Real-world J_self values** for the same hull family at full
  scale.

## See also

- `RESULTS-selfprop-VLM.md` — Wigley J_self.
- `RESULTS-selfprop-containership.md` — Containership J_self.
- `RESULTS-cb-vs-Jself.md` — sweep (I1).
- `RESULTS-Fr-scan-current.md` — wave resistance vs Fr (I5).
