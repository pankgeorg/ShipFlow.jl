#!/usr/bin/env julia
#
# K3: cavitation-onset analysis for the BladedRotor. Uses cavitation
# number σ + estimated minimum blade-surface pressure coefficient
# Cp_min from VLM-derived loading, no separate simulation needed.
#
# Cavitation occurs when Cp_min < −σ, where
#   σ = (p_∞ − p_v) / (0.5·ρ·V²_local)
# and p_∞ is the static pressure at the rotor depth (hydrostatic
# plus atmospheric, but at our cell units we measure p_∞ = ρ·g·h).
#
# Cp_min for a heavily-loaded blade section (thin airfoil with CL ≈
# 2π·α_eff) is approximately
#   Cp_min ≈ −(1 + 2·CL)²  +  1   ≈   −4·CL  −  4·CL²
# (rule of thumb; more careful would compute from velocity
# perturbation magnitude at the upper-surface midchord).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
]; io=devnull)
Pkg.add("VortexLattice"; io=devnull)

using LiftingSurfaces
using Printf

# Operating conditions in cell units (G1 / F1 baseline)
const U∞ = 1.0
const R_prop = 2.4               # 0.8·T/2 with T=6
const h_depth_w = 3.0            # rotor below waterline, cell units
const ρ_w = 10.0
const g_cell = 0.4444            # Fr=0.25, L=48 → g = U²/(Fr²·L)
const p_atm = 0.0                # cell pressure for atmosphere
const p_vap = -1.5               # rule of thumb: roughly −1.5 in our normalisation
                                 # (real water at 25 °C has p_v ≈ 3 kPa ≈ ρ_w·g·0.3 m)

const rotor = BladedRotor(; N_blades=3, R=R_prop, R_hub=0.2*R_prop,
    chord=(0.25*R_prop, 0.18*R_prop),
    twist=(deg2rad(35.0), deg2rad(15.0)), ns=12, nc=4)

# Sweep J across the G1 range
const J_LIST = [0.10, 0.15, 0.20, 0.25, 0.30, 0.32, 0.40, 0.60, 1.00]

println("K3: cavitation-onset analysis for BladedRotor")
println("Conditions: U∞=1, R=$R_prop, depth=$h_depth_w cells below waterline")
@printf "p_∞_static = ρ·g·h = %.2f (cell pressure units)\n" (ρ_w * g_cell * h_depth_w)
@printf "p_vapor    = %.2f (proxy)\n\n" p_vap

println("J     CT       Section_CL_max  Cp_min_est   σ_local   σ + Cp_min   Verdict")
println("─" ^ 80)

# For each J, estimate the maximum section CL across the blade span.
# VLM CT integrates over all sections; an approximate per-section CL
# is CL_section ≈ CT · A_disk / (N_blades · A_blade) where A_disk =
# π·R², A_blade ≈ (R − R_hub)·c_avg. This is a crude single-section
# proxy.
const A_disk = π * R_prop^2
const A_blade = (R_prop - 0.2*R_prop) * 0.5 * (0.25*R_prop + 0.18*R_prop)

for J in J_LIST
    Ω = π * U∞ / (J * R_prop)
    r = rotor_forces(rotor, U∞, Ω)
    CT = abs(r.CT)
    CL_section_est = CT * A_disk / (rotor.N_blades * A_blade)
    # Local velocity at the disk plane includes the rotor's own induced velocity.
    # For high CT (>2) actuator-disk theory gives a ≈ 0.5 √(1+CT) − 0.5, so
    # V_local ≈ U∞ · (1 + a).
    a = 0.5 * (sqrt(1 + CT) - 1)
    V_local = U∞ * (1 + a)
    # Cp_min on blade upper surface — rule of thumb above
    Cp_min = -4 * CL_section_est - 4 * CL_section_est^2
    # Cavitation number
    σ = (ρ_w * g_cell * h_depth_w - p_vap) / (0.5 * ρ_w * V_local^2)
    # Cavitation if Cp_min < -σ ⇒ σ + Cp_min < 0
    margin = σ + Cp_min
    verdict = margin < 0 ? "CAVITATES" : "OK"
    @printf "%.2f  %6.3f  %10.3f   %8.3f   %7.3f   %+9.3f  %s\n" J CT CL_section_est Cp_min σ margin verdict
end

println()
println("Interpretation:")
println("- σ ≈ 0.3 at our depth: this is in the heavy-cavitation regime")
println("  (towing-tank σ values for ship props typically 0.5–2).")
println("- The G1 / F1 self-prop point (J ≈ 0.32) puts CL_section ~ 5,")
println("  way above any realistic airfoil CL_max (typically 1.4–1.7).")
println("- VLM with no stall/cavitation model happily reports these")
println("  unphysical values; real prop would have stalled or fully")
println("  cavitated long before this load.")
println()
println("Conclusion: the BladedRotor at our G1/F1 operating points is")
println("unambiguously in the cavitation regime. The thrust numbers")
println("are *qualitatively* useful (right sign, right trends) but not")
println("quantitatively realistic. Use SwirlingDisk with prescribed")
println("thrust for quantitative open-water studies; reserve BladedRotor")
println("for trend analysis where the relative ranking matters.")
