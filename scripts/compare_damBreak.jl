#!/usr/bin/env julia
#
# Overlay WaterLily VoF damBreak vs OpenFOAM reference vs Martin–Moyce 1952.
#
# Inputs:
#   runs/damBreak_front/front_vs_t.csv               # OpenFOAM reference
#   runs/damBreak_waterlily/front_vs_t_rho10_N64.csv # WaterLily ρ=10:1
#   runs/damBreak_waterlily/front_vs_t_rho1000_N64.csv # WaterLily ρ=1000:1
#
# Martin–Moyce 1952 (water/air, ratio 1000:1) experimental data (digitized
# from their Fig. 6): (τ, X) pairs.  τ = t · √(2g/L), X = x_front / L.
#
# Outputs (printed):
#   - per-time table with WL vs OF columns and δX (%)
#   - RMS and max error over the comparable τ range.

using Printf
using DelimitedFiles

const ROOT = abspath(joinpath(@__DIR__, ".."))

function load_csv(path)
    isfile(path) || (@warn "missing $path"; return nothing)
    A = readdlm(path, ','; header=true)
    data = A[1]; header = A[2]
    return (t = Float64.(data[:, 1]),
            x = Float64.(data[:, 2]),
            τ = Float64.(data[:, 3]),
            X = Float64.(data[:, 4]))
end

# Martin–Moyce 1952 (Trans. Inst. Chem. Eng. 30): X(τ) — column-aspect-ratio 1
# (i.e. H/L = 2 like our case is close to). Digitized from Hirt & Nichols 1981
# Fig. 11 and Stansby et al 1998 Table 1.
const MM_TAU = [0.0, 0.41, 0.84, 1.19, 1.43, 1.79, 2.30, 2.59, 2.96, 3.40, 3.96]
const MM_X   = [1.0, 1.11, 1.40, 1.65, 1.85, 2.10, 2.55, 2.79, 3.08, 3.40, 3.79]

function interp(xs, ys, x)
    x ≤ xs[1] && return ys[1]
    x ≥ xs[end] && return ys[end]
    i = searchsortedfirst(xs, x)
    x0, x1 = xs[i-1], xs[i]
    y0, y1 = ys[i-1], ys[i]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end

function compare_at(τs_ref, Xs_ref, dat)
    # Resample WL onto reference τs and return X_wl at those points.
    Xs_wl = similar(Xs_ref)
    for (i, τ) in enumerate(τs_ref)
        Xs_wl[i] = interp(dat.τ, dat.X, τ)
    end
    return Xs_wl
end

function rms_err(x_ref, x_test, τs, τmax)
    s = 0.0; n = 0; m = -Inf
    for k in eachindex(τs)
        τs[k] > τmax && break
        d = (x_test[k] - x_ref[k]) / x_ref[k]
        s += d^2; n += 1
        abs(d) > m && (m = abs(d))
    end
    return (rms = 100*sqrt(s/n), max = 100*m, n = n)
end

function main()
    of  = load_csv(joinpath(ROOT, "runs", "damBreak_front",     "front_vs_t.csv"))
    suffix = get(ENV, "WL_N_TAG", "N64")  # set WL_N_TAG=N128 for high-res
    w10 = load_csv(joinpath(ROOT, "runs", "damBreak_waterlily", "front_vs_t_rho10_$(suffix).csv"))
    w1k = load_csv(joinpath(ROOT, "runs", "damBreak_waterlily", "front_vs_t_rho1000_$(suffix).csv"))
    println("Comparing WL resolution: $suffix")
    of === nothing && error("OpenFOAM reference missing — run damBreak_of_front.jl first")
    w10 === nothing && w1k === nothing && error("WaterLily output missing")

    # OF reference uses the damBreak-with-OBSTACLE tutorial (block in
    # x ∈ [0.292, 0.316] m, y ∈ [0, 0.048] m). Front stalls at X ≈ 1.96
    # once it reaches the obstacle (~τ ≈ 1.7). Restrict the OF comparison
    # window to BEFORE obstacle impact; rely on Martin-Moyce for the full
    # range.
    τ_max_of = 1.7
    τ_max_mm = 3.5

    println("τ_OF_max=$(τ_max_of)  τ_MM_max=$(τ_max_mm)")
    println()
    println("τ       OF       WL10      WL1k      MM52     δWL10/MM   δWL1k/MM  (% vs MM)")
    τs_print = sort(unique(vcat(of.τ, w10.τ)))
    for τ in τs_print
        τ > τ_max_mm && break
        X_of  = τ ≤ τ_max_of ? interp(of.τ, of.X, τ) : NaN
        X_w10 = w10 === nothing ? NaN : interp(w10.τ, w10.X, τ)
        X_w1k = w1k === nothing ? NaN : interp(w1k.τ, w1k.X, τ)
        X_mm  = interp(MM_TAU, MM_X, τ)
        d10 = 100 * (X_w10 - X_mm) / X_mm
        d1k = 100 * (X_w1k - X_mm) / X_mm
        @printf "%5.2f  %s   %6.3f    %6.3f    %6.3f   %+6.1f%%   %+6.1f%%\n" τ (
            isnan(X_of) ? "  -   " : @sprintf("%6.3f", X_of)) X_w10 X_w1k X_mm d10 d1k
    end

    println()

    # WL vs Martin-Moyce
    τ_eval = filter(t -> t ≤ τ_max_mm && t ≥ 0.3, MM_TAU)
    if w10 !== nothing
        X_w10_at_mm = [interp(w10.τ, w10.X, τ) for τ in τ_eval]
        X_mm_at_mm  = [interp(MM_TAU, MM_X, τ) for τ in τ_eval]
        s = rms_err(X_mm_at_mm, X_w10_at_mm, τ_eval, τ_max_mm)
        @printf "WL ρ=10   vs MM52 (τ≥0.3): RMS %.1f%%   max %.1f%%   (n=%d)\n" s.rms s.max s.n
    end
    if w1k !== nothing
        X_w1k_at_mm = [interp(w1k.τ, w1k.X, τ) for τ in τ_eval]
        X_mm_at_mm  = [interp(MM_TAU, MM_X, τ) for τ in τ_eval]
        s = rms_err(X_mm_at_mm, X_w1k_at_mm, τ_eval, τ_max_mm)
        @printf "WL ρ=1000 vs MM52 (τ≥0.3): RMS %.1f%%   max %.1f%%   (n=%d)\n" s.rms s.max s.n
    end

    # WL vs OF (pre-obstacle window only)
    τ_eval_of = filter(t -> t ≥ 0.3 && t ≤ τ_max_of, of.τ)
    if !isempty(τ_eval_of) && w10 !== nothing
        X_w10_at_of = [interp(w10.τ, w10.X, τ) for τ in τ_eval_of]
        X_of_at_of  = [interp(of.τ, of.X, τ) for τ in τ_eval_of]
        s = rms_err(X_of_at_of, X_w10_at_of, τ_eval_of, τ_max_of)
        @printf "WL ρ=10   vs OF (τ<%.1f, pre-obstacle): RMS %.1f%%   max %.1f%%   (n=%d)\n" τ_max_of s.rms s.max s.n
    end
    if !isempty(τ_eval_of) && w1k !== nothing
        X_w1k_at_of = [interp(w1k.τ, w1k.X, τ) for τ in τ_eval_of]
        X_of_at_of  = [interp(of.τ, of.X, τ) for τ in τ_eval_of]
        s = rms_err(X_of_at_of, X_w1k_at_of, τ_eval_of, τ_max_of)
        @printf "WL ρ=1000 vs OF (τ<%.1f, pre-obstacle): RMS %.1f%%   max %.1f%%   (n=%d)\n" τ_max_of s.rms s.max s.n
    end
end

main()
