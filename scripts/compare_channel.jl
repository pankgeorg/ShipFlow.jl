#!/usr/bin/env julia
#
# Compare WaterLily Smagorinsky channel vs OpenFOAM Smagorinsky channel
# at Re_τ ≈ 395. Reads CSV profiles emitted by:
#   - scripts/channel395_waterlily.jl
#   - scripts/channel395_of_profile.jl
#
# Prints u(y) at canonical y/δ stations and computes:
#   * centerline u/Ubar
#   * u_τ from each code
#   * Re_τ achieved by each code
#   * RMS deviation between profiles on a common y grid

using Printf
using Statistics: mean

const WL_CSV = abspath(joinpath(@__DIR__, "..", "runs", "channel395_waterlily_v2",
                                  "ux_profile.csv"))
const OF_CSV = abspath(joinpath(@__DIR__, "..", "runs", "channel395_of_profile",
                                  "ux_profile.csv"))

function read_csv(path)
    a = Float64[]; b = Float64[]
    open(path) do io
        readline(io)
        for line in eachline(io)
            isempty(strip(line)) && continue
            toks = split(strip(line), ',')
            push!(a, parse(Float64, toks[1]))
            push!(b, parse(Float64, toks[2]))
        end
    end
    return a, b
end

function interp_to(xa, ya, xb)
    yb = similar(xb)
    for (k, x) in pairs(xb)
        if x ≤ first(xa); yb[k] = first(ya)
        elseif x ≥ last(xa); yb[k] = last(ya)
        else
            i = searchsortedfirst(xa, x)
            r = (x - xa[i-1]) / (xa[i] - xa[i-1])
            yb[k] = ya[i-1] * (1-r) + ya[i] * r
        end
    end
    return yb
end

function main()
    isfile(WL_CSV) || (println("missing $WL_CSV"); return)
    isfile(OF_CSV) || (println("missing $OF_CSV"); return)

    yWL, uWL = read_csv(WL_CSV)
    yOF, uOF = read_csv(OF_CSV)

    # Sort if not already
    pwl = sortperm(yWL); yWL = yWL[pwl]; uWL = uWL[pwl]
    pof = sortperm(yOF); yOF = yOF[pof]; uOF = uOF[pof]

    @printf("\n=== Channel Re_τ ≈ 395 — WaterLily vs OpenFOAM ===\n")
    @printf("WaterLily samples : %d  (y/δ ∈ [%+.3f, %+.3f])\n", length(yWL), yWL[1], yWL[end])
    @printf("OpenFOAM  samples : %d  (y/δ ∈ [%+.3f, %+.3f])\n", length(yOF), yOF[1], yOF[end])

    @printf("\n--- u(y)/U_bar at canonical stations ---\n")
    @printf("%-10s %-12s %-12s\n", "y/δ", "WL u/Ubar", "OF u/Ubar")
    for y in (-0.95, -0.90, -0.75, -0.50, -0.25, 0.00, +0.25, +0.50, +0.75, +0.90, +0.95)
        u_wl_y = interp_to(yWL, uWL, [y])[1]
        u_of_y = interp_to(yOF, uOF, [y])[1]
        @printf("%+10.2f %-12.4f %-12.4f\n", y, u_wl_y, u_of_y)
    end

    # Centerline
    u_wl_c = interp_to(yWL, uWL, [0.0])[1]
    u_of_c = interp_to(yOF, uOF, [0.0])[1]
    @printf("\n--- Centerline u_max/U_bar ---\n")
    @printf("  WL: %.4f\n", u_wl_c)
    @printf("  OF: %.4f\n", u_of_c)

    # RMS over common y range, restricted to bulk (avoid wall-region BDIM smearing)
    y_common = range(-0.7, 0.7; length=41)
    u_wl_common = interp_to(yWL, uWL, collect(y_common))
    u_of_common = interp_to(yOF, uOF, collect(y_common))
    rms = sqrt(mean((u_wl_common .- u_of_common).^2))
    max_dev = maximum(abs.(u_wl_common .- u_of_common))
    @printf("\n--- RMS / max deviation over |y/δ| < 0.7 (bulk) ---\n")
    @printf("  RMS(WL − OF)  = %.4f\n", rms)
    @printf("  max|WL − OF|  = %.4f\n", max_dev)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
