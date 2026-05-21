#!/usr/bin/env julia
#
# Compare time-averaged spatial profiles between WaterLily and OpenFOAM
# for the cylinder Re=100 case. Reads CSVs emitted by:
#   - scripts/cylinder_Re100_profiles.jl  (WL)
#   - scripts/of_sample_profiles.jl       (OF)
#
# Prints a tabular comparison and computes L₂ errors between the two
# datasets at the published spatial stations.

using Printf
using Statistics: mean

const WL_DIR = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_waterlily_profiles"))
const OF_DIR = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_of_profiles"))

function read_two_col(path)
    a = Float64[]; b = Float64[]
    open(path) do io
        readline(io)              # header
        for line in eachline(io)
            isempty(strip(line)) && continue
            x, y = split(strip(line), ',')
            push!(a, parse(Float64, x))
            push!(b, parse(Float64, y))
        end
    end
    return a, b
end

# Interpolate (xa, ya) onto coord vector xb, returning ya values at xb.
function interp_to(xa, ya, xb)
    yb = similar(xb)
    for (k, x) in pairs(xb)
        if x ≤ first(xa)
            yb[k] = first(ya)
        elseif x ≥ last(xa)
            yb[k] = last(ya)
        else
            i = searchsortedfirst(xa, x)
            r = (x - xa[i-1]) / (xa[i] - xa[i-1])
            yb[k] = ya[i-1] * (1 - r) + ya[i] * r
        end
    end
    return yb
end

function l2_error(y_ref, y_test)
    n = length(y_ref)
    sqrt(sum((y_ref .- y_test).^2) / n)
end

function compare_line(label, wl_path, of_path; tol)
    isfile(wl_path) || (return @warn "missing $wl_path")
    isfile(of_path) || (return @warn "missing $of_path")
    x_wl, y_wl = read_two_col(wl_path)
    x_of, y_of = read_two_col(of_path)
    # Sample both on the WL grid (it's regular)
    y_of_on_wl = interp_to(x_of, y_of, x_wl)
    err = l2_error(y_of_on_wl, y_wl)
    max_dev = maximum(abs, y_of_on_wl .- y_wl)
    @printf("\n%s\n", label)
    @printf("  RMS(WL − OF)    = %.4f\n", err)
    @printf("  max |WL − OF|   = %.4f\n", max_dev)
    @printf("  WL min/max      = %.4f / %.4f\n", minimum(y_wl), maximum(y_wl))
    @printf("  OF min/max      = %.4f / %.4f\n", minimum(y_of_on_wl), maximum(y_of_on_wl))
    @printf("  pass (RMS < %g): %s\n", tol, err < tol ? "YES" : "NO")
    return err, max_dev
end

function main()
    println()
    println("==============================================================")
    println("Spatial profile comparison — WaterLily vs OpenFOAM (Re=100)")
    println("==============================================================")

    # Wake profile cross-stream
    for xD in (2, 5, 10)
        compare_line(
            "u_x(y) at x = $(xD)D",
            joinpath(WL_DIR, "ux_profile_x$(xD)D.csv"),
            joinpath(OF_DIR, "ux_profile_x$(xD)D.csv");
            tol = 0.10,
        )
    end

    # Wake centerline streamwise
    compare_line(
        "u_x(x) along wake centerline (y=0)",
        joinpath(WL_DIR, "ux_centerline.csv"),
        joinpath(OF_DIR, "ux_centerline.csv");
        tol = 0.10,
    )

    # Surface pressure coefficient
    compare_line(
        "Cp(θ) on cylinder surface (just outside)",
        joinpath(WL_DIR, "Cp_surface.csv"),
        joinpath(OF_DIR, "Cp_surface.csv");
        tol = 0.15,
    )

    # Side-by-side print of a few key stations
    println("\n--- u_x(y=0) at downstream stations ---")
    @printf("%-12s %-12s %-12s\n", "x/D", "WL u_x/U", "OF u_x/U")
    for xD in (1, 2, 3, 5, 7, 10)
        xc_path_wl = joinpath(WL_DIR, "ux_centerline.csv")
        xc_path_of = joinpath(OF_DIR, "ux_centerline.csv")
        x_wl, y_wl = read_two_col(xc_path_wl)
        x_of, y_of = read_two_col(xc_path_of)
        wl_val = interp_to([Float64(xD)], y_wl, x_wl)[1]   # invert: ya@x
        # Easier:
        i_wl = argmin(abs.(x_wl .- xD))
        i_of = argmin(abs.(x_of .- xD))
        @printf("%-12.1f %-12.4f %-12.4f\n", xD, y_wl[i_wl], y_of[i_of])
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
