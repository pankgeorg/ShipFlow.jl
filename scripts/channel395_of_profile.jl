#!/usr/bin/env julia
#
# OpenFOAM channel395 — extract the time-averaged + (x,z)-averaged
# streamwise velocity profile u_x(y) for cross-validation against
# Turbulence.jl Smagorinsky.

using ShipFlow.Harness: read_volVectorField
using Printf

const CASE   = abspath(joinpath(@__DIR__, "..", "runs", "channel395"))
const OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "channel395_of_profile"))
mkpath(OUTDIR)

# Time-averaging window in physical seconds: average over the second half
# of the OF run so transients have decayed.
const T_AVG_START = 500.0
const T_AVG_END   = Inf

# Channel geometry (from blockMesh): y ∈ [0, 2], δ = 1.
const DELTA = 1.0
const UBAR  = 0.1335    # target bulk velocity

function time_dirs(case)
    times = String[]
    for d in readdir(case)
        isdir(joinpath(case, d)) || continue
        try
            v = parse(Float64, d); v ≥ 0 && push!(times, d)
        catch; end
    end
    return sort(times; by = s -> parse(Float64, s))
end

function average_U(case; t_start, t_end)
    times = time_dirs(case)
    selected = filter(t -> t_start ≤ parse(Float64, t) ≤ t_end, times)
    isempty(selected) && error("no time directories in [$t_start, $t_end]")
    @info "Averaging U over $(length(selected)) snapshots: $(first(selected)) → $(last(selected))"
    accum = nothing
    for t in selected
        n, vals = read_volVectorField(joinpath(case, t, "U"))
        if accum === nothing
            accum = deepcopy(vals)
        else
            for i in eachindex(accum)
                accum[i] = accum[i] .+ vals[i]
            end
        end
    end
    for i in eachindex(accum)
        accum[i] = accum[i] ./ length(selected)
    end
    return accum
end

function read_cell_centers(case)
    path = joinpath(case, "0", "C")
    n, vecs = read_volVectorField(path)
    return [v[1] for v in vecs], [v[2] for v in vecs], [v[3] for v in vecs]
end

function bin_average_u_vs_y(U_avg, ys, n_bins)
    # Use unique y-rows directly (the OF mesh is structured: every cell at the
    # same y has the same coordinate to within rounding).
    rounded_y = round.(ys; digits=5)
    unique_ys = sort(unique(rounded_y))
    sums    = zeros(Float64, length(unique_ys))
    counts  = zeros(Int,     length(unique_ys))
    idx_map = Dict(y => i for (i, y) in pairs(unique_ys))
    for (k, y) in pairs(rounded_y)
        i = idx_map[y]
        sums[i]   += U_avg[k][1]
        counts[i] += 1
    end
    return unique_ys, sums ./ counts, counts
end

function main()
    @info "Reading OF channel395 mesh + averaging U…"
    xs, ys, zs = read_cell_centers(CASE)
    U_avg = average_U(CASE; t_start=T_AVG_START, t_end=T_AVG_END)
    @info "Cells: $(length(ys)), y range: $(minimum(ys)) … $(maximum(ys))"

    y_c, u_c, counts = bin_average_u_vs_y(U_avg, ys, 0)
    n_bins = length(y_c)
    @info "Resolved $n_bins distinct y-rows from the mesh"

    csv = joinpath(OUTDIR, "ux_profile.csv")
    open(csv, "w") do io
        println(io, "y_over_delta,u_over_Ubar,n_cells")
        for k in eachindex(y_c)
            y_over_delta = (y_c[k] - DELTA) / DELTA
            u_over_Ubar  = u_c[k] / UBAR
            @printf io "%.4f,%.6f,%d\n" y_over_delta u_over_Ubar counts[k]
        end
    end
    @info "wrote $csv"

    # u_τ estimate: linear wall shear from the first cell (which sits at y_c[1])
    u_first = u_c[1]
    NU = 2e-5
    u_tau = sqrt(NU * u_first / y_c[1])
    Re_tau_est = u_tau * DELTA / NU
    @printf("\n=== OpenFOAM channel395 summary ===\n")
    @printf("  ⟨u⟩ from bin 1 (near wall) = %.4f (u/Ubar = %.4f)\n", u_first, u_first/UBAR)
    @printf("  ⟨u⟩ centre bin              = %.4f (u/Ubar = %.4f)\n",
            u_c[div(n_bins, 2)], u_c[div(n_bins, 2)]/UBAR)
    @printf("  u_τ estimate                = %.4e\n", u_tau)
    @printf("  Re_τ estimate               = %.1f (target 395)\n", Re_tau_est)
    @printf("===================================\n")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
