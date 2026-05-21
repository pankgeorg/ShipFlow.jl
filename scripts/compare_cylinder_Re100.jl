#!/usr/bin/env julia
#
# Side-by-side comparison of the 2D cylinder Re=100 flow as simulated by
#   • OpenFOAM   (runs/cylinder_Re100/postProcessing/forceCoeffs1/0/forceCoeffs.dat)
#   • WaterLily  (runs/cylinder_waterlily/cylinder_waterlily.csv)
#
# Computes mean drag coefficient (Cd) over the second half of each
# simulation and the Strouhal number St via zero-crossing counting on
# the lift coefficient (Cl). Robust to noise; no FFTW dependency.
#
# Run:
#   julia --project=. scripts/compare_cylinder_Re100.jl

using Printf
using Statistics: mean

const OF_FORCES = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_Re100",
                                    "postProcessing", "forceCoeffs1", "0",
                                    "forceCoeffs.dat"))
const WL_CSV    = abspath(joinpath(@__DIR__, "..", "runs", "cylinder_waterlily",
                                    "cylinder_waterlily.csv"))

# OpenFOAM physical parameters from the case
const OF_D  = 0.001
const OF_U  = 1.5

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

# forceCoeffs.dat layout (OpenFOAM 11):
#   t   Cm   Cd   Cl   Cl(f)   Cl(r)
# Comment lines start with #.
function read_openfoam_forces(path)
    t = Float64[]; Cd = Float64[]; Cl = Float64[]
    for line in eachline(path)
        sline = strip(line)
        (isempty(sline) || startswith(sline, '#')) && continue
        toks = split(sline)
        length(toks) ≥ 4 || continue
        push!(t,  parse(Float64, toks[1]))
        push!(Cd, parse(Float64, toks[3]))
        push!(Cl, parse(Float64, toks[4]))
    end
    return (t = t, Cd = Cd, Cl = Cl)
end

# WaterLily CSV layout:
#   t,Cd,Cl
function read_waterlily_csv(path)
    t = Float64[]; Cd = Float64[]; Cl = Float64[]
    open(path) do io
        header = readline(io)
        for line in eachline(io)
            isempty(strip(line)) && continue
            a, b, c = split(strip(line), ',')
            push!(t,  parse(Float64, a))
            push!(Cd, parse(Float64, b))
            push!(Cl, parse(Float64, c))
        end
    end
    return (t = t, Cd = Cd, Cl = Cl)
end

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

# Mean of x over the second half of its range, returned with the index window.
function tail_mean(x::AbstractVector)
    n = length(x)
    half = max(1, n ÷ 2)
    return mean(@view x[half:end])
end

# Strouhal estimate via zero-crossings of Cl over the second half:
#   St = f * D / U
#
# The trace is high-pass-filtered by subtracting the mean. We count
# *up-going* zero-crossings only (one per period). With M crossings in
# time window Δt, frequency f = (M-1)/Δt → St = f * D / U.
#
# Returns (St, n_crossings). NaN if too few crossings.
function strouhal_from_lift(t::AbstractVector, Cl::AbstractVector;
                            D::Real, U::Real)
    n = length(t); half = max(1, n ÷ 2)
    t_w  = @view t[half:end]
    Cl_w = @view Cl[half:end]
    Cl_m = mean(Cl_w)
    crossings = Int[]
    for i in 2:length(Cl_w)
        if Cl_w[i-1] < Cl_m && Cl_w[i] ≥ Cl_m   # up-going
            # linear interpolation in time for the crossing
            r = (Cl_m - Cl_w[i-1]) / (Cl_w[i] - Cl_w[i-1])
            push!(crossings, i - 1)             # store sample-index
            crossings_t = nothing
        end
    end
    M = length(crossings)
    M < 2 && return (NaN, M)
    Δt = t_w[crossings[end]] - t_w[crossings[1]]
    f  = (M - 1) / Δt
    St = f * D / U
    return (St, M)
end

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

function main()
    isfile(OF_FORCES) || error("OpenFOAM forceCoeffs not found: $OF_FORCES")
    isfile(WL_CSV)    || error("WaterLily CSV not found: $WL_CSV")

    of = read_openfoam_forces(OF_FORCES)
    wl = read_waterlily_csv(WL_CSV)

    # Normalize OpenFOAM time to t U/D for direct comparison.
    of_t_nondim = of.t .* (OF_U / OF_D)

    @printf("\n=== Cylinder Re=100: WaterLily vs OpenFOAM ===\n\n")
    @printf("%-24s %14s %14s\n", "Metric", "WaterLily", "OpenFOAM")
    @printf("%-24s %14s %14s\n", "------", "---------", "--------")

    # Means
    Cd_wl = tail_mean(wl.Cd)
    Cd_of = tail_mean(of.Cd)
    Cl_wl = tail_mean(wl.Cl)
    Cl_of = tail_mean(of.Cl)
    @printf("%-24s %14.4f %14.4f\n", "Cd (mean, t>=t_mid)", Cd_wl, Cd_of)
    @printf("%-24s %14.4f %14.4f\n", "Cl (mean, t>=t_mid)", Cl_wl, Cl_of)

    # Cl peak-to-peak (proxy for shedding amplitude)
    Cl_pp_wl = let n = length(wl.Cl); half = max(1, n÷2)
        maximum(@view wl.Cl[half:end]) - minimum(@view wl.Cl[half:end])
    end
    Cl_pp_of = let n = length(of.Cl); half = max(1, n÷2)
        maximum(@view of.Cl[half:end]) - minimum(@view of.Cl[half:end])
    end
    @printf("%-24s %14.4f %14.4f\n", "Cl peak-to-peak",      Cl_pp_wl, Cl_pp_of)

    # Strouhal — WaterLily times are already in tU/D (sim_time), OpenFOAM
    # gets the dimensional D, U.
    St_wl, nx_wl = strouhal_from_lift(wl.t,         wl.Cl; D=1.0, U=1.0)
    St_of, nx_of = strouhal_from_lift(of_t_nondim, of.Cl; D=1.0, U=1.0)
    @printf("%-24s %14.4f %14.4f\n", "Strouhal (St)",        St_wl, St_of)
    @printf("%-24s %14d %14d\n",     "Cl crossings counted", nx_wl, nx_of)

    # Time-window of each tail
    twl_lo, twl_hi = let n = length(wl.t); half = max(1, n÷2)
        wl.t[half], wl.t[end]
    end
    tof_lo, tof_hi = let n = length(of_t_nondim); half = max(1, n÷2)
        of_t_nondim[half], of_t_nondim[end]
    end
    @printf("%-24s %14s %14s\n", "tU/D window",
            @sprintf("%5.1f..%5.1f", twl_lo, twl_hi),
            @sprintf("%5.1f..%5.1f", tof_lo, tof_hi))

    # Reference (Williamson 1996 / standard published values)
    println()
    @printf("Reference (Williamson 1996, Park-Kwon-Choi 1998):\n")
    @printf("  Cd ≈ 1.32–1.40\n")
    @printf("  St ≈ 0.165\n")
    @printf("  Cl peak-to-peak ≈ 0.60–0.66\n")
    println()

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
