#!/usr/bin/env julia
#
# Plot J_self vs Cb from runs/cb_sweep/scans.csv.

import Pkg
Pkg.activate(; temp=true)
Pkg.add(["CairoMakie", "DelimitedFiles"]; io=devnull)

using CairoMakie, DelimitedFiles

const CSV = abspath(joinpath(@__DIR__, "..", "runs", "cb_sweep", "scans.csv"))
const OUT = abspath(joinpath(@__DIR__, "..", "runs", "cb_sweep", "Jself_vs_Cb.png"))
isfile(CSV) || error("$CSV not found — run containership_cb_sweep.jl first.")

raw, hdr = readdlm(CSV, ','; header=true)
hdr = vec(string.(hdr))
i_par = findfirst(==("par_frac"), hdr)
i_Cb  = findfirst(==("Cb"), hdr)
i_J   = findfirst(==("J"), hdr)
i_TmD = findfirst(==("T-D"), hdr)

# Group by Cb, then linear-interp the J at which T-D crosses zero.
Cbs = sort(unique(Float64.(raw[:, i_Cb])))
Jselfs = Float64[]
for cb in Cbs
    mask = Float64.(raw[:, i_Cb]) .≈ cb
    Js = Float64.(raw[mask, i_J])
    TmDs = Float64.(raw[mask, i_TmD])
    perm = sortperm(Js)
    Js = Js[perm]; TmDs = TmDs[perm]
    j_self = NaN
    for i in 1:length(Js) - 1
        if TmDs[i] * TmDs[i+1] < 0
            t = -TmDs[i] / (TmDs[i+1] - TmDs[i])
            j_self = Js[i] + t * (Js[i+1] - Js[i])
            break
        end
    end
    push!(Jselfs, j_self)
end

# Linear fit J_self ≈ a + b·Cb (least squares)
A = hcat(ones(length(Cbs)), Cbs)
β = A \ Jselfs
a, b = β[1], β[2]

fig = Figure(size=(680, 420))
ax = Axis(fig[1, 1];
    xlabel="Block coefficient Cb", ylabel="Self-prop J",
    title="J_self vs Cb (Containership family, Fr=0.25)")
scatter!(ax, Cbs, Jselfs; color=:steelblue, markersize=14)
Cb_grid = range(extrema(Cbs)..., length=80)
lines!(ax, Cb_grid, a .+ b .* Cb_grid;
    color=:tomato, linestyle=:dash,
    label="fit: $(round(a; digits=3)) + ($(round(b; digits=3)))·Cb")
axislegend(ax; position=:lb)
save(OUT, fig)
println("Wrote $OUT")
