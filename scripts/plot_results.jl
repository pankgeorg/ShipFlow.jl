#!/usr/bin/env julia
#
# Plot the two headline curves: self-propulsion (T, D vs C_T) and
# wave-resistance (drag vs Fr).  Writes PNGs to runs/.

import Pkg
Pkg.activate(; temp=true)
Pkg.add(["CairoMakie", "DelimitedFiles"]; io=devnull)

using CairoMakie
using DelimitedFiles

const ROOT = abspath(joinpath(@__DIR__, ".."))

# ------------------ self-propulsion ------------------
sp_csv = joinpath(ROOT, "runs", "wigley_selfprop_scan", "scan.csv")
sp = readdlm(sp_csv, ','; header=true)[1]
CTs, Ts, Ds, Diffs = sp[:,1], sp[:,2], sp[:,3], sp[:,4]

fig1 = Figure(size=(800, 500))
ax = Axis(fig1[1, 1],
    xlabel="C_T (thrust coefficient)",
    ylabel="force (cell-units)",
    title="Self-propulsion scan — Wigley + actuator disk + VoF + WALE")
lines!(ax, CTs, Ts, label="thrust (T)", color=:steelblue, linewidth=2.5)
scatter!(ax, CTs, Ts, color=:steelblue, markersize=8)
lines!(ax, CTs, Ds, label="hull drag (D)", color=:tomato, linewidth=2.5)
scatter!(ax, CTs, Ds, color=:tomato, markersize=8)
# Mark self-propulsion C_T ≈ 2.27 (intersection)
sp_CT = 2.27
sp_force = 70.0
scatter!(ax, [sp_CT], [sp_force], color=:black, markersize=14, marker=:star5,
    label="self-propulsion (C_T≈2.27)")
axislegend(ax; position=:lt)
fig1_path = joinpath(ROOT, "runs", "wigley_selfprop_scan", "self_propulsion.png")
save(fig1_path, fig1)
println("Wrote ", fig1_path)

# ------------------ Fr scan ------------------
fr_csv = joinpath(ROOT, "runs", "wigley_resistance_Fr", "drag_vs_Fr.csv")
fr = readdlm(fr_csv, ','; header=true)[1]
Frs, dr, CTr = fr[:,1], fr[:,2], fr[:,5]

fig2 = Figure(size=(800, 500))
ax2 = Axis(fig2[1, 1],
    xlabel="Froude number Fr = U/√(gL)",
    ylabel="C_T (= 2D / (ρ U² A_wet))",
    title="Wigley bare-hull resistance vs Fr — wave-making hump")
lines!(ax2, Frs, CTr, color=:darkgreen, linewidth=2.5)
scatter!(ax2, Frs, CTr, color=:darkgreen, markersize=10)
# Annotate the hump
peak_i = argmax(CTr)
peak_Fr = Frs[peak_i]
peak_C = CTr[peak_i]
scatter!(ax2, [peak_Fr], [peak_C], color=:black, markersize=14, marker=:star5,
    label="hump (Fr=$peak_Fr)")
axislegend(ax2; position=:lt)
fig2_path = joinpath(ROOT, "runs", "wigley_resistance_Fr", "drag_vs_Fr.png")
save(fig2_path, fig2)
println("Wrote ", fig2_path)
