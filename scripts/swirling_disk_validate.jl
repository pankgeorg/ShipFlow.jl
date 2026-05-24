#!/usr/bin/env julia
#
# Validate SwirlingDisk: run in a uniform free stream (no hull, no
# free surface), then plot u_θ vs r at a downstream cross-section.
# Expect monotonic non-zero tangential velocity behind the disk and
# zero swirl far upstream.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path=joinpath(@__DIR__, "..", "..", "Propellers.jl")),
]; io=devnull)
Pkg.add("CairoMakie"; io=devnull)
println("packages added"); flush(stdout)

using WaterLily, Propellers
using Propellers: SVector, StaticArrays
using Printf, CairoMakie

const NX = parse(Int, get(ENV, "WL_NX", "96"))
const NY = parse(Int, get(ENV, "WL_NY", "48"))
const NZ = parse(Int, get(ENV, "WL_NZ", "48"))
const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "150"))

const U∞ = 1f0
const Re = 5000f0
const L_c = 16f0
const ν_c = U∞ * L_c / Re

const disk_xc = Float32(NX/4)
const disk_yc = Float32(NY/2)
const disk_zc = Float32(NZ/2)
const R_d = 6f0
const W_d = 2f0
const C_T = 0.6f0
const thrust = 0.5f0 * π * R_d^2 * U∞^2 * C_T
# Pick torque such that Ω·R/U ≈ 0.5 (typical advance coefficient).
# Tangential force scales as ω·r ⇒ torque ≈ ∫ ω r · 2π r dr = 2π/3 ω R³.
# Choose torque = thrust * R/2 as a reasonable swirl magnitude.
const torque = thrust * R_d / 2

@printf "=== SwirlingDisk validation ===\n"
@printf "  Grid     = %d × %d × %d\n" NX NY NZ
@printf "  Disk     = R=%.1f, w=%.1f at (%.1f, %.1f, %.1f)\n" R_d W_d disk_xc disk_yc disk_zc
@printf "  Thrust   = %.3f (C_T=%.2f)\n" thrust C_T
@printf "  Torque   = %.3f\n" torque
@printf "  Steps    = %d\n\n" NSTEPS
flush(stdout)

disk = SwirlingDisk(
    center = SVector(disk_xc, disk_yc, disk_zc),
    axis   = SVector(1f0, 0f0, 0f0),
    R      = R_d, w = W_d,
    thrust = thrust, torque = torque,
)

sim = WaterLily.Simulation((NX, NY, NZ),
    (U∞, 0f0, 0f0), L_c;
    T = Float32, ν = ν_c, Δt = 0.5f0,
    exitBC = true, U = U∞,
)

disk_udf = (flow, t; kwargs...) -> disk(flow, t; kwargs...)

@info "Running $NSTEPS steps…"; flush(stdout)
for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; udf=disk_udf,
        pois_tol=1f-6, pois_itmx=50)
    if step % 25 == 0
        @printf "  step=%3d  |u|=%.3f\n" step maximum(abs, sim.flow.u)
        flush(stdout)
    end
end

# Sample u_θ on a y-z cross-section downstream of the disk.
function sample_swirl(u, x_slice::Int, yc, zc)
    NX, NY, NZ, _ = size(u)
    rs = Float32[]; uθs = Float32[]
    for j in 2:NY-1, k in 2:NZ-1
        y = Float32(j - 1.5) - yc
        z = Float32(k - 1.5) - zc
        r = sqrt(y*y + z*z)
        # Tangential direction (axis=x, so ê_θ = (0, -z/r, y/r))
        if r < 1e-3; continue; end
        # u at this cell: u_y and u_z faces — interpolate to cell center.
        uy = 0.5 * (u[x_slice, j, k, 2] + u[x_slice, j+1, k, 2])
        uz = 0.5 * (u[x_slice, j, k, 3] + u[x_slice, j, k+1, 3])
        uθ = (-z * uy + y * uz) / r
        push!(rs, r); push!(uθs, uθ)
    end
    return rs, uθs
end

x_down = round(Int, disk_xc + R_d)         # 1 R behind the disk
rs_dn, uθ_dn = sample_swirl(sim.flow.u, x_down, disk_yc, disk_zc)
x_up   = round(Int, disk_xc - R_d)         # 1 R ahead of the disk
rs_up, uθ_up = sample_swirl(sim.flow.u, x_up, disk_yc, disk_zc)

# Bin by radius for a clean profile
function bin_radial(rs, vs; nbins=20, rmax=10f0)
    edges = range(0, rmax; length=nbins+1)
    centers = Float32[]; means = Float32[]; ns = Int[]
    for i in 1:nbins
        mask = (rs .>= edges[i]) .& (rs .< edges[i+1])
        n = count(mask)
        if n > 0
            push!(centers, Float32(0.5*(edges[i]+edges[i+1])))
            push!(means, mean(vs[mask]))
            push!(ns, n)
        end
    end
    return centers, means, ns
end

using Statistics
r_dn_c, uθ_dn_c, _ = bin_radial(rs_dn, uθ_dn; rmax=Float32(NY/2))
r_up_c, uθ_up_c, _ = bin_radial(rs_up, uθ_up; rmax=Float32(NY/2))

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "swirling_disk"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "uθ_vs_r.png")

fig = Figure(size=(800, 500))
ax = Axis(fig[1,1];
    xlabel="r (cells)", ylabel="⟨u_θ⟩ (cells/Δt)",
    title=@sprintf("SwirlingDisk swirl profile — C_T=%.2f, Q/T=%.2f, Re=%.0f", C_T, torque/thrust, Re))
lines!(ax, r_dn_c, uθ_dn_c; label="1·R downstream", color=:tomato, linewidth=2.5)
lines!(ax, r_up_c, uθ_up_c; label="1·R upstream",   color=:steelblue, linewidth=2.5, linestyle=:dash)
vlines!(ax, [R_d]; color=:black, linestyle=:dot, label="disk R")
hlines!(ax, [0]; color=:grey)
axislegend(ax, position=:rb)

save(out, fig)
@printf "\nWrote %s\n" out

# Quick numerical summary
peak_dn, idx = findmax(abs, uθ_dn_c)
peak_up = maximum(abs, uθ_up_c)
@printf "Peak |u_θ| downstream = %.4f at r=%.1f\n" peak_dn r_dn_c[idx]
@printf "Peak |u_θ| upstream   = %.4f  (should be ~0)\n" peak_up
@printf "Sign of peak u_θ_dn   = %s  (positive = right-hand swirl about +x)\n" (uθ_dn_c[idx] > 0 ? "+ ✓" : "−")
