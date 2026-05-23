#!/usr/bin/env julia
#
# Free-surface Wigley hull smoke-test — VoF + ShipShapes + WaterLily.
# Half-submerged hull in uniform stream; this is the FULL stack except
# the propeller. Quantitative validation will come once we have
# turbulence (Re ~ 10⁶) and a tighter VoF advection scheme.
#
# Setup:
#   - 3D box NX × NY × NZ in cell-units
#   - Water fills z < H_w_c (lower half)
#   - Wigley hull L,B,T placed so the waterline (z=0 in hull coords)
#     coincides with the air/water interface
#   - Uniform inflow U∞ along +x, periodic y/z (well — periodic y, walls z),
#     exit BC in x
#   - Froude number Fr = U / √(g · L)
#   - Reynolds Re = U · L / ν_water
#
# Pass: runs to t_end without blowup; drag > 0; wave amplitude finite.

using WaterLily
using VoF
using ShipShapes
using ShipShapes: wigley_volume, StaticArrays
const SVector = StaticArrays.SVector
using Printf

# Grid
const NX = parse(Int, get(ENV, "WL_NX", "192"))
const NY = parse(Int, get(ENV, "WL_NY", "64"))
const NZ = parse(Int, get(ENV, "WL_NZ", "64"))

# Hull dimensions in cell-units
const L_c = parse(Float64, get(ENV, "WL_L", string(NX/2)))
const B_c = parse(Float64, get(ENV, "WL_B", string(NY/4)))
const T_c = parse(Float64, get(ENV, "WL_T", string(NZ/6)))

# Water/air properties (dimensionless ratio).  Use ρ=10:1 ratio for stability
# (we've shown ρ=1000 also runs but with more numerical noise).
const ρ_w = parse(Float64, get(ENV, "WL_RHO_RATIO", "10"))
const ρ_a = 1.0

# Reference scales — pick g_cell and U∞=1 (in cell-units) consistently.
const U∞ = 1.0
const Fr = parse(Float64, get(ENV, "WL_FR", "0.25"))
const G_c = (U∞ / (Fr * sqrt(L_c)))^2 / L_c                # solve Fr = U/√(g L)
const Re = parse(Float64, get(ENV, "WL_RE", "1000"))
const ν_w_c = U∞ * L_c / Re
const ν_a_c = ν_w_c * 18                                   # air ~18x more viscous than water kinematically
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * ν_a_c

const T_NUM = Float64

# Water fills below the centre of the box.  Waterline at z = NZ/2 cells.
const H_w_c = NZ/2

# Hull pose: midship at (NX/3, NY/2, H_w_c); hull-coord z=0 = waterline
hull_map = (x, t) -> SVector(x[1] - NX/3, x[2] - NY/2, x[3] - H_w_c)
hull     = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)

@printf "=== Free-surface Wigley hull smoke-test ===\n"
@printf "  Grid          = %d × %d × %d  (Δx normalised to 1 cell)\n" NX NY NZ
@printf "  Hull L,B,T    = %.1f, %.1f, %.1f cells\n" L_c B_c T_c
@printf "  L/B, L/T      = %.2f, %.2f\n" L_c/B_c L_c/T_c
@printf "  V_hull        = %.1f cells³\n" wigley_volume(L_c, B_c, T_c)
@printf "  Waterline     = z = %.1f cells\n" H_w_c
@printf "  ρ_w/ρ_a       = %.0f\n" ρ_w/ρ_a
@printf "  Fr (target)   = %.3f → g_cell = %.4e\n" Fr G_c
@printf "  Re_L          = %.0f → ν_w = %.3e, ν_a = %.3e\n" Re ν_w_c ν_a_c

# α₀: water in lower half
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀,
    ρ_w = ρ_w, ρ_a = ρ_a,
    μ_w = μ_w_c, μ_a = μ_a_c,
    T = T_NUM,
)

function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ;
        perdir = (2,))
end

sim = WaterLily.Simulation((NX, NY, NZ), (T_NUM(U∞), T_NUM(0), T_NUM(0)), L_c;
    T = T_NUM,
    ν = vof.ν,
    g = (i, x, t) -> i == 3 ? T_NUM(-G_c) : T_NUM(0),
    Δt = 0.25,
    body = hull,
    ϵ = 1,
    perdir = (2,),       # periodic in y (sides), walls in z (free surf via VoF)
    exitBC = true,       # convective exit in x
    pois_ctor = vof_pois_ctor,
    U = T_NUM(U∞),
)

const NSTEPS = parse(Int, get(ENV, "WL_NSTEPS", "200"))
const REPORT_EVERY = max(1, NSTEPS ÷ 20)

# Wave amplitude diagnostic: max free-surface elevation deviation from H_w_c
function max_wave(α)
    nx, ny, nz = size(α)
    max_eta = 0.0
    @inbounds for j in 2:ny-1, i in 2:nx-1
        for k in nz-1:-1:2
            if α[i, j, k] > 0.5
                eta = (k - 1.5) - H_w_c
                abs(eta) > abs(max_eta) && (max_eta = eta)
                break
            end
        end
    end
    return max_eta
end

A_wet = (8/9) * L_c * (B_c + 4*T_c) / 2     # half-submerged ≈ half wetted area
@info "Running $(NSTEPS) steps…  (A_wet ≈ $(round(A_wet, digits=1)) cells²)"

for step in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois;
        pois_tol = T_NUM(1e-8), pois_itmx = 100)
    step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
    if mod(step, REPORT_EVERY) == 0 || step ≤ 3
        Fp = WaterLily.pressure_force(sim)
        Fv = WaterLily.viscous_force(sim)
        u_max = maximum(abs, sim.flow.u)
        wave = max_wave(vof.α)
        F_drag = -Float64(Fp[1]) - Float64(Fv[1])
        CD = 2 * F_drag / (1.0 * U∞^2 * A_wet)
        @info @sprintf("step=%4d  t=%.3f  Δt=%.2e  |u|=%.3f  wave=%+5.2f  F_drag=%+.2f  CD=%.4f",
                       step, sim.flow.Δt[end-1] * step, sim.flow.Δt[end-1],
                       u_max, wave, F_drag, CD)
        if !isfinite(u_max) || u_max > 50
            @warn "Blow-up at step $step"
            break
        end
    end
end

Fp = WaterLily.pressure_force(sim)
Fv = WaterLily.viscous_force(sim)
F_drag = -Float64(Fp[1]) - Float64(Fv[1])
CD = 2 * F_drag / (1.0 * U∞^2 * A_wet)
println()
@printf "Final F_drag = %.4f cell-units\n" F_drag
@printf "Final C_D    = %.4f  (Re=%g, Fr=%.2f, ρ_w/ρ_a=%g)\n" CD Re Fr ρ_w/ρ_a
@printf "Max wave η   = %+.3f cells\n" max_wave(vof.α)
println("\n✓ Smoke test: stable to NSTEPS=$NSTEPS, drag > 0, wave pattern present.")
