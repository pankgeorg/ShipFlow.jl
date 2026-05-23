using Test
using ShipFlow
using ShipFlow.Harness

@testset "Integration: five-package stack on a tiny grid" begin
    # Smoke-test that the full stack — WaterLily + VoF + ShipShapes +
    # Turbulence + Propellers — boots and steps once. Uses a 24×16×16
    # grid so it runs in seconds. Asserts only finiteness; quantitative
    # validation lives in the RESULTS docs.
    using WaterLily
    using VoF
    using ShipShapes
    using ShipShapes: StaticArrays
    using Turbulence
    using Propellers
    SVector = Propellers.StaticArrays.SVector

    NX, NY, NZ = 24, 16, 16
    L_c, B_c, T_c = 12, 4, 3
    H_w_c = NZ/2
    xc, yc, zc = NX/3, NY/2, H_w_c
    Fr, Re = 0.25, 100
    U∞ = 1.0
    G_c = U∞^2 / (Fr^2 * L_c)
    ν_w_c = U∞ * L_c / Re
    ρ_w, ρ_a = 10.0, 1.0

    α₀(_i, x) = (x[3] ≤ H_w_c) ? 1.0 : 0.0
    map = (x, t) -> SVector(x[1] - xc, x[2] - yc, x[3] - zc)
    hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = map)
    vof = VoFFlow((NX, NY, NZ); α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a,
        μ_w = ρ_w * ν_w_c, μ_a = ρ_a * 18 * ν_w_c, T = Float64)
    turb = WALE((NX, NY, NZ); Cw = 0.5, ν₀ = 0.0)
    disk = ActuatorDisk(
        center = SVector(xc + L_c/2 + T_c/2, yc, H_w_c - T_c/2),
        axis   = SVector(1.0, 0.0, 0.0),
        R = T_c/2, w = 1.5, thrust = 0.1,
    )
    pois_ctor = flow -> begin
        L = similar(flow.μ₀)
        @inbounds for I in CartesianIndices(L)
            L[I] = flow.μ₀[I] * vof.L[I]
        end
        WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
    end
    sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0.0, 0.0), Float64(L_c);
        T = Float64, ν = turb.ν,
        g = (i, x, t) -> i == 3 ? -G_c : 0.0,
        Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
        pois_ctor = pois_ctor, U = U∞)

    disk_udf = (flow, t; kw...) -> disk(flow, t; kw...)
    for step in 1:3
        WaterLily.mom_step!(sim.flow, sim.pois;
            udf = disk_udf, pois_tol = 1e-6, pois_itmx = 50)
        step_vof!(vof, sim; dt = sim.flow.Δt[end-1])
        update_νt!(turb, sim.flow.u, vof.ν)
    end

    # Finiteness of every field used in the stack
    @test all(isfinite, vof.α)
    @test all(isfinite, vof.ν)
    @test all(isfinite, vof.L)
    @test all(isfinite, turb.ν)
    @test all(isfinite, sim.flow.u)
    @test all(isfinite, sim.flow.p)
    # Mass conservation (loose) — α should sum to roughly initial value
    α_sum = sum(vof.α[2:NX+1, 2:NY+1, 2:NZ+1])
    α0_sum = NX * NY * Int(H_w_c)   # initial water cells
    @test 0.9 * α0_sum ≤ α_sum ≤ 1.1 * α0_sum
    # Forces are finite
    @test all(isfinite, WaterLily.pressure_force(sim))
    @test all(isfinite, WaterLily.viscous_force(sim))
end

@testset "Harness" begin

    # The cavity run from iteration 6 left its output here; reuse it
    # for parser tests so we don't need to rerun Docker on every CI.
    cavity_dir = abspath(joinpath(@__DIR__, "..", "runs", "cavity"))

    @testset "latest_time" begin
        if isdir(cavity_dir)
            t = Harness.latest_time(cavity_dir)
            @test t == "10"   # endTime in the cavity controlDict
        else
            @info "Skipping latest_time test (cavity run not present)"
        end
    end

    @testset "read_volVectorField on cavity U" begin
        u_path = joinpath(cavity_dir, "10", "U")
        if isfile(u_path)
            n, vecs = read_volVectorField(u_path)
            @test n == 400               # 20×20 cavity mesh
            @test length(vecs) == 400
            # Lid-driven cavity at Re=10 (from tutorial) — max |U| should
            # be < 2 (the lid moves at U=1).
            speeds = [sqrt(v[1]^2 + v[2]^2 + v[3]^2) for v in vecs]
            @test maximum(speeds) < 2.0
            @test maximum(speeds) > 0.0  # not all zero
        else
            @info "Skipping U parser test (no $u_path)"
        end
    end

    @testset "read_volScalarField on cavity p" begin
        p_path = joinpath(cavity_dir, "10", "p")
        if isfile(p_path)
            n, vals = read_volScalarField(p_path)
            @test n == 400
            @test length(vals) == 400
            # Pressure field has finite range, not NaN
            @test all(isfinite, vals)
        else
            @info "Skipping p parser test (no $p_path)"
        end
    end

end
