using Test
using ShipFlow
using ShipFlow.Harness

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
