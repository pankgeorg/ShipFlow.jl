module Harness

export run_foam_tutorial, read_volVectorField, read_volScalarField

const FOAM_IMAGE = "openfoam/openfoam11-paraview510"
const FOAM_ENV   = "/opt/openfoam11/etc/bashrc"

"""
    run_foam_tutorial(tutorial_path; case_dir, image=FOAM_IMAGE, env=FOAM_ENV,
                      pre_steps=["blockMesh"], main="foamRun")

Run an OpenFOAM tutorial inside a Docker container and return the
absolute path to the populated case directory.

- `tutorial_path` is interpreted relative to the locally cloned OpenFOAM
  tutorials tree at `OpenFOAM/tutorials/`. Example: `"incompressibleFluid/cavity"`.
- `case_dir` is where the case will live on the host (will be created
  and made world-writable so the container's user can write to it).
- `image` is the Docker image (default: OpenFOAM 11 Foundation).
- `pre_steps` defaults to `["blockMesh"]`. Pass `[]` to skip.
- `main` is the solver command (default `foamRun`, which dispatches to
  the `solver` keyword in `controlDict`).

The harness streams stdout/stderr to log files inside the case dir
(`log.<step>` for each pre-step and `log.foamRun` for the main run)
and throws on non-zero exit. Caller is responsible for selecting which
time step to inspect.
"""
function run_foam_tutorial(tutorial_path::AbstractString;
                           case_dir::AbstractString,
                           tutorials_root::AbstractString = joinpath(
                               dirname(dirname(@__DIR__)), "..", "OpenFOAM", "tutorials"
                           ),
                           image::AbstractString = FOAM_IMAGE,
                           env::AbstractString = FOAM_ENV,
                           pre_steps = ["blockMesh"],
                           main::AbstractString = "foamRun")
    src = joinpath(tutorials_root, tutorial_path)
    isdir(src) || error("Tutorial not found: $src")
    mkpath(case_dir)
    # Mirror the tutorial files into case_dir. Use cp -RL to dereference
    # any symlinks the tutorial may use; ignore failures from already-existing files.
    run(`cp -RL $src/. $case_dir`)
    # Make the case world-writable so the container's non-root user can write logs.
    run(`chmod -R 777 $case_dir`)

    # Build the script to run inside the container.
    pre_cmds = join((string(s, " >log.", s, " 2>&1") for s in pre_steps), " && ")
    main_cmd = "$main >log.foamRun 2>&1"
    inner = "source $env && cd /case && " *
            (isempty(pre_steps) ? "" : pre_cmds * " && ") * main_cmd

    cmd = `docker run --rm --entrypoint /bin/bash -v $case_dir:/case $image -c $inner`
    run(cmd)
    return case_dir
end

"""
    read_volVectorField(path) -> (n, vectors::Vector{NTuple{3,Float64}})

Parse an OpenFOAM ASCII `volVectorField` file (e.g. `<case>/<time>/U`)
into a flat list of cell-centered vectors. Boundary patches are
ignored — only `internalField` is returned.
"""
function read_volVectorField(path::AbstractString)
    raw = read(path, String)
    return _parse_internal_vector_list(raw)
end

"""
    read_volScalarField(path) -> (n, scalars::Vector{Float64})

Parse an OpenFOAM ASCII `volScalarField` file into a flat list of
cell-centered scalars.
"""
function read_volScalarField(path::AbstractString)
    raw = read(path, String)
    return _parse_internal_scalar_list(raw)
end

const _NUMBER_REGEX = r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?|[-+]?\.\d+(?:[eE][-+]?\d+)?"

function _parse_internal_vector_list(raw)
    # Locate internalField, then the nonuniform List<vector> N ( ... ) block.
    # Strategy: skip the header up to the opening paren that follows N,
    # then scan numbers, grouping them into triples (skipping the parens
    # that delimit each vector via the regex's tolerance for non-numeric
    # whitespace between tokens).
    m = match(r"internalField\s+nonuniform\s+List<[^>]+>\s*(\d+)\s*\("s, raw)
    if m === nothing
        m_u = match(r"internalField\s+uniform\s*\(([^)]+)\)", raw)
        m_u === nothing && error("internalField not found or unsupported")
        nums = parse.(Float64, split(m_u.captures[1]))
        return (1, [length(nums) ≥ 3 ? (nums[1], nums[2], nums[3]) : (nums[1], nums[2], 0.0)])
    end
    n = parse(Int, m.captures[1])
    body_start = m.offsets[1] + length(m.captures[1])     # past the count
    body_start = something(findnext('(', raw, body_start), body_start) + 1
    # Collect numbers up to 3n
    nums = Float64[]
    for tok in eachmatch(_NUMBER_REGEX, SubString(raw, body_start))
        push!(nums, parse(Float64, tok.match))
        length(nums) == 3n && break
    end
    length(nums) == 3n || @warn "parsed $(length(nums)) values, expected $(3n) for $n vectors"
    vectors = Vector{NTuple{3,Float64}}(undef, n)
    for k in 1:n
        vectors[k] = (nums[3k-2], nums[3k-1], nums[3k])
    end
    return (n, vectors)
end

function _parse_internal_scalar_list(raw)
    m = match(r"internalField\s+nonuniform\s+List<[^>]+>\s*(\d+)\s*\("s, raw)
    if m === nothing
        m_u = match(r"internalField\s+uniform\s+([^\s;]+)\s*;", raw)
        m_u === nothing && error("internalField not found or unsupported")
        return (1, [parse(Float64, m_u.captures[1])])
    end
    n = parse(Int, m.captures[1])
    body_start = m.offsets[1] + length(m.captures[1])
    body_start = something(findnext('(', raw, body_start), body_start) + 1
    nums = Float64[]
    for tok in eachmatch(_NUMBER_REGEX, SubString(raw, body_start))
        push!(nums, parse(Float64, tok.match))
        length(nums) == n && break
    end
    length(nums) == n || @warn "parsed $(length(nums)) scalars, expected $n"
    return (n, nums)
end

"""
    latest_time(case_dir) -> String

Return the numerically-largest time directory in `case_dir`. Used to
locate the final result of a steady or transient run.
"""
function latest_time(case_dir::AbstractString)
    times = String[]
    for entry in readdir(case_dir)
        isdir(joinpath(case_dir, entry)) || continue
        try
            parse(Float64, entry)
            push!(times, entry)
        catch
            # not a time dir
        end
    end
    isempty(times) && error("no time directories in $case_dir")
    return last(sort(times; by = s -> parse(Float64, s)))
end

end # module
