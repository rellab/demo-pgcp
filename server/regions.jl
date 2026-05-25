"""
RegionState: per-preset in-memory state held by the webapp.

The state is split into a stable "preset header" (polygon, bbox, display name,
…) and a swappable `SolverArtifact` that holds the current Manager / BDD /
mesh / circle list. Initial load builds an artifact from the on-disk preset;
the `/replan` endpoint builds a new artifact from a user-supplied sensor list
and swaps it in place. `/evaluate` always queries `region.art`.

Cleanup note: each replan creates a fresh MiniCUDD.Manager and discards the
old one. We rely on Julia GC to release the C-side allocations; for an
event-day demo this is fine.
"""

using JSON
using MiniCUDD
using GeometryBasics: Point

# Solver, BDD-analysis layer, and CDT wrapper are vendored into this
# directory so the webapp builds without the pgcp-experiments tree alongside.
# `bdd_solver_triangle.jl` and `bdd_analysis.jl` are byte-identical copies of
# the canonical files in `../pgcp-experiments/scripts/`; sync them by hand
# (or via the script in webapp/scripts/) whenever the upstream changes.
include(joinpath(@__DIR__, "bdd_solver_triangle.jl"))
include(joinpath(@__DIR__, "bdd_analysis.jl"))
include(joinpath(@__DIR__, "triangulation.jl"))

using .BDDSolverTriangle
using .BDDSolverTriangle: Circle, Triangle3
using .BDDAnalysis

mutable struct SolverArtifact
    M::MiniCUDD.Manager
    vars::Dict{String,Int}
    circles::Vector{Circle}
    phi::MiniCUDD.BDDNode
    convergence::Bool
    vertices::Vector{Vector{Float64}}
    triangles::Vector{Vector{Int}}   # 0-based triples
end

mutable struct RegionState
    id::String
    display_name::String
    subtitle::String
    default_pk::Float64
    maxlevel::Int

    polygon_normalized::Vector{Vector{Float64}}
    polygon_latlon::Vector{Vector{Float64}}
    bbox_normalized::Vector{Float64}
    bbox_latlon::Dict{String,Float64}
    side_meters::Float64

    # Sensor specs are plain dicts so we can round-trip them with the
    # frontend (and reset by re-applying the original list).
    original_sensors::Vector{Dict{String,Any}}
    current_sensors::Vector{Dict{String,Any}}

    art::SolverArtifact

    # Cached analysis of the original artifact (for the GET payload's
    # `precomputed` block). The frontend reuses this for the "initial state"
    # baseline; after a replan the frontend uses the /replan response.
    precomputed::Dict{String,Any}

    # CUDD / MiniCUDD is not thread-safe: concurrent operations on the same
    # Manager (or even read-allocating ops like `birnbaum` / `minsol`) can
    # corrupt internal unique tables. We serialize every BDD-touching
    # operation on this region behind a per-region lock. Reentrant so a
    # replan that calls evaluate_region inside can acquire it twice.
    lock::ReentrantLock
end

"""
    build_artifact(polygon, sensors, maxlevel) -> SolverArtifact

Build a fresh (Manager, BDD, mesh) tuple for the given sensor list. Performs:

  1. constrained Delaunay triangulation of `polygon` with sensor centers
     as Steiner points,
  2. variable allocation in CUDD ordered by sensor x-coordinate,
  3. a full run of the adaptive symbolic refinement solver to `maxlevel`.

`sensors` is a vector of dicts/NamedTuples with keys `"name"`, `"x"`, `"y"`,
`"radius"`.
"""
function build_artifact(polygon::Vector{Vector{Float64}},
                         sensors::AbstractVector,
                         maxlevel::Int)::SolverArtifact
    # Build Circle list, sorted by x.
    circles = [Circle(String(s["name"]), Point(Float64(s["x"]), Float64(s["y"])),
                       Float64(s["radius"])) for s in sensors]
    sort!(circles, by = c -> c.center[1])

    vars = Dict{String,Int}()
    for (i, c) in enumerate(circles)
        vars[c.name] = i - 1
    end

    # CDT (in-process, via DelaunayTriangulation.jl).
    sensor_pts = [(x = Float64(s["x"]), y = Float64(s["y"])) for s in sensors]
    vertices, triangles = cdt_polygon_with_steiners(polygon, sensor_pts)

    # Triangle3 list (1-based indices into `vertices`).
    init_triangles = Triangle3[]
    for t in triangles
        p1 = Point(vertices[t[1] + 1][1], vertices[t[1] + 1][2])
        p2 = Point(vertices[t[2] + 1][1], vertices[t[2] + 1][2])
        p3 = Point(vertices[t[3] + 1][1], vertices[t[3] + 1][2])
        push!(init_triangles, Triangle3((p1, p2, p3)))
    end

    M = MiniCUDD.Manager(nvars = max(1, length(circles)))
    res = BDDSolverTriangle.triangle_bdd_solver(
        M, vars, circles, init_triangles;
        maxlevel = maxlevel, verbose = false, pb = nothing)

    return SolverArtifact(M, vars, circles, res.varphi1, res.conv,
                          vertices, triangles)
end

"""
    load_region(preset_dir) -> RegionState

Read the five JSON files in `preset_dir` and build the in-memory state,
including a fresh BDD solve. Uses the on-disk triangulation when present
(matches the published case study); otherwise re-triangulates.
"""
function load_region(preset_dir::AbstractString)::RegionState
    meta = JSON.parsefile(joinpath(preset_dir, "meta.json"))
    polygon = JSON.parsefile(joinpath(preset_dir, "polygon.json"))
    sensors_in = JSON.parsefile(joinpath(preset_dir, "sensors.json"))
    pre = JSON.parsefile(joinpath(preset_dir, "precomputed.json"))

    raw_circles = Vector{Dict{String,Any}}(sensors_in["circles"])
    maxlevel = Int(get(meta, "maxlevel", 11))
    polygon_n = Vector{Vector{Float64}}(polygon["polygon_normalized"])

    art = build_artifact(polygon_n, raw_circles, maxlevel)

    return RegionState(
        meta["id"], meta["display_name"], get(meta, "subtitle", ""),
        Float64(get(meta, "default_pk", 0.9)), maxlevel,
        polygon_n,
        Vector{Vector{Float64}}(polygon["polygon_latlon"]),
        Float64.(meta["bbox_normalized"]),
        Dict(string(k) => Float64(v) for (k, v) in meta["bbox_latlon"]),
        Float64(get(meta, "side_meters", 1.0)),
        raw_circles, deepcopy(raw_circles),
        art, pre,
        ReentrantLock(),
    )
end

function load_all_regions(presets_dir::AbstractString)::Dict{String,RegionState}
    out = Dict{String,RegionState}()
    for entry in sort!(readdir(presets_dir))
        sub = joinpath(presets_dir, entry)
        if isdir(sub) && isfile(joinpath(sub, "meta.json"))
            @info "Loading preset" entry
            t0 = time()
            state = load_region(sub)
            t1 = time()
            R = BDDSolverTriangle.prob(state.art.M, state.art.phi, default_pb(state))
            @info "  loaded" preset=entry n_sensors=length(state.art.circles) R=R conv=state.art.convergence elapsed_sec=round(t1 - t0; digits=3)
            out[entry] = state
        end
    end
    return out
end

default_pb(region::RegionState) = default_pb(region.art, region.default_pk)
function default_pb(art::SolverArtifact, pk::Float64)
    pb = Dict{Int,Float64}()
    for (_, v) in art.vars
        pb[v] = pk
    end
    return pb
end

"""
    evaluate_region(region, pk, disabled) -> (R, importance)

`pk` is the common reliability of every (active) sensor. `disabled` is a list
of sensor names treated as failed (p = 0). Returns R and the Birnbaum
importance vector ∂R/∂p_k of every sensor, computed at the same pb.
"""
function evaluate_region(region::RegionState, pk::Float64,
                          disabled::AbstractVector{<:AbstractString})
    return lock(region.lock) do
        art = region.art
        pb = Dict{Int,Float64}()
        for (_, v) in art.vars
            pb[v] = pk
        end
        for name in disabled
            if haskey(art.vars, String(name))
                pb[art.vars[String(name)]] = 0.0
            end
        end
        R = BDDSolverTriangle.prob(art.M, art.phi, pb)
        imp_by_var = BDDAnalysis.birnbaum(art.M, art.phi, pb)
        importance = Dict{String,Float64}()
        for (name, vidx) in art.vars
            importance[name] = get(imp_by_var, vidx, 0.0)
        end
        return (R, importance)
    end
end

"""
    set_summary(named_sets, total_count) -> Dict

Build the JSON-shaped summary for one family of minimal sets (path or cut).
Matches the schema used in `presets/*/precomputed.json` (count / smallest /
mean_cardinality / …) so the frontend treats the live-recomputed sets the
same way as the on-disk precomputed ones.
"""
function set_summary(named_sets::Vector{Vector{String}}, total_count::Int)::Dict{String,Any}
    cards = sort(length.(named_sets))
    smallest = sort(named_sets, by = s -> (length(s), s))
    n_show = min(10, length(smallest))
    return Dict{String,Any}(
        "count"            => total_count,
        "count_enumerated" => length(named_sets),
        "truncated"        => length(named_sets) < total_count,
        "min_cardinality"  => isempty(cards) ? 0 : cards[1],
        "max_cardinality"  => isempty(cards) ? 0 : cards[end],
        "mean_cardinality" => isempty(cards) ? 0.0 : sum(cards) / length(cards),
        "smallest"         => smallest[1:n_show],
    )
end

"""
    compute_minsets_summary(art) -> (min_path_sets, min_cut_sets)

Cheap analysis of the current BDD: minimal solutions of Φ (= min-path sets)
and of its dual (= min-cut sets), returned in the same shape as the on-disk
precomputed JSON. `ENUM_LIMIT` bounds the listing so an extreme configuration
with millions of sets cannot stall the request — the `count` field is still
exact because it comes from a memoized DAG count.
"""
const ENUM_LIMIT = 200_000

function compute_minsets_summary(art::SolverArtifact)
    M, phi, vars = art.M, art.phi, art.vars
    inv = Dict(v => k for (k, v) in vars)
    tonames(sets) = [sort([inv[i] for i in s]) for s in sets]
    mp_bdd = BDDAnalysis.minsol(M, phi)
    mc_bdd = BDDAnalysis.minsol(M, BDDAnalysis.dual(M, phi))
    n_path = BDDAnalysis.count_solutions(M, mp_bdd)
    n_cut  = BDDAnalysis.count_solutions(M, mc_bdd)
    path_named = tonames(BDDAnalysis.enumerate_solutions(M, mp_bdd; limit = ENUM_LIMIT))
    cut_named  = tonames(BDDAnalysis.enumerate_solutions(M, mc_bdd; limit = ENUM_LIMIT))
    return (
        min_path_sets = set_summary(path_named, n_path),
        min_cut_sets  = set_summary(cut_named, n_cut),
    )
end

"""
    replan_region!(region, new_sensors, pk, disabled) -> NamedTuple

Replace the region's solver artifact with one built from `new_sensors`, then
evaluate R / Birnbaum / min-path / min-cut at (pk, disabled). The previous
Manager is dropped (reachable only as a return value, which is also dropped
at the end of the caller) so Julia GC can release CUDD memory.
"""
function replan_region!(region::RegionState,
                         new_sensors::Vector{Dict{String,Any}},
                         pk::Float64,
                         disabled::AbstractVector{<:AbstractString})
    return lock(region.lock) do
        new_art = build_artifact(region.polygon_normalized, new_sensors, region.maxlevel)
        region.art = new_art
        region.current_sensors = deepcopy(new_sensors)
        # evaluate_region also takes the lock, hence the ReentrantLock.
        R, importance = evaluate_region(region, pk, disabled)
        msets = compute_minsets_summary(new_art)
        return (
            R = R,
            importance = importance,
            convergence = new_art.convergence,
            n_sensors = length(new_art.circles),
            vertices = new_art.vertices,
            triangles = new_art.triangles,
            min_path_sets = msets.min_path_sets,
            min_cut_sets  = msets.min_cut_sets,
        )
    end
end

"""
Build the full GET /api/regions/{id} payload from `region`. The payload
reflects whatever the current artifact is, so after a replan the new mesh
and sensor list are returned.
"""
function region_payload(region::RegionState)::Dict{String,Any}
    return Dict{String,Any}(
        "id" => region.id,
        "display_name" => region.display_name,
        "subtitle" => region.subtitle,
        "n_sensors" => length(region.art.circles),
        "default_pk" => region.default_pk,
        "side_meters" => region.side_meters,
        "bbox_normalized" => region.bbox_normalized,
        "bbox_latlon" => region.bbox_latlon,
        "polygon_normalized" => region.polygon_normalized,
        "polygon_latlon" => region.polygon_latlon,
        "original_sensors" => region.original_sensors,
        "sensors" => region.current_sensors,
        "triangulation" => Dict(
            "vertices" => region.art.vertices,
            "triangles" => region.art.triangles,
        ),
        "precomputed" => region.precomputed,
        "solver" => Dict(
            "convergence" => region.art.convergence,
            "phi_nodes" => Int(MiniCUDD.dag_size(region.art.phi)),
        ),
    )
end
