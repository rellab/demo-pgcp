#!/usr/bin/env julia
"""
Build a preset from an Overpass-API JSON file (one way / relation with
geometry). Replaces the Python pipeline of experiment5 (gen_sensors.py +
triangulate.py) with a Julia-only flow that reuses the in-process CDT and
solver already deployed in the webapp container.

Usage (inside the cudd-julia container):
  julia /app/webapp/scripts/build_preset.jl  \\
      /app/webapp/presets/_inputs/kasumi.raw.json  \\
      kasumi  \\
      "広島大学 霞キャンパス"  \\
      "Kasumi Campus, Hiroshima University"

Outputs five JSON files under webapp/presets/<preset_id>/ matching the schema
of higashi_hiroshima/.
"""

using JSON
using MiniCUDD
using Random
using GeometryBasics: Point

const HERE = @__DIR__
const SERVER = normpath(joinpath(HERE, "..", "server"))
include(joinpath(SERVER, "bdd_solver_triangle.jl"))
include(joinpath(SERVER, "bdd_analysis.jl"))
include(joinpath(SERVER, "triangulation.jl"))

using .BDDSolverTriangle
using .BDDSolverTriangle: Circle, Triangle3
using .BDDAnalysis

# ------- OSM polygon extraction -------

function load_osm_polygon(path::String)
    d = JSON.parsefile(path)
    els = d["elements"]
    isempty(els) && error("no elements in $path")
    el = els[1]
    geom = el["geometry"]
    return [(c["lat"], c["lon"]) for c in geom], get(el, "tags", Dict{String,Any}())
end

# ------- Projection (local equirectangular, isotropic-meter) -------

function project_to_normalized(polygon_latlon)
    lats = [p[1] for p in polygon_latlon]
    lons = [p[2] for p in polygon_latlon]
    min_lat, max_lat = minimum(lats), maximum(lats)
    min_lon, max_lon = minimum(lons), maximum(lons)
    ref_lat = (min_lat + max_lat) / 2
    MPD = 111319.0
    cos_ref = cos(ref_lat * pi / 180)
    width_m = (max_lon - min_lon) * MPD * cos_ref
    height_m = (max_lat - min_lat) * MPD
    side_meters = width_m
    h_ratio = height_m / width_m
    polygon_normalized = [
        [(lon - min_lon) * MPD * cos_ref / side_meters,
         (lat - min_lat) * MPD / side_meters]
        for (lat, lon) in polygon_latlon
    ]
    # Strip the closing duplicate (OSM ways usually close the ring).
    if length(polygon_normalized) > 1 && polygon_normalized[1] == polygon_normalized[end]
        polygon_normalized = polygon_normalized[1:end-1]
    end
    return polygon_normalized, side_meters, h_ratio,
           (min_lat = min_lat, max_lat = max_lat, min_lon = min_lon, max_lon = max_lon)
end

# ------- Geometry: point-in-polygon + PDS dart-throwing -------

function point_in_polygon(x::Float64, y::Float64, polygon::Vector{Vector{Float64}})
    n = length(polygon)
    inside = false
    j = n
    for i in 1:n
        xi, yi = polygon[i][1], polygon[i][2]
        xj, yj = polygon[j][1], polygon[j][2]
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

function pds_sample(polygon::Vector{Vector{Float64}}, h_ratio::Float64;
                     min_dist::Float64 = 0.089, max_failures::Int = 10000,
                     n_max::Int = 200, seed::Int = 42)
    rng = MersenneTwister(seed)
    samples = Tuple{Float64,Float64}[]
    failures = 0
    while failures < max_failures && length(samples) < n_max
        x = rand(rng)
        y = rand(rng) * h_ratio
        if !point_in_polygon(x, y, polygon)
            failures += 1
            continue
        end
        ok = true
        for s in samples
            if (s[1] - x)^2 + (s[2] - y)^2 < min_dist^2
                ok = false
                break
            end
        end
        if !ok
            failures += 1
            continue
        end
        push!(samples, (x, y))
        failures = 0
    end
    return samples
end

function assign_radii(n::Int; large_frac = 0.5, small = 0.089, large = 0.178, seed = 42)
    rng = MersenneTwister(seed + 1)
    n_large = min(round(Int, n * large_frac), n)
    radii = fill(small, n)
    perm = randperm(rng, n)
    for i in firstindex(perm):firstindex(perm) + n_large - 1
        radii[perm[i]] = large
    end
    return radii
end

# ------- Minimal-set summary (same shape as preset precomputed.json) -------

function set_summary(named_sets::Vector{Vector{String}}, total_count::Int)
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

# ------- Main -------

function main(input_path, preset_id, display_name, subtitle)
    polygon_latlon, _tags = load_osm_polygon(input_path)
    @info "Loaded $(length(polygon_latlon)) polygon vertices for $display_name"

    polygon_n, side_meters, h_ratio, bl =
        project_to_normalized(polygon_latlon)
    @info "side_meters = $(round(side_meters; digits=1)) m, h_ratio = $(round(h_ratio; digits=4))"

    samples = pds_sample(polygon_n, h_ratio; min_dist = 0.089)
    radii = assign_radii(length(samples))
    sensors = [Dict{String,Any}(
        "name" => string(i),
        "x" => samples[i][1],
        "y" => samples[i][2],
        "radius" => radii[i],
    ) for i in eachindex(samples)]
    @info "Generated $(length(sensors)) sensors via PDS"

    if isempty(sensors)
        error("PDS produced 0 sensors — polygon may be too small relative to min_dist")
    end

    vertices, triangles = cdt_polygon_with_steiners(
        polygon_n,
        [(x = s["x"], y = s["y"]) for s in sensors])
    @info "Mesh: $(length(vertices)) vertices, $(length(triangles)) triangles"

    circles = [Circle(s["name"], Point(s["x"], s["y"]), s["radius"]) for s in sensors]
    sort!(circles, by = c -> c.center[1])
    vars = Dict{String,Int}()
    for (i, c) in enumerate(circles)
        vars[c.name] = i - 1
    end

    init_triangles = Triangle3[]
    for t in triangles
        p1 = Point(vertices[t[1] + 1][1], vertices[t[1] + 1][2])
        p2 = Point(vertices[t[2] + 1][1], vertices[t[2] + 1][2])
        p3 = Point(vertices[t[3] + 1][1], vertices[t[3] + 1][2])
        push!(init_triangles, Triangle3((p1, p2, p3)))
    end

    M = MiniCUDD.Manager(nvars = max(1, length(circles)))
    res = nothing
    final_ml = 0
    for ml in 8:14
        @info "Solving at maxLevel=$ml ..."
        res = BDDSolverTriangle.triangle_bdd_solver(
            M, vars, circles, init_triangles;
            maxlevel = ml, verbose = false, pb = nothing)
        final_ml = ml
        if res.conv
            @info "  Converged at maxLevel = $ml"
            break
        end
    end
    if res === nothing || !res.conv
        @warn "Did not converge by maxLevel=14; using last attempt (bounds may not coincide)"
    end

    phi = res.varphi1
    pk = 0.9
    pb = Dict{Int,Float64}(v => pk for v in values(vars))
    R = BDDSolverTriangle.prob(M, phi, pb)
    @info "R($pk) = $R"

    imp = BDDAnalysis.importances(M, phi, pb)
    importance_list = []
    for (i, c) in enumerate(circles)
        idx = i - 1
        push!(importance_list, Dict{String,Any}(
            "name"          => c.name,
            "x"             => c.center[1],
            "y"             => c.center[2],
            "radius"        => c.radius,
            "birnbaum"      => get(imp.birnbaum, idx, 0.0),
            "criticality1"  => get(imp.crit1, idx, 0.0),
            "criticality0"  => get(imp.crit0, idx, 0.0),
            "structure"     => 0.0,
        ))
    end

    inv = Dict(v => k for (k, v) in vars)
    tonames(sets) = [sort([inv[i] for i in s]) for s in sets]
    mp_bdd = BDDAnalysis.minsol(M, phi)
    mc_bdd = BDDAnalysis.minsol(M, BDDAnalysis.dual(M, phi))
    n_path = BDDAnalysis.count_solutions(M, mp_bdd)
    n_cut  = BDDAnalysis.count_solutions(M, mc_bdd)
    path_named = tonames(BDDAnalysis.enumerate_solutions(M, mp_bdd; limit = 200_000))
    cut_named  = tonames(BDDAnalysis.enumerate_solutions(M, mc_bdd; limit = 200_000))
    @info "min-path sets: $n_path,  min-cut sets: $n_cut"

    essential = String[s[1] for s in cut_named if length(s) == 1]

    base = joinpath(HERE, "..", "presets", preset_id)
    mkpath(base)

    meta = Dict{String,Any}(
        "id"               => preset_id,
        "display_name"     => display_name,
        "subtitle"         => subtitle,
        "side_meters"      => side_meters,
        "bbox_normalized"  => [0.0, 0.0, 1.0, h_ratio],
        "bbox_latlon"      => Dict("min_lat" => bl.min_lat, "max_lat" => bl.max_lat,
                                   "min_lon" => bl.min_lon, "max_lon" => bl.max_lon),
        "n_sensors"        => length(sensors),
        "default_pk"       => pk,
        "maxlevel"         => final_ml,
    )

    poly_data = Dict{String,Any}(
        "polygon_normalized" => polygon_n,
        "polygon_latlon"     => [[lat, lon] for (lat, lon) in polygon_latlon],
    )

    sensors_data = Dict{String,Any}(
        "circles"     => sensors,
        "reliability" => pk,
        "maxlevel"    => final_ml,
    )

    tri_data = Dict{String,Any}(
        "vertices"  => vertices,
        "triangles" => triangles,
    )

    pre_data = Dict{String,Any}(
        "R"                  => R,
        "default_pk"         => pk,
        "n_sensors"          => length(sensors),
        "convergence"        => res === nothing ? false : res.conv,
        "importance"         => importance_list,
        "essential_sensors"  => essential,
        "min_cut_sets"       => set_summary(cut_named, n_cut),
        "min_path_sets"      => set_summary(path_named, n_path),
        "pk_sweep"           => [],
    )

    open(joinpath(base, "meta.json"), "w")          do io; JSON.print(io, meta,          2); end
    open(joinpath(base, "polygon.json"), "w")       do io; JSON.print(io, poly_data,     2); end
    open(joinpath(base, "sensors.json"), "w")       do io; JSON.print(io, sensors_data,  2); end
    open(joinpath(base, "triangulation.json"), "w") do io; JSON.print(io, tri_data,      2); end
    open(joinpath(base, "precomputed.json"), "w")   do io; JSON.print(io, pre_data,      2); end

    @info "wrote 5 files to $base"
    @info "  R = $R, n_sensors = $(length(sensors)), n_min_cut = $n_cut, n_min_path = $n_path"
    @info "  essential = $essential"
end

if length(ARGS) < 4
    error("Usage: julia build_preset.jl <input.json> <preset_id> <display_name> <subtitle>")
end
main(ARGS[1], ARGS[2], ARGS[3], ARGS[4])
