"""
PGCP webapp server.

Loads all presets under `../presets/`, holds their converged BDDs in memory,
and exposes 4 HTTP endpoints (see README). Static files are served from
`../static/`.
"""

using HTTP
using Oxygen
using JSON
using JSON3

include("regions.jl")

const STATIC_DIR  = normpath(joinpath(@__DIR__, "..", "static"))
const PRESETS_DIR = normpath(joinpath(@__DIR__, "..", "presets"))

const REGIONS = load_all_regions(PRESETS_DIR)

const MIME_BY_EXT = Dict(
    ".html" => "text/html; charset=utf-8",
    ".js"   => "application/javascript; charset=utf-8",
    ".css"  => "text/css; charset=utf-8",
    ".json" => "application/json; charset=utf-8",
    ".svg"  => "image/svg+xml",
    ".png"  => "image/png",
    ".ico"  => "image/x-icon",
)

function mime_for(path::AbstractString)::String
    _, ext = splitext(path)
    return get(MIME_BY_EXT, lowercase(ext), "application/octet-stream")
end

function serve_static(req_path::AbstractString)
    # Default to index.html.
    rel = isempty(req_path) || req_path == "/" ? "index.html" : lstrip(req_path, '/')
    # Prevent path traversal.
    candidate = normpath(joinpath(STATIC_DIR, rel))
    if !startswith(candidate, STATIC_DIR) || !isfile(candidate)
        return HTTP.Response(404, "not found")
    end
    body = read(candidate)
    return HTTP.Response(200, ["Content-Type" => mime_for(candidate)], body=body)
end

# ============================== Routes ==============================

@get "/" function(req::HTTP.Request)
    return serve_static("index.html")
end

@get "/index.html" function(req::HTTP.Request)
    return serve_static("index.html")
end

@get "/static/*" function(req::HTTP.Request)
    # URI like /static/app.js  -> read static/app.js
    uri = HTTP.URI(req.target).path
    rel = replace(uri, r"^/static/?" => "")
    return serve_static(rel)
end

@get "/api/regions" function(req::HTTP.Request)
    list = [
        Dict(
            "id" => r.id,
            "display_name" => r.display_name,
            "subtitle" => r.subtitle,
            "n_sensors" => length(r.art.circles),
            "bbox_latlon" => r.bbox_latlon,
        )
        for r in values(REGIONS)
    ]
    sort!(list, by = x -> x["id"])
    return list
end

@get "/api/regions/{id}" function(req::HTTP.Request, id::String)
    if !haskey(REGIONS, id)
        return HTTP.Response(404, "unknown region: $id")
    end
    return region_payload(REGIONS[id])
end

@post "/api/regions/{id}/evaluate" function(req::HTTP.Request, id::String)
    if !haskey(REGIONS, id)
        return HTTP.Response(404, "unknown region: $id")
    end
    body = JSON3.read(String(req.body))
    pk = haskey(body, :pk) ? Float64(body.pk) : REGIONS[id].default_pk
    # Force the element type so an empty `disabled` array doesn't end up as Any[].
    disabled = String[String(s) for s in (haskey(body, :disabled) ? body.disabled : [])]
    t0 = time()
    R, importance = evaluate_region(REGIONS[id], pk, disabled)
    elapsed_ms = (time() - t0) * 1000
    return Dict(
        "R" => R,
        "pk" => pk,
        "disabled" => disabled,
        "importance" => importance,
        "elapsed_ms" => round(elapsed_ms; digits = 3),
    )
end

# Replan: caller submits a full sensor list (possibly with adds, moves, or
# deletes relative to the current state). We re-triangulate, re-solve, and
# swap the active BDD. Subsequent /evaluate hits the new BDD. Concurrency
# safety is handled inside `replan_region!` via the per-region ReentrantLock
# (CUDD/MiniCUDD is not thread-safe; the lock also covers evaluate).

@post "/api/regions/{id}/replan" function(req::HTTP.Request, id::String)
    if !haskey(REGIONS, id)
        return HTTP.Response(404, "unknown region: $id")
    end
    body = JSON3.read(String(req.body))
    if !haskey(body, :sensors)
        return HTTP.Response(400, "missing 'sensors'")
    end
    sensors = Dict{String,Any}[]
    for s in body.sensors
        push!(sensors, Dict{String,Any}(
            "name"   => String(s.name),
            "x"      => Float64(s.x),
            "y"      => Float64(s.y),
            "radius" => Float64(s.radius),
        ))
    end
    pk = haskey(body, :pk) ? Float64(body.pk) : REGIONS[id].default_pk
    disabled = String[String(s) for s in (haskey(body, :disabled) ? body.disabled : [])]

    if isempty(sensors)
        return HTTP.Response(400, "sensors[] is empty")
    end

    t0 = time()
    result = replan_region!(REGIONS[id], sensors, pk, disabled)
    elapsed_ms = (time() - t0) * 1000

    return Dict(
        "R" => result.R,
        "pk" => pk,
        "disabled" => disabled,
        "importance" => result.importance,
        "convergence" => result.convergence,
        "n_sensors" => result.n_sensors,
        "triangulation" => Dict(
            "vertices" => result.vertices,
            "triangles" => result.triangles,
        ),
        "sensors" => sensors,
        "min_path_sets" => result.min_path_sets,
        "min_cut_sets"  => result.min_cut_sets,
        "elapsed_ms" => round(elapsed_ms; digits = 3),
    )
end

# Health check.
@get "/api/health" function(::HTTP.Request)
    return Dict("ok" => true, "n_regions" => length(REGIONS))
end

@info "PGCP webapp starting" host="0.0.0.0" port=8080 n_regions=length(REGIONS)
serve(host = "0.0.0.0", port = 8080, async = false)
