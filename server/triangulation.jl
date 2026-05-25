"""
Constrained Delaunay triangulation of a polygon with sensor centers as
Steiner points, using DelaunayTriangulation.jl. Mirrors what
`pgcp-experiments/experiment5/scripts/triangulate.py` does with Shewchuk's
`triangle` library, but stays in-process and avoids the Python toolchain.
"""

using DelaunayTriangulation

"""
    cdt_polygon_with_steiners(polygon, sensors)

Triangulate `polygon` (a Vector of [x, y]) with `sensors` (a Vector of
NamedTuple-like with `.x` and `.y`) as Steiner points. Returns
`(vertices, triangles)` where:

  - `vertices :: Vector{Vector{Float64}}` — all vertex coordinates ([x, y]),
    starting with polygon vertices in input order, followed by sensor centers.
  - `triangles :: Vector{Vector{Int}}` — 0-based vertex indices of each
    interior triangle (`each_solid_triangle`), in CCW order. 0-based to match
    the existing mesh JSON schema in `presets/*/triangulation.json`.
"""
function cdt_polygon_with_steiners(polygon::Vector{Vector{Float64}},
                                    sensors::AbstractVector)
    # GeoJSON-style polygons often repeat the first vertex at the end to
    # close the ring. DelaunayTriangulation.jl rejects duplicate points, so
    # strip the trailing repeat. The boundary_nodes loop closes the ring.
    poly = collect(polygon)
    if length(poly) > 1 && poly[1] == poly[end]
        poly = poly[1:end-1]
    end
    n_poly = length(poly)

    points = Vector{NTuple{2,Float64}}(undef, n_poly + length(sensors))
    for i in 1:n_poly
        points[i] = (Float64(poly[i][1]), Float64(poly[i][2]))
    end
    for (j, s) in enumerate(sensors)
        points[n_poly + j] = (Float64(s.x), Float64(s.y))
    end

    # Closed polygon ring as the single outer boundary.
    boundary_nodes = [vcat(collect(1:n_poly), [1])]

    tri = DelaunayTriangulation.triangulate(points;
        boundary_nodes = boundary_nodes)

    triangles = Vector{Vector{Int}}()
    for T in DelaunayTriangulation.each_solid_triangle(tri)
        i, j, k = DelaunayTriangulation.triangle_vertices(T)
        # Convert to 0-based indices, matching the mesh JSON schema.
        push!(triangles, [i - 1, j - 1, k - 1])
    end

    vertices = [Float64[p[1], p[2]] for p in points]
    return (vertices, triangles)
end
