module BDDSolverTriangle

using MiniCUDD
using Printf
using JSON
using ProgressMeter
using GeometryBasics: Point

# ===================== Types =====================
@inline getx(p::Point) = @inbounds p[1]
@inline gety(p::Point) = @inbounds p[2]

struct Circle
    name::String
    center::Point
    radius::Float64
end

struct Triangle3
    vertices::NTuple{3,Point}
end

@inline v1(t::Triangle3) = @inbounds t.vertices[1]
@inline v2(t::Triangle3) = @inbounds t.vertices[2]
@inline v3(t::Triangle3) = @inbounds t.vertices[3]

isinside(c::Circle, p::Point) =
    (getx(p) - getx(c.center))^2 + (gety(p) - gety(c.center))^2 <= c.radius^2

@inline midpoint(a::Point, b::Point) =
    Point((getx(a) + getx(b)) / 2, (gety(a) + gety(b)) / 2)

# 1-to-4 edge-midpoint subdivision.
#
#        v1                       v1
#        /\                       /\
#       /  \           →        m31--m12
#      /    \                   / \  / \
#     /      \                 /   \/   \
#    v3-------v2              v3---m23---v2
#
# Children: (v1, m12, m31), (v2, m23, m12), (v3, m31, m23), (m12, m23, m31).
function subdivide(t::Triangle3)
    a, b, c = v1(t), v2(t), v3(t)
    mab = midpoint(a, b)
    mbc = midpoint(b, c)
    mca = midpoint(c, a)
    return (
        Triangle3((a,   mab, mca)),
        Triangle3((b,   mbc, mab)),
        Triangle3((c,   mca, mbc)),
        Triangle3((mab, mbc, mca)),
    )
end

# ===================== BDD utilities =====================
function bdd_and(M::MiniCUDD.Manager, vs::Vector{MiniCUDD.BDDNode})
    isempty(vs) && return const1(M)
    acc = vs[1]
    for i in 2:length(vs)
        acc = MiniCUDD.bdd_and(M, acc, vs[i])
    end
    return acc
end

function bdd_or(M::MiniCUDD.Manager, vs::Vector{MiniCUDD.BDDNode})
    isempty(vs) && return const0(M)
    acc = vs[1]
    for i in 2:length(vs)
        acc = MiniCUDD.bdd_or(M, acc, vs[i])
    end
    return acc
end

# OR of a sorted list of variable indices; const0 if empty.
function or_of_vars(M::MiniCUDD.Manager, ids::Vector{Int})
    isempty(ids) && return MiniCUDD.const0(M)
    acc = var(M, ids[1])
    for i in 2:length(ids)
        acc = MiniCUDD.bdd_or(M, acc, var(M, ids[i]))
    end
    return acc
end

# Sorted-vector intersection (mutates inputs).
function set_intersection_vec(a::Vector{T}, b::Vector{T}) where T
    sort!(a); sort!(b)
    out = Vector{T}()
    i = 1; j = 1
    while i <= length(a) && j <= length(b)
        if a[i] < b[j]
            i += 1
        elseif b[j] < a[i]
            j += 1
        else
            push!(out, a[i]); i += 1; j += 1
        end
    end
    return out
end

# Memoized probability evaluation of a BDD.
function prob(M::MiniCUDD.Manager, node::MiniCUDD.BDDNode, pb::Dict{Int,Float64})
    memo = Dict{UInt, Float64}()
    oneptr = const1(M).ptr
    function go(n::MiniCUDD.BDDNode)
        key = UInt(n.ptr)
        if haskey(memo, key)
            return memo[key]
        end
        if isconstant(n)
            v = (n.ptr == oneptr) ? 1.0 : 0.0
            memo[key] = v
            return v
        end
        idx0 = Int(node_index(n))           # 0-based index in the variable order
        pvar = pb[idx0]
        t = bdd_then(M, n)
        e = bdd_else(M, n)
        val = pvar * go(t) + (1.0 - pvar) * go(e)
        memo[key] = val
        return val
    end
    return go(node)
end

# ===================== Per-triangle phi_1, phi_2 =====================
function disks_at_vertex(p::Point, circles::Vector{Circle}, vars::Dict{String,Int})
    idxs = Int[]
    for c in circles
        if isinside(c, p) && haskey(vars, c.name)
            push!(idxs, vars[c.name])
        end
    end
    return idxs
end

# phi_1(T) = AND over u in V(T) of OR over k in D(u) of x_k.
function phi1_of_triangle(M::MiniCUDD.Manager, t::Triangle3,
        circles::Vector{Circle}, vars::Dict{String,Int})
    ors = MiniCUDD.BDDNode[]
    for p in t.vertices
        push!(ors, or_of_vars(M, disks_at_vertex(p, circles, vars)))
    end
    return bdd_and(M, ors)
end

# phi_2(T) = OR over k in D(v1) ∩ D(v2) ∩ D(v3) of x_k.
function phi2_of_triangle(M::MiniCUDD.Manager, t::Triangle3,
        circles::Vector{Circle}, vars::Dict{String,Int})
    a = disks_at_vertex(v1(t), circles, vars)
    b = disks_at_vertex(v2(t), circles, vars)
    c = disks_at_vertex(v3(t), circles, vars)
    common = set_intersection_vec(copy(a), copy(b))
    common = set_intersection_vec(common, copy(c))
    return or_of_vars(M, common)
end

# ===================== Solver =====================
struct LevelSnapshot
    level::Int
    n_reached::Int
    n_resolved::Int
    n_refined::Int
    area_resolved::Float64
    nodes_phi1::Int
    nodes_phi2::Int
    peak_nodes::Int
    time_sec::Float64
    R_lower::Float64
    R_upper::Float64
end

struct SolverResult
    varphi1::MiniCUDD.BDDNode
    varphi2::MiniCUDD.BDDNode
    conv::Bool
    level::Int
    resolvedByLevel::Vector{Tuple{Int,Int}}   # (level, n_resolved)
    snapshots::Vector{LevelSnapshot}
    # Run totals (summed over all levels), for the separate-run-per-maxLevel
    # protocol. `resolved_equiv_total` counts triangles accepted by the local
    # or global equivalence test (Algorithm lines 8/14); `forced_total` counts
    # triangles committed only because they hit the depth limit (line 21).
    reached_total::Int
    resolved_equiv_total::Int
    forced_total::Int
    refined_total::Int
end

mutable struct LevelStat
    reached::Int
    resolved::Int    # all commits: equivalence-accepted + depth-limit forced
    refined::Int
    forced::Int      # subset of `resolved` committed only because Level == maxlevel
    area::Float64
end
LevelStat() = LevelStat(0, 0, 0, 0, 0.0)

@inline function triangle_area(t::Triangle3)
    a, b, c = v1(t), v2(t), v3(t)
    return 0.5 * abs(
        (getx(b) - getx(a)) * (gety(c) - gety(a)) -
        (getx(c) - getx(a)) * (gety(b) - gety(a)))
end

"""
    triangle_bdd_solver(M, vars, circles, init_triangles; maxlevel=10)

Adaptive symbolic refinement over a triangular partition. Mirrors
Algorithm 1 of the parent paper, with:

- vertex BDDs over 3 corners (not 4),
- 4-way edge-midpoint subdivision (not 2x2 grid).

Each initial triangle is pushed at level 1. When a triangle's local
phi_1 / phi_2 disagree and global merging into the running Phi1 / Phi2
also fails the equivalence test, the triangle is subdivided into its
four midpoint children and they are re-queued at level + 1.
"""
function triangle_bdd_solver(M::MiniCUDD.Manager, vars::Dict{String,Int},
        circles::Vector{Circle}, init_triangles::Vector{Triangle3};
        maxlevel::Int = 10, verbose::Bool = true,
        pb::Union{Nothing,Dict{Int,Float64}} = nothing)

    varphi1 = const1(M)
    varphi2 = const1(M)

    q = Vector{Tuple{Triangle3,Int}}()
    for t in init_triangles
        push!(q, (t, 1))
    end
    total = length(q); processed = 0
    level_stats = Dict{Int,LevelStat}()
    snapshots = LevelSnapshot[]
    level = 1
    prev_level = 0
    t0 = time()

    pm = verbose ?
        Progress(total; dt=0.3,
                 desc=@sprintf("Solving (Level %d, %d/%d)", level, processed, total)) :
        nothing
    last_total = total

    # Certified bounds at end-of-level-L cover the FULL partition at that
    # moment: resolved triangles (already merged into varphi1/varphi2) plus
    # pending triangles still in the queue (children of refined parents,
    # which together cover the same area as their parents). `extra` is the
    # triangle that was just popped at the level transition but not yet
    # processed; it must also be counted as pending.
    function take_snapshot!(L::Int, extra::Union{Nothing,Triangle3} = nothing)
        s = get(level_stats, L, LevelStat())
        if pb === nothing
            lb = NaN; ub = NaN
        else
            cur_phi1 = varphi1
            cur_phi2 = varphi2
            if extra !== nothing
                p1 = phi1_of_triangle(M, extra, circles, vars)
                p2 = phi2_of_triangle(M, extra, circles, vars)
                cur_phi1 = MiniCUDD.bdd_and(M, cur_phi1, p1)
                cur_phi2 = MiniCUDD.bdd_and(M, cur_phi2, p2)
            end
            for (tq, _) in q
                p1 = phi1_of_triangle(M, tq, circles, vars)
                p2 = phi2_of_triangle(M, tq, circles, vars)
                cur_phi1 = MiniCUDD.bdd_and(M, cur_phi1, p1)
                cur_phi2 = MiniCUDD.bdd_and(M, cur_phi2, p2)
            end
            lb = prob(M, cur_phi2, pb)
            ub = prob(M, cur_phi1, pb)
        end
        push!(snapshots, LevelSnapshot(
            L, s.reached, s.resolved, s.refined, s.area,
            Int(dag_size(varphi1)), Int(dag_size(varphi2)),
            Int(peak_node_count(M)), time() - t0, lb, ub))
    end

    while !isempty(q)
        t, curLevel = popfirst!(q)

        # Snapshot at level transition (end of previous level). The
        # just-popped triangle `t` belongs to curLevel and has not yet
        # been processed, so it counts as pending for prev_level.
        if curLevel != prev_level
            if prev_level > 0
                take_snapshot!(prev_level, t)
            end
            prev_level = curLevel
        end

        level = curLevel
        processed += 1

        p1 = phi1_of_triangle(M, t, circles, vars)
        p2 = phi2_of_triangle(M, t, circles, vars)

        v1d = MiniCUDD.bdd_and(M, varphi1, p1)
        v2d = MiniCUDD.bdd_and(M, varphi2, p2)
        same_local  = (p1.ptr  == p2.ptr)
        same_global = (v1d.ptr == v2d.ptr)

        s = get!(level_stats, curLevel, LevelStat())
        s.reached += 1

        if same_local || same_global || curLevel == maxlevel
            varphi1 = v1d
            varphi2 = v2d
            s.resolved += 1
            s.area    += triangle_area(t)
            if !(same_local || same_global)
                s.forced += 1          # committed at the depth limit (Alg. line 21)
            end
        else
            s.refined += 1
            for child in subdivide(t)
                push!(q, (child, curLevel + 1))
                total += 1
            end
        end

        if verbose
            if total != last_total
                pm.n = total
                last_total = total
            end
            update!(pm, processed;
                    desc=@sprintf("Solving (Level %d, %d/%d)",
                                  level, processed, total))
        end
    end
    verbose && finish!(pm)

    # final snapshot at the last level seen
    if prev_level > 0
        take_snapshot!(prev_level)
    end

    by_level = sort(collect(level_stats); by = first)
    pairs = [(Int(k), Int(v.resolved)) for (k,v) in by_level]
    reached_total = sum(v.reached  for (_, v) in level_stats; init = 0)
    commits_total = sum(v.resolved for (_, v) in level_stats; init = 0)
    forced_total  = sum(v.forced   for (_, v) in level_stats; init = 0)
    refined_total = sum(v.refined  for (_, v) in level_stats; init = 0)
    return SolverResult(varphi1, varphi2,
                        varphi1.ptr == varphi2.ptr, level, pairs, snapshots,
                        reached_total, commits_total - forced_total,
                        forced_total, refined_total)
end

# ===================== JSON I/O =====================
function load_sensors_json(path::String)
    cfg = JSON.parse(read(path, String))
    reliability = haskey(cfg, "reliability") ? Float64(cfg["reliability"]) : 0.9
    maxlevel    = haskey(cfg, "maxlevel")    ? Int(cfg["maxlevel"])        : 10
    circles = Circle[]
    for c in cfg["circles"]
        push!(circles, Circle(
            string(c["name"]),
            Point(Float64(c["x"]), Float64(c["y"])),
            Float64(c["radius"]),
        ))
    end
    return circles, reliability, maxlevel, cfg
end

function load_triangulation_json(path::String)
    d = JSON.parse(read(path, String))
    V = [Point(Float64(p[1]), Float64(p[2])) for p in d["vertices"]]
    tris = Triangle3[]
    for t in d["triangles"]
        @inbounds push!(tris, Triangle3((V[t[1]+1], V[t[2]+1], V[t[3]+1])))
    end
    return tris, V, d
end

# ===================== Square-area smoke test =====================
"""
    smoke_test_square(; case = :four_corners, maxlevel = 8) -> NamedTuple

Self-contained sanity check on the unit square split into 2 triangles
along the y = x diagonal.

`case`:
  - :one_big       — single disk of radius 0.8 at (0.5, 0.5).
                     Covers both triangles; Phi1 and Phi2 should match
                     immediately (no refinement).
  - :four_corners  — four disks of radius 0.6 at the corners.
                     Bounds should be non-trivial; refinement happens.
  - :two_disjoint  — two disks of radius 0.3 at (0.2, 0.5), (0.8, 0.5).
                     Coverage is impossible (region is not covered by
                     any single disk choice); R_lower should be 0 but
                     R_upper > 0; should stop at maxlevel.
  - :two_bands     — two disks of radius 0.8 at (0.5, 0) and (0.5, 1).
                     Together they cover the whole unit square.
                     phi_1 / phi_2 disagree at level 1 but converge to
                     (x_1 AND x_2) at level 2 → R = p_k^2 = 0.81.
"""
function smoke_test_square(; case::Symbol = :four_corners, maxlevel::Int = 8)
    ll = Point(0.0, 0.0); lr = Point(1.0, 0.0)
    ur = Point(1.0, 1.0); ul = Point(0.0, 1.0)
    init_triangles = [
        Triangle3((ll, lr, ur)),
        Triangle3((ll, ur, ul)),
    ]

    if case === :one_big
        circles = [Circle("1", Point(0.5, 0.5), 0.8)]
    elseif case === :four_corners
        r = 0.6
        circles = [
            Circle("1", Point(0.0, 0.0), r),
            Circle("2", Point(1.0, 0.0), r),
            Circle("3", Point(1.0, 1.0), r),
            Circle("4", Point(0.0, 1.0), r),
        ]
    elseif case === :two_disjoint
        circles = [
            Circle("1", Point(0.2, 0.5), 0.3),
            Circle("2", Point(0.8, 0.5), 0.3),
        ]
    elseif case === :two_bands
        circles = [
            Circle("1", Point(0.5, 0.0), 0.8),
            Circle("2", Point(0.5, 1.0), 0.8),
        ]
    else
        error("unknown case: $case")
    end

    p_k = 0.9
    vars = Dict{String,Int}()
    pb = Dict{Int,Float64}()
    for (i,c) in enumerate(circles)
        vars[c.name] = i - 1  # 0-based
        pb[i - 1] = p_k
    end

    M = Manager(nvars = max(1, length(circles)))
    t0 = @elapsed begin
        res = triangle_bdd_solver(M, vars, circles, init_triangles;
                                  maxlevel = maxlevel, verbose = false)
    end
    lb = prob(M, res.varphi2, pb)
    ub = prob(M, res.varphi1, pb)
    n_initial = length(init_triangles)
    n_resolved = sum(v for (_,v) in res.resolvedByLevel)
    nodes_phi1 = dag_size(res.varphi1)
    nodes_phi2 = dag_size(res.varphi2)

    out = (
        case         = case,
        n_circles    = length(circles),
        n_initial    = n_initial,
        n_resolved   = n_resolved,
        max_level    = res.level,
        resolvedByLevel = res.resolvedByLevel,
        converged    = res.conv,
        R_lower      = lb,
        R_upper      = ub,
        nodes_phi1   = nodes_phi1,
        nodes_phi2   = nodes_phi2,
        time_sec     = t0,
    )

    quit(M)
    return out
end

function print_smoke_result(r)
    println("------------------------------------------------------------")
    @printf("case            = %s\n", r.case)
    @printf("n_circles       = %d\n", r.n_circles)
    @printf("init triangles  = %d\n", r.n_initial)
    @printf("resolved        = %d\n", r.n_resolved)
    @printf("max level seen  = %d\n", r.max_level)
    @printf("resolved/level  = %s\n", r.resolvedByLevel)
    @printf("converged       = %s\n", r.converged)
    @printf("R_lower (Phi2)  = %.12f\n", r.R_lower)
    @printf("R_upper (Phi1)  = %.12f\n", r.R_upper)
    @printf("nodes phi1/phi2 = %d / %d\n", r.nodes_phi1, r.nodes_phi2)
    @printf("time (sec)      = %.4f\n", r.time_sec)
end

# ===================== CLI entry point for the smoke test =====================
function main(argv)
    cases = isempty(argv) ?
        [:one_big, :four_corners, :two_disjoint, :two_bands] :
        [Symbol(a) for a in argv]
    for c in cases
        r = smoke_test_square(case = c, maxlevel = 8)
        print_smoke_result(r)
    end
end

end  # module

if abspath(PROGRAM_FILE) == @__FILE__
    BDDSolverTriangle.main(ARGS)
end
