module BDDAnalysis

# BDD-based reliability analysis on the MiniCUDD backend.
#
# This module ports the analysis layer of FaultTree.jl (which uses the DD.jl
# backend) to the MiniCUDD (CUDD) BDDs produced by the coverage solvers
# `bdd_solver.jl` and `bdd_solver_triangle.jl`. It operates directly on a
# `(MiniCUDD.Manager, MiniCUDD.BDDNode)` pair, so it is independent of which
# solver built the structure function.
#
# Provided analyses:
#   * minsol           — minimal-solutions BDD (the BDD whose 1-paths are the
#                        minimal solutions of a coherent function).
#   * min_path_sets    — minimal sensor sets that guarantee coverage.
#   * min_cut_sets     — minimal sensor-failure sets that destroy coverage
#                        (minimal solutions of the dual function).
#   * birnbaum         — Birnbaum importance dR/dp_k of every component.
#   * criticality1/0   — criticality importances.
#   * structure_importance — Birnbaum importance at p = 1/2.
#
# The algorithms are functional recurrences over the Shannon decomposition,
# so they are correct regardless of CUDD's complement-edge representation:
# `bdd_then` / `bdd_else` return true functional cofactors, and `bdd_ite`
# rebuilds nodes functionally. Nodes are identified by `UInt(n.ptr)`, which
# carries the complement bit, so a function and its negation are distinguished.

using MiniCUDD

export minsol, dual, enumerate_solutions, count_solutions, min_cardinality
export min_path_sets, min_cut_sets, min_path_sets_named, min_cut_sets_named
export prob_value, birnbaum, importances, birnbaum_named, structure_importance
export smoke_test

@inline _id(n::MiniCUDD.BDDNode) = UInt(n.ptr)

# ===================== minimal solutions =====================

# _without(f, g): the function whose solutions are the solutions of `f` that
# are NOT supersets of any solution of `g` (set subtraction modulo
# subsumption). Ported from FaultTree.jl/src/_without.jl. The DD.jl `level`
# test (root has the larger level) becomes a `node_index` test (root has the
# smaller index), so the comparison direction is flipped.
function _without(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode, g::MiniCUDD.BDDNode,
        one_ptr::UInt, cache::Dict{Tuple{UInt,UInt},MiniCUDD.BDDNode})
    fc = isconstant(f)
    gc = isconstant(g)
    if fc && gc
        if _id(f) != one_ptr            # f == 0
            return const0(M)
        elseif _id(g) == one_ptr        # f == 1, g == 1
            return const0(M)
        else                            # f == 1, g == 0
            return f
        end
    elseif fc                           # f terminal, g non-terminal
        return f
    elseif gc                           # f non-terminal, g terminal
        return (_id(g) == one_ptr) ? const0(M) : f
    end
    key = (_id(f), _id(g))
    haskey(cache, key) && return cache[key]
    idf = Int(node_index(f))
    idg = Int(node_index(g))
    if idf < idg                        # f's variable is nearer the root
        n0 = _without(M, bdd_else(M, f), g, one_ptr, cache)
        n1 = _without(M, bdd_then(M, f), g, one_ptr, cache)
        res = bdd_ite(M, var(M, idf), n1, n0)
    elseif idf > idg                    # g's variable is nearer the root
        res = _without(M, f, bdd_else(M, g), one_ptr, cache)
    else                                # same variable
        n0 = _without(M, bdd_else(M, f), bdd_else(M, g), one_ptr, cache)
        n1 = _without(M, bdd_then(M, f), bdd_then(M, g), one_ptr, cache)
        res = bdd_ite(M, var(M, idf), n1, n0)
    end
    cache[key] = res
    return res
end

function _minsol(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode, one_ptr::UInt,
        c1::Dict{UInt,MiniCUDD.BDDNode}, c2::Dict{Tuple{UInt,UInt},MiniCUDD.BDDNode})
    isconstant(f) && return f
    key = _id(f)
    haskey(c1, key) && return c1[key]
    g = bdd_then(M, f)
    h = bdd_else(M, f)
    k = _minsol(M, g, one_ptr, c1, c2)
    u = _without(M, k, h, one_ptr, c2)          # then-child: minsol(g) \ h
    v = _minsol(M, h, one_ptr, c1, c2)          # else-child: minsol(h)
    res = bdd_ite(M, var(M, Int(node_index(f))), u, v)
    c1[key] = res
    return res
end

"""
    minsol(M, f) -> BDDNode

Return the minimal-solutions BDD of `f`. For a coherent (monotone) `f`, the
1-paths of the result are exactly the minimal solutions of `f`. Ported from
`minsol` in FaultTree.jl/src/_without.jl.
"""
function minsol(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode)
    one_ptr = _id(const1(M))
    c1 = Dict{UInt,MiniCUDD.BDDNode}()
    c2 = Dict{Tuple{UInt,UInt},MiniCUDD.BDDNode}()
    return _minsol(M, f, one_ptr, c1, c2)
end

"""
    enumerate_solutions(M, g; limit = typemax(Int)) -> Vector{Vector{Int}}

Enumerate the 1-paths of a minimal-solutions BDD `g`, each as the sorted set
of variable indices taken on the `then` branch. Stops after `limit` sets.
Ported from `extractpath` in FaultTree.jl/src/_mcs.jl.
"""
function enumerate_solutions(M::MiniCUDD.Manager, g::MiniCUDD.BDDNode;
        limit::Int = typemax(Int))
    one_ptr = _id(const1(M))
    out = Vector{Vector{Int}}()
    path = Int[]
    function rec(n::MiniCUDD.BDDNode)
        length(out) >= limit && return
        if isconstant(n)
            (_id(n) == one_ptr) && push!(out, sort(copy(path)))
            return
        end
        idx = Int(node_index(n))
        push!(path, idx)
        rec(bdd_then(M, n))
        pop!(path)
        rec(bdd_else(M, n))
    end
    rec(g)
    return out
end

"""
    count_solutions(M, g) -> Int

Number of 1-paths of `g` (the number of minimal solutions when `g` is a
minimal-solutions BDD). Memoized DAG traversal, so it does not enumerate.
"""
function count_solutions(M::MiniCUDD.Manager, g::MiniCUDD.BDDNode)
    one_ptr = _id(const1(M))
    memo = Dict{UInt,Int}()
    function go(n::MiniCUDD.BDDNode)
        k = _id(n)
        haskey(memo, k) && return memo[k]
        v = isconstant(n) ? (k == one_ptr ? 1 : 0) :
            go(bdd_then(M, n)) + go(bdd_else(M, n))
        memo[k] = v
        return v
    end
    return go(g)
end

"""
    min_cardinality(M, g) -> Int

Fewest `then` edges on any 1-path of `g` (the smallest minimal-solution size
when `g` is a minimal-solutions BDD). Returns `-1` if `g` has no 1-path.
"""
function min_cardinality(M::MiniCUDD.Manager, g::MiniCUDD.BDDNode)
    one_ptr = _id(const1(M))
    INF = typemax(Int) ÷ 2
    memo = Dict{UInt,Int}()
    function go(n::MiniCUDD.BDDNode)
        k = _id(n)
        haskey(memo, k) && return memo[k]
        v = isconstant(n) ? (k == one_ptr ? 0 : INF) :
            min(1 + go(bdd_then(M, n)), go(bdd_else(M, n)))
        memo[k] = v
        return v
    end
    r = go(g)
    return r >= INF ? -1 : r
end

# ===================== duality (min-cut sets) =====================

_not(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode) = bdd_xor(M, f, const1(M))

# negate_inputs(f) = f(not x): swap then/else at every node.
function _negate_inputs(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode,
        cache::Dict{UInt,MiniCUDD.BDDNode})
    isconstant(f) && return f
    key = _id(f)
    haskey(cache, key) && return cache[key]
    t = _negate_inputs(M, bdd_then(M, f), cache)
    e = _negate_inputs(M, bdd_else(M, f), cache)
    res = bdd_ite(M, var(M, Int(node_index(f))), e, t)   # children swapped
    cache[key] = res
    return res
end

"""
    dual(M, f) -> BDDNode

The dual function `f^d(x) = not f(not x)`. For a coherent coverage function
`f`, the minimal solutions of `dual(f)` are the minimal cut sets of `f`.
"""
function dual(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode)
    g = _negate_inputs(M, f, Dict{UInt,MiniCUDD.BDDNode}())
    return _not(M, g)
end

# ===================== named set extraction =====================

_invert(vars::Dict{String,Int}) = Dict{Int,String}(v => k for (k, v) in vars)

"""
    min_path_sets(M, f; limit) -> Vector{Vector{Int}}

Minimal solutions of `f` (for the coverage function, the minimal sensor sets
that guarantee coverage).
"""
min_path_sets(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode; limit::Int = typemax(Int)) =
    enumerate_solutions(M, minsol(M, f); limit = limit)

"""
    min_cut_sets(M, f; limit) -> Vector{Vector{Int}}

Minimal solutions of `dual(f)` (for the coverage function, the minimal
sensor-failure sets that destroy coverage).
"""
min_cut_sets(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode; limit::Int = typemax(Int)) =
    enumerate_solutions(M, minsol(M, dual(M, f)); limit = limit)

function min_path_sets_named(M, f, vars::Dict{String,Int}; limit::Int = typemax(Int))
    inv = _invert(vars)
    [sort([inv[i] for i in s]) for s in min_path_sets(M, f; limit = limit)]
end

function min_cut_sets_named(M, f, vars::Dict{String,Int}; limit::Int = typemax(Int))
    inv = _invert(vars)
    [sort([inv[i] for i in s]) for s in min_cut_sets(M, f; limit = limit)]
end

# ===================== probability and importance =====================

"""
    prob_table(M, f, pb) -> Dict{UInt,Float64}

Memoized probability of every node reachable from `f`, keyed by `UInt(ptr)`.
`pb` maps a variable index to its component reliability.
"""
function prob_table(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode, pb::Dict{Int,Float64})
    one_ptr = _id(const1(M))
    memo = Dict{UInt,Float64}()
    function go(n::MiniCUDD.BDDNode)
        k = _id(n)
        haskey(memo, k) && return memo[k]
        if isconstant(n)
            v = (k == one_ptr) ? 1.0 : 0.0
            memo[k] = v
            return v
        end
        p = pb[Int(node_index(n))]
        v = p * go(bdd_then(M, n)) + (1.0 - p) * go(bdd_else(M, n))
        memo[k] = v
        return v
    end
    go(f)
    return memo
end

"""
    prob_value(M, f, pb) -> Float64

Probability that `f` evaluates to 1 under component reliabilities `pb`.
"""
prob_value(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode, pb::Dict{Int,Float64}) =
    prob_table(M, f, pb)[_id(f)]

# Topological order of the BDD DAG, root first, every parent before its
# children. Ported from FaultTree.jl/src/_tsort.jl.
function _tsort(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode)
    state = Dict{UInt,Int}()
    result = MiniCUDD.BDDNode[]
    function visit(n::MiniCUDD.BDDNode)
        k = _id(n)
        s = get(state, k, 0)
        s == 2 && return
        s == 1 && error("BDD DAG has a cycle")
        state[k] = 1
        if !isconstant(n)
            visit(bdd_then(M, n))
            visit(bdd_else(M, n))
        end
        state[k] = 2
        pushfirst!(result, n)
    end
    visit(f)
    return result
end

"""
    birnbaum(M, f, pb) -> Dict{Int,Float64}

Birnbaum importance dPr(f=1)/dp_k of every component, computed in one adjoint
pass over the BDD. Ported from `grad` in FaultTree.jl/src/_grad.jl. Components
on which `f` does not depend are absent from the result (importance 0).
"""
function birnbaum(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode, pb::Dict{Int,Float64})
    pr = prob_table(M, f, pb)
    weight = Dict{UInt,Float64}(_id(f) => 1.0)
    imp = Dict{Int,Float64}()
    for n in _tsort(M, f)
        isconstant(n) && continue
        w = get(weight, _id(n), 0.0)
        idx = Int(node_index(n))
        p = pb[idx]
        t = bdd_then(M, n)
        e = bdd_else(M, n)
        kt = _id(t)
        ke = _id(e)
        weight[kt] = get(weight, kt, 0.0) + w * p
        weight[ke] = get(weight, ke, 0.0) + w * (1.0 - p)
        imp[idx] = get(imp, idx, 0.0) + w * (pr[kt] - pr[ke])
    end
    return imp
end

"""
    importances(M, f, pb) -> NamedTuple

Reliability `R`, and the Birnbaum, criticality-1 and criticality-0 importance
of every component variable index, all from the one converged BDD `f`.
"""
function importances(M::MiniCUDD.Manager, f::MiniCUDD.BDDNode, pb::Dict{Int,Float64})
    pr = prob_table(M, f, pb)
    R = pr[_id(f)]
    weight = Dict{UInt,Float64}(_id(f) => 1.0)
    bi = Dict{Int,Float64}()
    for n in _tsort(M, f)
        isconstant(n) && continue
        w = get(weight, _id(n), 0.0)
        idx = Int(node_index(n))
        p = pb[idx]
        t = bdd_then(M, n)
        e = bdd_else(M, n)
        kt = _id(t)
        ke = _id(e)
        weight[kt] = get(weight, kt, 0.0) + w * p
        weight[ke] = get(weight, ke, 0.0) + w * (1.0 - p)
        bi[idx] = get(bi, idx, 0.0) + w * (pr[kt] - pr[ke])
    end
    c1 = Dict{Int,Float64}(k => v * pb[k] / R for (k, v) in bi)
    c0 = Dict{Int,Float64}(k => v * (1.0 - pb[k]) / (1.0 - R) for (k, v) in bi)
    return (R = R, birnbaum = bi, crit1 = c1, crit0 = c0)
end

"""
    birnbaum_named(M, f, pb, vars) -> Dict{String,Float64}

Birnbaum importance keyed by component name; names absent from the BDD get 0.
"""
function birnbaum_named(M, f, pb::Dict{Int,Float64}, vars::Dict{String,Int})
    bi = birnbaum(M, f, pb)
    Dict{String,Float64}(name => get(bi, idx, 0.0) for (name, idx) in vars)
end

"""
    structure_importance(M, f, vars) -> Dict{String,Float64}

Structure importance (Birnbaum importance evaluated at p = 1/2 for every
component), keyed by component name. Ported from `smeas`.
"""
function structure_importance(M, f, vars::Dict{String,Int})
    pb = Dict{Int,Float64}(idx => 0.5 for idx in values(vars))
    birnbaum_named(M, f, pb, vars)
end

# ===================== smoke test =====================

# Compare two collections of index sets as sets-of-sets.
function _sameset(a::Vector{Vector{Int}}, b::Vector{Vector{Int}})
    norm(x) = Set(Tuple(sort(s)) for s in x)
    return norm(a) == norm(b)
end

"""
    smoke_test() -> Bool

Self-contained correctness checks on hand-computable BDDs. Returns `true` if
all checks pass; prints a report.
"""
function smoke_test()
    ok = true
    report(name, pass) = (println("  ", pass ? "PASS" : "FAIL", "  ", name);
                          ok &= pass)

    M = Manager(nvars = 3)
    x0 = var(M, 0)
    x1 = var(M, 1)
    x2 = var(M, 2)
    p = Dict{Int,Float64}(0 => 0.9, 1 => 0.8, 2 => 0.7)

    # --- f = x0 AND x1 ---
    fand = bdd_and(M, x0, x1)
    report("AND: min_path_sets == {{0,1}}",
           _sameset(min_path_sets(M, fand), [[0, 1]]))
    report("AND: min_cut_sets == {{0},{1}}",
           _sameset(min_cut_sets(M, fand), [[0], [1]]))
    bi = birnbaum(M, fand, p)
    report("AND: Birnbaum(x0) == p1", isapprox(bi[0], p[1]; atol = 1e-12))
    report("AND: Birnbaum(x1) == p0", isapprox(bi[1], p[0]; atol = 1e-12))

    # --- f = x0 OR x1 ---
    forr = bdd_or(M, x0, x1)
    report("OR: min_path_sets == {{0},{1}}",
           _sameset(min_path_sets(M, forr), [[0], [1]]))
    report("OR: min_cut_sets == {{0,1}}",
           _sameset(min_cut_sets(M, forr), [[0, 1]]))
    bi = birnbaum(M, forr, p)
    report("OR: Birnbaum(x0) == 1 - p1", isapprox(bi[0], 1.0 - p[1]; atol = 1e-12))

    # --- f = 2-out-of-3 ---
    f2 = bdd_or(M, bdd_or(M, bdd_and(M, x0, x1), bdd_and(M, x0, x2)),
                bdd_and(M, x1, x2))
    report("2oo3: min_path_sets == {{0,1},{0,2},{1,2}}",
           _sameset(min_path_sets(M, f2), [[0, 1], [0, 2], [1, 2]]))
    report("2oo3: min_cut_sets == {{0,1},{0,2},{1,2}}",
           _sameset(min_cut_sets(M, f2), [[0, 1], [0, 2], [1, 2]]))
    report("2oo3: min_cardinality == 2",
           min_cardinality(M, minsol(M, f2)) == 2)
    report("2oo3: count_solutions == 3",
           count_solutions(M, minsol(M, f2)) == 3)
    bi = birnbaum(M, f2, p)
    # Birnbaum(x0) of 2oo3 = Pr(x1 OR x2) - Pr(x1 AND x2) = p1 + p2 - 2 p1 p2
    b0 = p[1] + p[2] - 2 * p[1] * p[2]
    report("2oo3: Birnbaum(x0) closed form", isapprox(bi[0], b0; atol = 1e-12))

    # --- Birnbaum cross-check by central finite difference ---
    eps = 1e-6
    pp = copy(p); pp[1] += eps
    pm = copy(p); pm[1] -= eps
    fd = (prob_value(M, f2, pp) - prob_value(M, f2, pm)) / (2eps)
    report("2oo3: Birnbaum(x1) matches finite difference",
           isapprox(bi[1], fd; atol = 1e-6))

    # --- minsol soundness: OR-of-ANDs of min-path sets re-equals f ---
    sets = min_path_sets(M, f2)
    acc = const0(M)
    for s in sets
        term = const1(M)
        for i in s
            term = bdd_and(M, term, var(M, i))
        end
        acc = bdd_or(M, acc, term)
    end
    report("2oo3: OR-of-ANDs of min-path sets re-equals f",
           UInt(acc.ptr) == UInt(f2.ptr))

    quit(M)
    println(ok ? "ALL PASS" : "FAILURES PRESENT")
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    smoke_test()
end

end  # module
