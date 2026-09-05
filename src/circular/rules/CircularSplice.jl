# circular/rules/CircularSplice.jl — SPLICE a morphism `f : w → w`
# into a cut curve of `k` edges. It only INSERTS, it never cuts a window out
# (the collapse barrier concerns the cut-out surgery, not this).
#
# THE PICTURE. The `k` edges (the `parent_edge` chain of a region `R`,
# colours = the region word `w`) form a cut curve. Each edge is cut
# once; into the resulting band the morphism `f` is inserted. `f` comes as a
# CircularMorphismGraph — the CUTS carry which leaves are bottom and which are
# top (⚠ `CircularWord` normalises to a canonical ROTATION, so bare leaf
# numbers 1..2k mean nothing).
#
# ORIENTATION IS NOT GUESSED (house style, the `_fds_braid2`/`expand_gbraid`
# precedent): which end of a cut edge attaches to the bottom and which to the
# top side of `f` is a choice PER EDGE (2^k assignments). All are tried in a
# fixed order and the FIRST one that yields `euler == 2`, an empty
# `check_wiring` and the expected region count wins. `nothing` if none does;
# whether the cut curve is always planar-insertable is open — a hit rate
# below 100% is not itself an error.
#
# EXPECTED REGION COUNT. Splicing the IDENTITY on `w` (k parallel strands,
# `k + 1` regions as a circular graph) must reproduce `g` exactly — so a general
# insert has to end at `region_count(g) + region_count(f) - (k + 1)`.
#
# LABELS: version 1 is label-free — this function works on the bare
# CircularGraph; callers with non-trivial region labels must refuse first.

"""
    circular_splice_identity(w::Vector{Int}) -> CircularMorphismGraph

The identity on `w` as a splice-ready morphism (`k` parallel strands) — the
neutral element of [`circular_splice`](@ref). The cuts compensate the canonical
rotation `CircularWord` applies (file head).
"""
function circular_splice_identity(w::Vector{Int})
    k = length(w)
    raw = vcat(w, reverse(w))
    cw = CircularWord(raw)
    n = length(raw)
    n == 0 && return CircularMorphismGraph(CircularGraph(cw, CircularNode[], Edge[]), 0, 0)
    stored = letters(cw)
    r = findfirst(r -> all(stored[i] == raw[mod1(i + r, n)] for i in 1:n), 0:(n - 1)) - 1
    σ(p) = mod1(p - r, n)              # raw position -> stored position
    g = CircularGraph(cw, CircularNode[],
                      Edge[Edge(w[j], Leaf(σ(j)), Leaf(σ(2k + 1 - j))) for j in 1:k])
    return CircularMorphismGraph(g, mod(σ(1) - 1, n), mod(σ(k), n))
end

"""
    circular_splice(g::CircularGraph, cut::Vector{Int}, fm::CircularMorphismGraph;
                    sides = nothing) -> Union{Nothing, NamedTuple}

Cuts the edges `cut = e_1..e_k` of `g` once each and splices `fm` (a `w → w`
morphism, `w` = the colours of the cut) into the band. Returns `(graph,
sides)` — `sides[j] = true` means end `b` of `e_j` sits on the BOTTOM side —
or `nothing` if no assignment passes the three checks (`euler == 2`, empty
`check_wiring`, expected region count).

`sides` PINS the assignment instead of searching it: several morphisms spliced
into the SAME curve must be inserted the same way round, or the terms of one
relation would not describe one and the same local replacement (anchors the
assignment on the identity, which has to reproduce `g`). A pinned assignment
is still checked, not trusted — `nothing` if it fails.

Throws on malformed input (wrong boundary words, multi-hit leaves): those are
caller errors, not measurement outcomes.
"""
function circular_splice(g::CircularGraph, cut::Vector{Int}, fm::CircularMorphismGraph;
                         sides::Union{Nothing,AbstractVector{Bool}} = nothing)
    k = length(cut)
    k == 0 && throw(ArgumentError("circular_splice: empty cut curve"))
    allunique(cut) || throw(ArgumentError("circular_splice: cut edges must be distinct"))
    all(1 <= e <= length(g.edges) for e in cut) ||
        throw(ArgumentError("circular_splice: cut edge index out of range"))
    w = [g.edges[e].colour for e in cut]
    (bottom(fm) == w && top(fm) == w) || throw(ArgumentError(
        "circular_splice: f must be a morphism $w -> $w, is $(bottom(fm)) -> $(top(fm))"))
    f = fm.graph
    bl = _bottom_leaves(fm)                       # bottom leaf of strand j: bl[j]
    tl = _top_leaves(fm)                          # top word j sits at tl[k - j + 1]
    # every boundary leaf of f must be hit by exactly one edge end
    deg = Dict(i => 0 for i in vcat(bl, tl))
    for e in f.edges, p in (e.a, e.b)
        p isa Leaf && (deg[p.k] += 1)
    end
    all(==(1), values(deg)) || throw(ArgumentError(
        "circular_splice: every boundary leaf of f needs exactly one edge end"))
    strand = Dict{Int,Tuple{Symbol,Int}}()        # leaf -> (:bot/:top, j)
    for (j, lf) in enumerate(bl); strand[lf] = (:bot, j); end
    for (t, lf) in enumerate(tl); strand[lf] = (:top, k - t + 1); end

    basis = length(g.nodes)
    nodes = vcat(g.nodes, f.nodes)
    keep  = [e for (i, e) in enumerate(g.edges) if !(i in cut)]
    expected = region_count(g) + region_count(f) - (k + 1)

    masks = sides === nothing ? collect(0:(2^k - 1)) :
            (length(sides) == k ? [sum(sides[j] ? 1 << (j - 1) : 0 for j in 1:k)] :
             throw(ArgumentError("circular_splice: `sides` must have $k entries")))

    for mask in masks
        side = [((mask >> (j - 1)) & 1) == 1 for j in 1:k]
        bot = Dict{Int,Any}(); top_ = Dict{Int,Any}()
        for j in 1:k
            e = g.edges[cut[j]]
            bot[j]  = side[j] ? e.b : e.a
            top_[j] = side[j] ? e.a : e.b
        end
        remap(p) = p isa Leaf ?
                   (strand[p.k][1] === :bot ? bot[strand[p.k][2]] : top_[strand[p.k][2]]) :
                   NodePort(p.node + basis, p.slot)
        cand = try
            edges = vcat(keep, Edge[Edge(e.colour, remap(e.a), remap(e.b)) for e in f.edges])
            CircularGraph(g.word, nodes, edges)
        catch
            continue                       # cell trace refused — not this assignment
        end
        (euler(cand) == 2 && isempty(check_wiring(cand)) &&
         region_count(cand) == expected) || continue
        return (graph = cand, sides = side)
    end
    return nothing
end

"""
    circular_splice_sides(g::CircularGraph, cut::Vector{Int}, fm::CircularMorphismGraph)
        -> Vector{Vector{Bool}}

ALL assignments under which `fm` splices into `cut` — the same three checks as
[`circular_splice`](@ref), in mask order; `circular_splice` returns the first of
them.

This is what pins the assignment across the terms of ONE relation: the caller
intersects these sets over all terms and takes the first survivor
(`_crf_common_sides`, circular/rules/CircularRexFusion.jl). The identity splices
under all `2^k` of them, which is exactly why it cannot serve as the anchor.
"""
function circular_splice_sides(g::CircularGraph, cut::Vector{Int},
                               fm::CircularMorphismGraph)
    k = length(cut)
    out = Vector{Bool}[]
    for mask in 0:(2^k - 1)
        side = [((mask >> (j - 1)) & 1) == 1 for j in 1:k]
        circular_splice(g, cut, fm; sides = side) === nothing || push!(out, side)
    end
    return out
end
