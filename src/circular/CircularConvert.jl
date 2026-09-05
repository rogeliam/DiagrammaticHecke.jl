# circular/CircularConvert.jl — WordGraph → CircularGraph.
#
# One plain node ↦ one new node, SLOT BY SLOT: `arms[s]` is the colour of the
# plain slot `s`, read via `_slot_colour`. This keeps `word`, `edges` and hence ALL
# ports (Leaf/NodePort/Circle) and node indices literally unchanged — the
# rotation system is identical, and so are the cells. That is exactly what
# makes the converter trivially correct and turns the test
# `face_count(circular(g)) == face_count(g)` into a genuine check rather than a
# tautology.
#
# `circular` DOES NOT MERGE: merging is a reduction rule, not a constructor
# obligation — unmerged CircularGraphs must stay representable.

"The CircularNode kind a single colour belongs to: 2 ⇒ `:mono`, 1/3 ⇒ `:mixed`."
_circular_kind_of_colour(c::Int) =
    c == 2 ? :mono :
    (c == 1 || c == 3) ? :mixed :
    throw(ArgumentError("colour $c belongs to no CircularNode kind (expected 1, 2 or 3)"))

"""
    circular(nd::Node) -> CircularNode

Translates a `Node` into the open node type:

| plain | new |
|---|---|
| `dot(c)` | `[c]` |
| `trivalent(c)` | `[c,c,c]` |
| `braid(1,3; m=2)` | `CircularNode(:mixed, [1,3,1,3])` |
| `braid(s,t; m=3)` | `CircularNode(:braid, [s,t,s,t,s,t])` |
| `braid(s,t; m)` , m ≥ 3 | `CircularNode(:braid, 2m arms)` — e.g. m = 5, a 10-armed braid node |

The m=2 commutation is NOT `:braid`: it is just a 4-armed node of
the 1/3 world, and that identification is one of the reasons for the new type
(1/3 combinations that the plain kinds could not see as isomorphic).
"""
function circular(nd::Node)
    deg = _node_degree(nd)
    cols = [_slot_colour(nd, s) for s in 1:deg]
    if nd.kind === :braid
        if nd.m >= 3
            # m ≥ 3 ⇒ a `:braid` node with 2m arms. The arm count is an
            # ATTRIBUTE of `:braid` (2k arms, k ≥ 3), so m = 5 (a 10-armed
            # braid node) needs no new type — just this branch, and lets such
            # diagrams use the ORDINARY renderer instead of a dedicated one.
            # Nothing about A₃ changes: m = 3 takes this very path, m = 2
            # keeps its own branch below (`_circular_kind_of_colour` is
            # untouched — it is only consulted for degree < 3 nodes, i.e.
            # dots).
            return _unchecked_circularnode(:braid, cols)
        elseif nd.m == 2
            # The only m=2 commutation in A₃ is {1,3}. Anything else aborts
            # LOUDLY here instead of silently slipping through as :mixed.
            Set(nd.colours) == Set((1, 3)) || throw(ArgumentError(
                "circular: m=2 braid on colours $(nd.colours) is not in the 1/3 world"))
            return _unchecked_circularnode(:mixed, cols)
        else
            throw(ArgumentError("circular: braid with m = $(nd.m) is not modelled"))
        end
    end
    # For low-degree plain nodes (dots/beads) construct unchecked CircularNodes so
    # that the converter never fails; the reduction rules handle their removal.
    if deg < 3
        return _unchecked_circularnode(_circular_kind_of_colour(cols[1]), cols)
    end
    return circular_node(cols)
end

"""
    circular(g::WordGraph) -> CircularGraph

Converts a whole diagram. `word` and `edges` are taken over VERBATIM, only the
node objects are replaced — ports and node indices thus stay valid.
Consequences pinned down by `test/circular.jl`: same cell count, same
`gap_cell`/`sector_cell`, same boundary word.
"""
circular(g::WordGraph) = CircularGraph(g.word, [circular(nd) for nd in g.nodes], copy(g.edges))

"""
    circular(m::MorphismGraph) -> CircularGraph

Converts the underlying diagram. The cuts are lost in the process — use
[`circular_morphism`](@ref) if `cut1`/`cut2` should be preserved (e.g. for the
bottom/top boundary marking in the renderer).
"""
circular(m::MorphismGraph) = circular(m.graph)

# ---- CircularMorphismGraph: the bottom/top view for CircularGraph -----------------------
#
# Literal analogue of `MorphismGraph` (morphism/MorphismGraph.jl): the same two
# cut fields `cut1`/`cut2`, just with `CircularGraph` instead of `WordGraph`. The cut
# logic itself (`_bottom_leaves`/`_top_leaves`/`_left_mark`/`_right_mark`)
# depends only on `word`/`cut1`/`cut2` — never on the node type — and is
# therefore NOT duplicated below, but reused via a tiny adapter layer.

"""
    CircularMorphismGraph

A `CircularGraph` with two cuts `cut1`/`cut2` (same gap convention as
[`MorphismGraph`](@ref)) — turns the diagram into a bottom→top view, so the
renderer can draw the same black/grey boundary marking as for an ordinary
`MorphismGraph`.
"""
struct CircularMorphismGraph
    graph::CircularGraph
    cut1::Int
    cut2::Int
end

"Bottom-and-top view of a finished `CircularGraph`: the whole boundary word as bottom, ε as top."
circular_morphism_graph(g::CircularGraph) = CircularMorphismGraph(g, 0, 0)
circular_morphism_graph(g::CircularGraph, cut1::Int, cut2::Int) = CircularMorphismGraph(g, cut1, cut2)

# Regions under a solid boundary (the second cell notion, see diagram/Faces.jl):
# the cuts don't change that, so just pass through — as for MorphismGraph.
region_count(fm::CircularMorphismGraph) = region_count(fm.graph)
boundary_regions(fm::CircularMorphismGraph) = boundary_regions(fm.graph)
regions(fm::CircularMorphismGraph) = regions(fm.graph)
# Likewise the Euler and convention checks (diagram/WiringCheck.jl) — they don't depend on the cuts.
euler(fm::CircularMorphismGraph) = euler(fm.graph)
check_wiring(fm::CircularMorphismGraph) = check_wiring(fm.graph)

"""
    circular_morphism(m::MorphismGraph) -> CircularMorphismGraph

Like `circular(m::MorphismGraph)`, but `cut1`/`cut2` are preserved — the
counterpart the `circular` docstring above flags as missing.
"""
circular_morphism(m::MorphismGraph) = CircularMorphismGraph(circular(m.graph), m.cut1, m.cut2)

# `_bottom_leaves`/`_top_leaves` (morphism/MorphismGraph.jl) depend only on
# `word`/`cut1`/`cut2`, never on the node type — the same gap-walk logic,
# reimplemented directly on `CircularMorphismGraph` here (instead of constructing a
# `WordGraph` with empty nodes/edges just for the cut computation, which would
# needlessly trigger `_trace_cells`).
function _bottom_leaves(m::CircularMorphismGraph)
    n = length(m.graph.word)
    n == 0 && return Int[]
    m.cut1 == m.cut2 == _EMPTY_BOTTOM_CUT && return Int[]
    m.cut1 == m.cut2 && return collect(1:n)
    ks = Int[]; k = mod1(m.cut1 + 1, n)
    while true
        push!(ks, k)
        k == mod1(m.cut2, n) && break
        k = mod1(k + 1, n)
    end
    return ks
end
function _top_leaves(m::CircularMorphismGraph)
    n = length(m.graph.word)
    n == 0 && return Int[]
    m.cut1 == m.cut2 == _EMPTY_BOTTOM_CUT && return collect(1:n)
    m.cut1 == m.cut2 && return Int[]
    ks = Int[]; k = mod1(m.cut2 + 1, n)
    while true
        push!(ks, k)
        k == mod1(m.cut1, n) && break
        k = mod1(k + 1, n)
    end
    return ks
end

"The bottom (domain) word of a CircularMorphismGraph, as a colour vector."
bottom(m::CircularMorphismGraph) = [leaf_colour(m.graph, k) for k in _bottom_leaves(m)]

"The top (codomain) word of a CircularMorphismGraph, as a colour vector (reversed top arc)."
top(m::CircularMorphismGraph) = reverse([leaf_colour(m.graph, k) for k in _top_leaves(m)])

"""
    _circular_left_mark(m::CircularMorphismGraph) -> MarkSpec
    _circular_right_mark(m::CircularMorphismGraph) -> MarkSpec

Literal analogue of `_left_mark`/`_right_mark` (render/tutte/Display.jl), just
over `CircularMorphismGraph.cut1`/`cut2` instead of `MorphismGraph.cut1`/`cut2`.
"""
function _circular_left_mark(m::CircularMorphismGraph)
    n = length(m.graph.word)
    n == 0 && return :left
    return mod(m.cut1, n)
end
function _circular_right_mark(m::CircularMorphismGraph)
    n = length(m.graph.word)
    n == 0 && return :right
    m.cut1 == m.cut2 && return (:opposite, mod(m.cut1, n))
    return mod(m.cut2, n)
end
