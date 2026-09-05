# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/Graph.jl  —  the WordGraph type
#
# A circular word together with a planar diagram inside the disc it bounds, read
# as a (local picture of a) morphism in the diagrammatic Hecke category of A₃.
# REPRESENTATION. The truth of the diagram is:
#   - `word`  :: the boundary CircularWord i₁…iₙ  (leaf k has colour iₖ);
#   - `nodes` :: the internal vertices, each a `Node` (kind + colour(s));
#   - `edges` :: coloured edges, each connecting two ENDPOINTS.
# An endpoint (`Port`) is either a boundary leaf `Leaf(k)` (position k in the word)
# or a node connector `NodePort(node_index, slot)`. Every edge has a colour and its
# two ports must carry that colour. We do NOT yet store a full rotation system
# (cyclic order at each node) — for the generators and a single join the planarity
# is unambiguous; we add the rotation system when a move actually needs it.
# Diagrams may be DISCONNECTED.

# ---- ports (edge endpoints) -------------------------------------------------

abstract type Port end

"A boundary leaf: position `k` (1-based) in the circular word."
struct Leaf <: Port
    k::Int
end

"A connector on internal node `node` (1-based index into `nodes`), slot `slot`."
struct NodePort <: Port
    node::Int
    slot::Int
end

"""
A closed-loop marker: a free monochrome CIRCLE inside the disc, touching no leaf and
no node. The model has no bare port-less edge, so a free circle of colour `c` is the
edge `Edge(c, Circle(c), Circle(c))` — both endpoints are this marker. Such a circle
evaluates to 0 (the needle/circle relation, R7), so `reduce_graph` sends any diagram
carrying one to the empty combination (K1).
"""
struct Circle <: Port
    colour::Int
end

# ---- internal nodes ---------------------------------------------------------
#
# kind ∈ (:dot, :trivalent, :braid). `colours` is:
#   :dot        one colour  [c]                 (degree 1: 1 slot)
#   :trivalent  one colour  [c]                 (degree 3: slots 1,2,3)
#   :braid      two colours [s,t], with m       (degree 2m: slots 1..2m,
#                                                colour of slot j is s,t,s,t,…)
#
# PLANAR SLOT CONVENTION. A braid node carries no free planar
# orientation beyond ONE binary choice (which parity is which colour). We PIN it
# down by convention rather than a stored field: **the slot indices 1,2,…,2m ARE
# the planar cyclic order around the node** (say, clockwise). Consequences the
# wiring / renderer / rules may rely on:
#   * neighbours in the cyclic order are slots j and j+1 (mod 2m);
#   * slot j sits DIRECTLY OPPOSITE slot j+m (mod 2m) — see `opposite_slot`.
# This is what lets K1 / R5-local decide "what faces what". The parity→colour
# choice (odd=colours[1]) is the remaining binary freedom, already fixed above.
# `is_planar_braid` checks a WordGraph respects this.
struct Node
    kind::Symbol
    colours::Vector{Int}
    m::Int                      # only meaningful for :braid (m_{st}); else 0
end

Base.:(==)(a::Node, b::Node) = a.kind == b.kind && a.colours == b.colours && a.m == b.m
Base.hash(nd::Node, h::UInt) = hash((nd.kind, nd.colours, nd.m), h)

_node_degree(nd::Node) =
    nd.kind === :dot        ? 1 :
    nd.kind === :trivalent  ? 3 :
    nd.kind === :braid      ? 2 * nd.m :
    error("unknown node kind $(nd.kind)")

"""
    opposite_slot(nd::Node, slot::Int) -> Int

The slot directly across the braid node from `slot`, under the planar slot
convention (slots `1..2m` are the cyclic order): `slot ↦ slot + m (mod 2m)`.
Its colour is the OTHER colour (since `m` shifts parity iff `m` is odd — for the
even-shift case both share a colour, which is exactly the `13↔31` commutation
`m = 2`). Only meaningful for `:braid`.
"""
function opposite_slot(nd::Node, slot::Int)
    nd.kind === :braid || error("opposite_slot only defined for :braid nodes")
    twom = 2 * nd.m
    return mod1(slot + nd.m, twom)
end

"Colour carried by slot `slot` of node `nd`."
function _slot_colour(nd::Node, slot::Int)
    if nd.kind === :braid
        # alternating s,t,s,t,… starting with colours[1]
        return isodd(slot) ? nd.colours[1] : nd.colours[2]
    else
        return nd.colours[1]
    end
end

# ---- an edge (coloured, two ports) ------------------------------------------

struct Edge
    colour::Int
    a::Port
    b::Port
end

Base.:(==)(x::Edge, y::Edge) = x.colour == y.colour && x.a == y.a && x.b == y.b

# ---- the diagram ------------------------------------------------------------

"""
    Cells

The cells (regions) of a diagram, computed outside-in by face tracing (see
diagram/Faces.jl) and **always carried along**. A pure function of
`(word, nodes, edges)`, hence irrelevant to the graph's `==`/`hash`.

- `ncells`           number of inner cells (each free circle counts as one).
- `gap_cell[k]`      cell id of the gap between leaf `k` and `k+1` (boundary arc).
- `sector_cell[i]`   per node `i`: cell id of the sector clockwise AFTER slot `s`
  (length = degree of the node; 0 = slot not wired).
- `port_cell[p]`     `(left, right)` cell id of the edge at port `p`
  (0 = outer face).
"""
struct Cells
    ncells::Int
    gap_cell::Vector{Int}
    sector_cell::Vector{Vector{Int}}
    port_cell::Dict{Port, NTuple{2, Int}}
end

"""
    WordGraph

A planar word diagram: a boundary `word :: CircularWord`, internal `nodes`,
and coloured `edges` between ports (leaves or node connectors). Internals WILL
change.

The `cells` field is filled automatically by the 3-argument constructor
(diagram/Faces.jl), so `WordGraph(word, nodes, edges)` keeps working everywhere.
"""
struct WordGraph
    word::CircularWord
    nodes::Vector{Node}
    edges::Vector{Edge}
    cells::Cells
end

boundary(g::WordGraph) = g.word
Base.length(g::WordGraph) = length(g.word)      # number of boundary leaves

"""
    is_planar_braid(g::WordGraph) -> Bool

Sanity-check that every `:braid` node in `g` has exactly its `2m` slots wired
(each slot `1..2m` used exactly once across the edges), so the planar slot
convention above is meaningful. Does NOT check global planarity of the drawing —
only that the slot indices are the complete, non-repeated `1..2m`.
"""
function is_planar_braid(g::WordGraph)
    for (ni, nd) in enumerate(g.nodes)
        nd.kind === :braid || continue
        used = Int[]
        for e in g.edges
            for p in (e.a, e.b)
                p isa NodePort && p.node == ni && push!(used, p.slot)
            end
        end
        sort(used) == collect(1:(2 * nd.m)) || return false
    end
    return true
end

# ---------------------------------------------------------------------------
# Generators — each is a complete valid diagram on a small circular word.
# ---------------------------------------------------------------------------

"""
    dot(i) -> WordGraph

The dot generator on the 1-letter word `(i)`: a single leaf of colour `i` joined
to a degree-1 dot node.
"""
function dot(i::Integer)
    word = CircularWord([i])
    nodes = [Node(:dot, [Int(i)], 0)]
    edges = [Edge(Int(i), Leaf(1), NodePort(1, 1))]
    return WordGraph(word, nodes, edges)
end

"""
    trivalent(i) -> WordGraph

The trivalent generator on `(i i i)`: three leaves of colour `i` joined to the
three slots of one trivalent node.
"""
function trivalent(i::Integer)
    c = Int(i)
    word = CircularWord([c, c, c])
    nodes = [Node(:trivalent, [c], 0)]
    # ORIENTATION (diagram/Faces.jl): boundary ccw, arms cw — a node at several
    # leaves runs the leaf ring BACKWARDS. Slot 1 holds leaf 1 fixed.
    edges = [Edge(c, Leaf(mod1(2 - s, 3)), NodePort(1, s)) for s in 1:3]
    return WordGraph(word, nodes, edges)
end

"""
    braid(s, t; m = 3) -> WordGraph

The 2m-valent braid generator for the two distinct colours `s ≠ t`: the boundary
word is the alternating `(s t s t …)` of length `2m`, each leaf wired to the
matching slot of one braid node. In A₃: `braid(1,3; m=2)` (word 1313),
`braid(1,2)` / `braid(2,3)` (m=3, words 121212 / 232323).
"""
function braid(s::Integer, t::Integer; m::Integer = 3)
    s == t && error("braid needs two distinct colours")
    s, t, m = Int(s), Int(t), Int(m)
    letters = [isodd(k) ? s : t for k in 1:(2m)]
    word = CircularWord(letters)
    nodes = [Node(:braid, [s, t], m)]
    # arms run cw against the ccw leaf ring (diagram/Faces.jl): slot s holds leaf
    # mod1(2 - s, 2m). Colours stay consistent because 2m is even, so the
    # reflection preserves parity.
    edges = [Edge(letters[mod1(2 - s, 2m)], Leaf(mod1(2 - s, 2m)), NodePort(1, s))
             for s in 1:(2m)]
    return WordGraph(word, nodes, edges)
end

# ---------------------------------------------------------------------------
# A restricted JOIN: connect two same-colour trivalents.
# ---------------------------------------------------------------------------

"""
    join_trivalents(i) -> WordGraph

Two trivalent nodes of the same colour `i`, joined by one internal edge (slot 3 of
each), leaving 4 boundary leaves of colour `i`. Boundary word `(i i i i)`.

Picture (colour i):   leaf1 leaf2          leaf3 leaf4
                          \\  /                \\  /
                        [node1]==(internal)==[node2]
This is the associativity/"H" shape — two merges glued along one strand. A first,
hand-wired example of a diagram with an internal edge and two components' worth of
structure joined into one.
"""
function join_trivalents(i::Integer)
    c = Int(i)
    word = CircularWord([c, c, c, c])
    nodes = [Node(:trivalent, [c], 0), Node(:trivalent, [c], 0)]
    edges = [
        # arms cw, leaf ring ccw (diagram/Faces.jl): at each node the leaves run
        # backwards; the internal edge closes the ring.
        Edge(c, Leaf(2), NodePort(1, 1)),      # node1: slots 1,2 -> leaves 2,1
        Edge(c, Leaf(1), NodePort(1, 2)),
        Edge(c, Leaf(4), NodePort(2, 1)),      # node2: slots 1,2 -> leaves 4,3
        Edge(c, Leaf(3), NodePort(2, 2)),
        Edge(c, NodePort(1, 3), NodePort(2, 3)), # the internal join edge
    ]
    return WordGraph(word, nodes, edges)
end

# ---------------------------------------------------------------------------
# Console (Unicode) rendering — first pass
# ---------------------------------------------------------------------------
#
# The boundary is cyclic; for a readable first pass we "cut the circle open" and
# list the boundary word on top, then describe the nodes and edges below in a
# compact, colour-tagged form. A true planar picture comes later — this pass is to
# SEE the structure (which leaves/nodes each edge connects) at a glance.

# per-colour ANSI 256 codes
const _WG_COLOR = Dict(1 => 33, 2 => 203, 3 => 28)   # blue, tomato, spinach

_wg_paint(c::Int, s; on = true) = on ? "\e[38;5;$(get(_WG_COLOR, c, 7))m$s\e[0m" : s

_port_str(p::Leaf) = "leaf$(p.k)"
_port_str(p::NodePort) = "n$(p.node).$(p.slot)"
_port_str(p::Circle) = "○$(p.colour)"          # free-circle marker

_node_str(idx, nd::Node) =
    nd.kind === :dot        ? "n$idx = dot($(nd.colours[1]))" :
    nd.kind === :trivalent  ? "n$idx = trivalent($(nd.colours[1]))" :
    "n$idx = braid($(nd.colours[1]),$(nd.colours[2]); m=$(nd.m), $(2nd.m)-valent)"

"""
    show_wordgraph(g; color = true)

A first console rendering: the boundary word (leaf index + colour), the internal
nodes, and the coloured edges (which two ports each connects). Not yet a planar
picture — a readable structural dump.
"""
function show_wordgraph(g::WordGraph; color::Bool = true)
    io = IOBuffer()
    lets = letters(g.word)
    println(io, "WordGraph  (boundary length ", length(lets), ")")
    # boundary row
    print(io, "  boundary: ")
    if isempty(lets)
        print(io, "ε")
    else
        for (k, c) in enumerate(lets)
            print(io, _wg_paint(c, string(c); on = color))
            k < length(lets) && print(io, " ")
        end
    end
    println(io)
    print(io, "            ")
    for k in 1:length(lets)
        print(io, "^", k < length(lets) ? " " : "")
    end
    println(io, "   (leaf indices)")
    # nodes
    println(io, "  nodes:")
    isempty(g.nodes) && println(io, "    (none)")
    for (idx, nd) in enumerate(g.nodes)
        cstr = nd.kind === :braid ?
            _wg_paint(nd.colours[1], string(nd.colours[1]); on = color) * "/" *
            _wg_paint(nd.colours[2], string(nd.colours[2]); on = color) :
            _wg_paint(nd.colours[1], string(nd.colours[1]); on = color)
        println(io, "    ", _node_str(idx, nd), "   [deg ", _node_degree(nd),
                ", colour ", cstr, "]")
    end
    # edges
    println(io, "  edges:")
    isempty(g.edges) && println(io, "    (none)")
    for e in g.edges
        line = "    " * _port_str(e.a) * " ── " * _port_str(e.b)
        println(io, _wg_paint(e.colour, line; on = color), "   (colour ",
                _wg_paint(e.colour, string(e.colour); on = color), ")")
    end
    return String(take!(io))
end

Base.show(io::IO, g::WordGraph) =
    print(io, "WordGraph(", isempty(letters(g.word)) ? "ε" : join(letters(g.word)),
          ", ", length(g.nodes), " nodes, ", length(g.edges), " edges)")
