# morphism/MorphismGraph.jl — a WordGraph read as a bottom → top morphism (two cuts).
# The two conversions graph↔path live in GraphToPath.jl / PathToGraph.jl.

# ---- leaf colour (the truth: the edge at a leaf, not the normalized word) --

"Colour of the edge attached to boundary `Leaf(k)` in `g` (0 if none)."
function leaf_colour(g::WordGraph, k::Int)
    for e in g.edges
        e.a isa Leaf && e.a.k == k && return e.colour
        e.b isa Leaf && e.b.k == k && return e.colour
    end
    return 0
end

# ---- MorphismGraph: a diagram read as bottom → top -------------------------

# sentinel cut value (see MorphismGraph docstring below): cut1 == cut2 == this means
# "bottom = ε, top = whole boundary", the OTHER degenerate case from cut1==cut2==0.
# Never a real gap position (those are always 0..n-1 via mod/mod1).
const _EMPTY_BOTTOM_CUT = -1

"""
    MorphismGraph

A `WordGraph` with two cuts making it a morphism. `cut1`/`cut2` are gap positions
(the gap sits between leaf `cut` and leaf `cut+1`, cyclically). The bottom word is
leaves `cut1+1 … cut2` (forward); the top word is leaves `cut2+1 … cut1` REVERSED.
`cut1 == cut2 == 0` gives bottom = whole boundary, top = ε (the "reduce to ε"
picture). `cut1 == cut2 == _EMPTY_BOTTOM_CUT` (`== -1`) is the OTHER degenerate
reading: bottom = ε, top = whole boundary (the "unit" picture), which `flip` and
`double_leaf` need: a raw `cut1 == cut2` alone
cannot distinguish the two (the walk that reads out a nonempty arc always visits
≥ 1 leaf, so it can only ever mean "whole", never "ε", for a NONEMPTY boundary), so
the sentinel is the one extra bit needed to tell them apart. Never produced by
`mod`-normalised cut arithmetic elsewhere (that always lands on plain `0`), so it
only ever appears via `flip`.
"""
struct MorphismGraph
    graph::WordGraph
    cut1::Int
    cut2::Int
end

"Bottom-and-top view of a plain graph: whole boundary as bottom, ε as top."
morphism_graph(g::WordGraph) = MorphismGraph(g, 0, 0)
morphism_graph(g::WordGraph, cut1::Int, cut2::Int) = MorphismGraph(g, cut1, cut2)

# leaf indices of the bottom arc (cut1+1 … cut2, cyclic), forward
function _bottom_leaves(m::MorphismGraph)
    n = length(m.graph.word)
    n == 0 && return Int[]
    m.cut1 == m.cut2 == _EMPTY_BOTTOM_CUT && return Int[]     # bottom = ε (unit picture)
    m.cut1 == m.cut2 && return collect(1:n)          # whole boundary
    ks = Int[]; k = mod1(m.cut1 + 1, n)
    while true
        push!(ks, k)
        k == mod1(m.cut2, n) && break
        k = mod1(k + 1, n)
    end
    return ks
end

# leaf indices of the top arc (cut2+1 … cut1, cyclic), forward (NOT yet reversed)
function _top_leaves(m::MorphismGraph)
    n = length(m.graph.word)
    n == 0 && return Int[]
    m.cut1 == m.cut2 == _EMPTY_BOTTOM_CUT && return collect(1:n)  # top = whole boundary
    m.cut1 == m.cut2 && return Int[]
    ks = Int[]; k = mod1(m.cut2 + 1, n)
    while true
        push!(ks, k)
        k == mod1(m.cut1, n) && break
        k = mod1(k + 1, n)
    end
    return ks
end

"The bottom (domain) word of the morphism, as a colour vector."
bottom(m::MorphismGraph) = [leaf_colour(m.graph, k) for k in _bottom_leaves(m)]

"The top (codomain) word of the morphism, as a colour vector (reversed top arc)."
top(m::MorphismGraph) = reverse([leaf_colour(m.graph, k) for k in _top_leaves(m)])

function Base.show(io::IO, m::MorphismGraph)
    b = bottom(m); t = top(m)
    bs = isempty(b) ? "ε" : join(b, " ")
    ts = isempty(t) ? "ε" : join(t, " ")
    print(io, "morphism  ", bs, "  →  ", ts, "        (bottom → top)")
end

"""
    show_morphism(m::MorphismGraph; moves = false) -> String

Console rendering: `bottom → top`, and with `moves = true` the peeled reduction
moves (from `graph_to_path`) one per line.
"""
function show_morphism(m::MorphismGraph; moves::Bool = false)
    io = IOBuffer()
    print(io, m)
    if moves
        d, status = graph_to_path(m.graph)
        println(io)
        if status === :stalled
            println(io, "  (moves: reduction STALLED — graph not fully peelable here)")
        end
        for (i, v) in enumerate(d.vertices)
            kind = v.kind === VDot ? "dot" : v.kind === VTrivalent ? "trivalent" :
                   v.kind === VComm ? "comm" : "braid"
            win = join(d.layers[i][[mod1(v.pos + j, length(d.layers[i])) for j in 0:(v.pl-1)]], "")
            rep = isempty(v.repl) ? "ε" : join(v.repl, "")
            println(io, "  via: ", kind, "(", win, "→", rep, ")@", v.pos)
        end
    end
    return String(take!(io))
end

# ---- Soergel grading (degree) on diagrams -----------------------------------
#
# The diagrammatic Soergel grading: each dot has degree +1, each trivalent vertex
# has degree −1 (EMTW convention), `:braid` nodes are degree 0 (they are
# isomorphisms/relations, not generators that shift degree). This is checked
# against the COMBINATORIAL grading `defect(word, e)` (LightLeaves.jl) on every
# light leaf — see the `double_leaves` test: the two gradings must
# agree, since they are two ways of computing the same thing (EMTW ch. 5-6).

"Number of `:dot` nodes in `g`."
dot_count(g::WordGraph) = count(nd -> nd.kind === :dot, g.nodes)

"Number of `:trivalent` nodes in `g`."
trivalent_count(g::WordGraph) = count(nd -> nd.kind === :trivalent, g.nodes)

"""
    degree(g::WordGraph) -> Int

The Soergel grading of the diagram `g`: `dot_count(g) − trivalent_count(g)`.
`:braid` nodes contribute 0 (they carry no grading shift)."
"""
degree(g::WordGraph) = dot_count(g) - trivalent_count(g)

"""
    degree(m::MorphismGraph) -> Int

The Soergel grading of the morphism `m`, i.e. of its underlying graph.
"""
degree(m::MorphismGraph) = degree(m.graph)

"""
    degree(a::SoergelPoly, m) -> Int

The grading of a polynomial-scaled morphism `a ⊗ m`: `2*degree(a) + degree(m)`
(each `α_i` has ring-degree 2, matching `degree(::SoergelPoly)`, Ring.jl).
"""
degree(a::SoergelPoly, m) = 2 * degree(a) + degree(m)

