# reduce/PathDiagram.jl  —  the guided position-faithful reduction path/vertex
# model, for drawing and `reduction_path_string`/`show_reduction_path`.
#
# A diagram here is the LAYER/VERTEX SEQUENCE
# produced by the guided reducer (`_reduction_layers` in render/Morphism.jl): a
# boundary word at the bottom, then one vertex per layer taking each word to the
# next.

# ---- the vertex model -------------------------------------------------------
#
# A vertex is what one guided move does between two consecutive layers. Its kind is
# read off (window length pl, replacement length rl):
#   :dot        a dot generator      (pl,rl) with |pl-rl| via a dot rule (ii↔i / ii↔ε
#               / ε↔ii): degree-1 touch on one strand.
#   :trivalent  ii → i or i → ii      (merge/split, one colour)
#   :braid      121↔212 / 232↔323 / 1313-type commutation 13↔31 (length-neutral)
# We store the raw window data so a rule can inspect/modify it and rebuild the path.

@enum VKind VDot VTrivalent VBraid VComm

"""
    Vertex

One generator on the reduction path: its `kind`, the colours involved, and the raw
window it acts on (start `pos`, source length `pl`, replacement letters `repl`) on
the layer BELOW it. `ridx` = the indices, in the layer ABOVE, of the letters this
vertex produced (the replacement).
"""
struct Vertex
    kind::VKind
    colours::Vector{Int}
    pos::Int
    pl::Int
    repl::Vector{Int}
    ridx::Vector{Int}
end

# classify a (window, replacement) into a vertex kind + colours
function _classify(window::Vector{Int}, repl::Vector{Int})
    us = sort(unique(vcat(window, repl)))
    if length(us) == 1                                 # single colour: dot or trivalent
        # ii→i / i→ii  = trivalent;  ii→ε / ε→ii / i→ (nothing) = dot-ish
        (length(window) == 2 && length(repl) == 1) && return (VTrivalent, us)
        (length(window) == 1 && length(repl) == 2) && return (VTrivalent, us)
        return (VDot, us)
    else                                               # two colours: braid or comm
        length(window) == 2 && return (VComm, us)      # 13↔31
        return (VBraid, us)                            # 121↔212 etc.
    end
end

"""
    PathDiagram

A diagram as the guided-reduction layer/vertex sequence of a circular word:
`layers[1]` = boundary word (bottom), `vertices[i]` takes `layers[i]` to
`layers[i+1]`, up to ε (or a stall). Built by `path_diagram(w)`.
"""
struct PathDiagram
    layers::Vector{Vector{Int}}
    vertices::Vector{Vertex}
end

"""
    path_diagram(w::CircularWord) -> PathDiagram

Build the layer/vertex sequence for `w` via the guided position-faithful reducer
(the same one `radial_svg` draws). Each vertex is classified into dot / trivalent /
braid / commutation. Ends at ε if `w` reduces; otherwise stalls at the last layer.
"""
function path_diagram(w::CircularWord)
    rings, steps = _reduction_layers(w)
    verts = Vertex[]
    # `layers` always has one more entry than `vertices`: the word ABOVE each
    # vertex. Compute each result word from the window so the top layer (ε on a
    # full reduction) is present too.
    layers = [copy(r) for r in rings]
    for (i, (wins, _, repl_idx)) in enumerate(steps)
        s, pl, _, repl = wins[1]
        window = [rings[i][mod1(s + j, length(rings[i]))] for j in 0:(pl-1)]
        kind, cols = _classify(window, collect(repl))
        push!(verts, Vertex(kind, cols, s, pl, collect(repl), copy(repl_idx[1])))
    end
    # ensure layers[i+1] exists for the last vertex: when the reduction reached ε,
    # `rings` has exactly `length(steps)` entries, so append the empty top layer.
    if length(layers) == length(verts)
        push!(layers, Int[])
    end
    return PathDiagram(layers, verts)
end

Base.length(d::PathDiagram) = length(d.vertices)

# ---- console: a reduction path as a stack of words --------------------------
#
# Print a morphism/diagram/word as its REDUCTION PATH: one circular word per line,
# from the boundary (top) down to ε, e.g.
#     1221
#     1122
#     22
#     ε
# This is the "path = sequence of words" reading (paper1/01 §3): each line is a
# circular word, consecutive lines differ by one elementary move.

"Render `layers` (a vector of words) as a reduction path, one word per line."
function _path_lines(layers::Vector{Vector{Int}})
    io = IOBuffer()
    for (i, L) in enumerate(layers)
        i == 1 || println(io)
        print(io, isempty(L) ? "ε" : join(L, ""))
    end
    return String(take!(io))
end

"""
    reduction_path_string(x) -> String
    show_reduction_path(x)

Show `x` as a REDUCTION PATH — a column of circular words from the boundary down to
`ε`, one elementary move between consecutive lines. `x` may be a `CircularWord` (its
guided boundary reduction, via `path_diagram`), a `PathDiagram`, or a `WordGraph`
(peeled outside-in via `graph_to_path`; prints a note if the peel stalls). Example:
`show_reduction_path(CircularWord([1,2,2,1]))` prints `1122 / 22 / ε`.
"""
reduction_path_string(d::PathDiagram) = _path_lines(d.layers)
reduction_path_string(w::CircularWord) = _path_lines(path_diagram(w).layers)
function reduction_path_string(g::WordGraph)
    d, status = graph_to_path(g)
    body = _path_lines(d.layers)
    return status === :stalled ?
        body * "\n(⚠ peel STALLED — not fully reducible from the graph this way)" : body
end

show_reduction_path(x) = (println(reduction_path_string(x)); nothing)
