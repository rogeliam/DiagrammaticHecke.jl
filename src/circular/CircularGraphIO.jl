# circular/CircularGraphIO.jl — save/load a CircularGraph/CircularMorphismGraph as versioned, git-diffable text
#
# The `.fwg`/`.fwgm` analogue of diagram/GraphIO.jl (`.wg`) and
# diagram/MorphismIO.jl (`.wgm`), for the PARALLEL diagram type with open arm
# count (circular/CircularGraph.jl). Layout, port tokens, readers/writers and the
# `---` container all come from diagram/IO.jl — only the NODE LINE changes, because
# a `CircularNode` carries an arbitrary number of arms instead of a fixed
# `colours`/`m` pair:
#
#     n <kind> <arm1>,<arm2>,...,<armk>
#
# e.g. `n mono 2,2,2,2` (a merged 4-armed mono node) or `n braid 1,2,1,2,1,2`.
# `kind` is stored ALONGSIDE the arms (not re-derived from them) and set on
# read via `_unchecked_circularnode` — that lets even transient node shapes
# forbidden by `circular_node`'s public constructor (2-arm "beads") round-trip, in
# case they ever show up in a test fixture; see CircularGraph.jl's header on
# `_unchecked_circularnode`.

const FWG_VERSION = 1

# ---- node <-> line (the only piece that is specific to `.fwg`) ----------------

_circularnode_to_line(nd::CircularNode) = "n $(nd.kind) $(join(nd.arms, ","))"

function _line_to_circularnode(line::AbstractString)
    parts = split(line)                        # ["n", kind, "a1,a2,..."]
    kind = Symbol(parts[2])
    arms = [parse(Int, a) for a in split(parts[3], ",")]
    return _unchecked_circularnode(kind, arms)
end

const FWG_FORMAT = GraphTextFormat{CircularGraph,CircularNode}(
    "fwg", [FWG_VERSION], FWG_VERSION, _circularnode_to_line, _line_to_circularnode,
    (word, nodes, edges) -> CircularGraph(word, nodes, edges))

# ---- the graph <-> text --------------------------------------------------------

"""
    circulargraph_to_string(g::CircularGraph) -> String

Serialise `g` to the versioned `.fwg` text format (see diagram/IO.jl). Inverse
of `circulargraph_from_string`; round-trips exactly (nodes via `kind`+`arms`, so even
transient/`_unchecked_circularnode`-only shapes survive).
"""
circulargraph_to_string(g::CircularGraph) = graph_to_text(FWG_FORMAT, g)

"""
    _circulargraph_from_string_with_cuts(s::AbstractString) -> (CircularGraph, Union{Nothing,Tuple{Int,Int}})

Parsing core shared by `circulargraph_from_string` (drops the cuts) and
`circularmorphismgraph_from_string` (requires them).
"""
_circulargraph_from_string_with_cuts(s::AbstractString) = graph_from_text(FWG_FORMAT, s)

"""
    circulargraph_from_string(s::AbstractString) -> CircularGraph

Parse a `.fwg` text back into a CircularGraph. Inverse of `circulargraph_to_string`. A
`cuts` line, if present, is IGNORED here — use `circularmorphismgraph_from_string` to
recover it.
"""
circulargraph_from_string(s::AbstractString) = first(_circulargraph_from_string_with_cuts(s))

"Write `g` to `path` in the `.fwg` text format."
save_circulargraph(path::AbstractString, g::CircularGraph) =
    open(io -> write(io, circulargraph_to_string(g)), path, "w")

"Read a CircularGraph from a `.fwg` file written by `save_circulargraph`."
load_circulargraph(path::AbstractString) = circulargraph_from_string(read(path, String))

# ---- a single CircularMorphismGraph <-> text ---------------------------------------

"""
    circularmorphismgraph_to_string(m::CircularMorphismGraph) -> String

Serialise `m` to `.fwg` text plus a trailing `cuts <cut1> <cut2>` line. Inverse of
`circularmorphismgraph_from_string`. Analogue of `morphismgraph_to_string`
(diagram/MorphismIO.jl).
"""
circularmorphismgraph_to_string(m::CircularMorphismGraph) = morphism_to_text(FWG_FORMAT, m)

"""
    circularmorphismgraph_from_string(s::AbstractString) -> CircularMorphismGraph

Parse text written by `circularmorphismgraph_to_string` back into a `CircularMorphismGraph`.
A `cuts` line is REQUIRED, exactly as in `morphismgraph_from_string`.
"""
circularmorphismgraph_from_string(s::AbstractString) =
    morphism_from_text(FWG_FORMAT, s, CircularMorphismGraph,
                       "circularmorphismgraph_from_string")

"Write `m` to `path` in the `.fwg` (+ `cuts`) text format."
save_circular_morphism(path::AbstractString, m::CircularMorphismGraph) =
    open(io -> write(io, circularmorphismgraph_to_string(m)), path, "w")

"Read a CircularMorphismGraph from a file written by `save_circular_morphism`."
load_circular_morphism(path::AbstractString) = circularmorphismgraph_from_string(read(path, String))

# ---- a container of many CircularMorphismGraphs (.fwgm) -----------------------------

"""
    CircularMorphismEntry

One record in a `.fwgm` container: a free-form metadata comment line (preserved
VERBATIM, not parsed) alongside the `CircularMorphismGraph` it describes. Analogue of
`MorphismEntry` (diagram/MorphismIO.jl); a typical use is one entry per
`Zamo(i,j)` term (`# zamo=1,8` etc.).
"""
struct CircularMorphismEntry
    meta::String
    morphism::CircularMorphismGraph
end

"""
    save_circular_morphisms(path::AbstractString, entries::Vector{CircularMorphismEntry})

Write several `CircularMorphismGraph`s to one `.fwgm` file, each preceded by its
metadata comment line and separated by a bare `---` line. Order is preserved.
"""
save_circular_morphisms(path::AbstractString, entries::Vector{CircularMorphismEntry}) =
    save_entry_file(path, entries, circularmorphismgraph_to_string)

"""
    load_circular_morphisms(path::AbstractString) -> Vector{CircularMorphismEntry}

Read a `.fwgm` file written by `save_circular_morphisms` back into its entries, in order.
"""
load_circular_morphisms(path::AbstractString) =
    load_entry_file(path, CircularMorphismEntry, circularmorphismgraph_from_string)

# ---- batch helper: save every Zamo(i,j) term as one .fwgm ----------------------

"""
    save_zamo_circular_morphisms(path::AbstractString; i_range = 1:14, j_range = 1:14)

Write `circular_morphism(Zamo(i,j))` for every `(i,j)` pair with `i != j` (`Zamo(i,i)`
is the whole 14-step cycle, not a "term") to one `.fwgm` batch file — the circular
analogue of what `save_morphisms` would do for a `Vector{MorphismEntry}` built
from `Zamo`. Metadata line records the pair as `zamo=i,j`.
"""
function save_zamo_circular_morphisms(path::AbstractString; i_range = 1:14, j_range = 1:14)
    entries = CircularMorphismEntry[]
    for i in i_range, j in j_range
        i == j && continue
        push!(entries, CircularMorphismEntry("zamo=$i,$j", circular_morphism(Zamo(i, j))))
    end
    save_circular_morphisms(path, entries)
    return entries
end
