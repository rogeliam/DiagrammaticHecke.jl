# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/GraphIO.jl  —  save/load a WordGraph as versioned, git-diffable text
#
# Same spirit as Checkpoint.jl (the .cwg circular-word-graph format): a plain-text,
# self-describing, one-record-per-line file so a diagram can be stored as test data,
# diffed in git, and eyeballed. A WordGraph is small (a boundary word + a handful of
# nodes and edges), so a few dependency-free lines suffice.
#
# The layout, the port tokens and the readers/writers live in diagram/IO.jl — this
# file only supplies the `.wg` NODE LINE and the public names. See that header for
# the file layout.
#
#     n dot 2                     :dot,       colour 2
#     n trivalent 1               :trivalent, colour 1
#     n braid 1,2 3               :braid,     colours 1,2  m=3
#
# VERSIONING. The writer emits `#wg 2`, the reader accepts `v in (1, 2)`. The
# version marks whether the format can carry a `cuts` line and the `C` port
# token: a v2 file round-trips both, a v1 file has neither and so loses nothing
# by being read as one.

const WG_VERSION = 2

# ---- node <-> line (the only piece that is specific to `.wg`) -----------------

function _node_to_line(nd::Node)
    cols = join(nd.colours, ",")
    nd.kind === :braid ? "n braid $cols $(nd.m)" : "n $(nd.kind) $cols"
end

function _line_to_node(line::AbstractString)
    parts = split(line)                       # ["n", kind, cols, (m)]
    kind = Symbol(parts[2])
    cols = [parse(Int, c) for c in split(parts[3], ",")]
    m = kind === :braid ? parse(Int, parts[4]) : 0
    return Node(kind, cols, m)
end

const WG_FORMAT = GraphTextFormat{WordGraph,Node}(
    "wg", [1, 2], WG_VERSION, _node_to_line, _line_to_node,
    (word, nodes, edges) -> WordGraph(word, nodes, edges))

# ---- the graph <-> text ------------------------------------------------------

"""
    wordgraph_to_string(g::WordGraph) -> String

Serialise `g` to the versioned `.wg` text format (see diagram/IO.jl). The inverse
of `wordgraph_from_string`; round-trips exactly.
"""
wordgraph_to_string(g::WordGraph) = graph_to_text(WG_FORMAT, g)

"""
    _wordgraph_from_string_with_cuts(s::AbstractString) -> (WordGraph, Union{Nothing,Tuple{Int,Int}})

The parsing core shared by `wordgraph_from_string` (drops the cuts) and
`morphismgraph_from_string` in diagram/MorphismIO.jl (requires them).
"""
_wordgraph_from_string_with_cuts(s::AbstractString) = graph_from_text(WG_FORMAT, s)

"""
    wordgraph_from_string(s::AbstractString) -> WordGraph

Parse a `.wg` text back into a WordGraph. Inverse of `wordgraph_to_string`. A
`cuts` line, if present (e.g. because `s` was written by
`morphismgraph_to_string`), is simply IGNORED here — use
`morphismgraph_from_string` to recover it.
"""
wordgraph_from_string(s::AbstractString) = first(_wordgraph_from_string_with_cuts(s))

"Write `g` to `path` in the `.wg` text format."
save_wordgraph(path::AbstractString, g::WordGraph) =
    open(io -> write(io, wordgraph_to_string(g)), path, "w")

"Read a WordGraph from a `.wg` file written by `save_wordgraph`."
load_wordgraph(path::AbstractString) = wordgraph_from_string(read(path, String))
