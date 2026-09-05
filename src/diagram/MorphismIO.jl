# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/MorphismIO.jl  —  save/load a MorphismGraph as versioned, git-diffable text
#
# A `MorphismGraph` is just a `WordGraph` (see morphism/MorphismGraph.jl) plus two
# cut positions, so the file format is the SAME `.wg` text with one extra line
# appended:
#
#     #wg 2
#     word 1212121
#     # nodes
#     ...
#     # edges
#     ...
#     cuts 0 3                  the morphism's cut1/cut2 (may be NEGATIVE — see
#                               `_EMPTY_BOTTOM_CUT`, morphism/MorphismGraph.jl)
#
# A single morphism uses `save_morphism`/`load_morphism` (`morphismgraph_to_string`/
# `_from_string` for the in-memory text). A `.wgm` ("word graph, many") file holds
# SEVERAL morphisms in the `---`-separated container of diagram/IO.jl, one free-form
# metadata comment line each (e.g. `# e=1,1,0 f=1,0,1 z=1 degree=2` — preserved
# VERBATIM, not parsed for meaning):
#
#     # e=1,1,0 f=1,0,1 z=1 degree=2
#     #wg 2
#     ...
#     cuts 0 3
#     ---
#     # e=0,1,1 f=1,1,0 z=1 degree=1
#     #wg 2
#     ...
#
# NOTE: `WordGraph`'s `cells` field is ALWAYS recomputed by the 3-argument
# constructor from `(word, nodes, edges)` (diagram/Faces.jl) — it is never written
# here and never read back. `graph_from_text` (diagram/IO.jl) already builds the
# graph that way, so this file only handles the `cuts` line.

# ---- a single morphism <-> text --------------------------------------------

"""
    morphismgraph_to_string(m::MorphismGraph) -> String

Serialise `m` to `.wg` text (see diagram/IO.jl) plus a trailing
`cuts <cut1> <cut2>` line. Inverse of `morphismgraph_from_string`.
"""
morphismgraph_to_string(m::MorphismGraph) = morphism_to_text(WG_FORMAT, m)

"""
    morphismgraph_from_string(s::AbstractString) -> MorphismGraph

Parse text written by `morphismgraph_to_string` back into a `MorphismGraph`. A
`cuts` line is REQUIRED — a string with none (e.g. a bare `.wg` file with no
morphism structure) is an ERROR, not silently read as `(0, 0)`, since that would
silently misrepresent "no cuts recorded" as the very meaningful "whole boundary →
ε" cut.
"""
morphismgraph_from_string(s::AbstractString) =
    morphism_from_text(WG_FORMAT, s, MorphismGraph, "morphismgraph_from_string")

"Write `m` to `path` in the `.wg` (+ `cuts`) text format."
save_morphism(path::AbstractString, m::MorphismGraph) =
    open(io -> write(io, morphismgraph_to_string(m)), path, "w")

"Read a MorphismGraph from a file written by `save_morphism`."
load_morphism(path::AbstractString) = morphismgraph_from_string(read(path, String))

# ---- a container of many morphisms (.wgm) ----------------------------------

"""
    MorphismEntry

One record in a `.wgm` container: a free-form metadata comment line (preserved
VERBATIM, not parsed) alongside the `MorphismGraph` it describes. The metadata is
typically something like `# e=1,1,0 f=1,0,1 z=1 degree=2` (a `double_leaf` recipe,
morphism/LightLeaves.jl) but the format does not care what it says.
"""
struct MorphismEntry
    meta::String
    morphism::MorphismGraph
end

"""
    save_morphisms(path::AbstractString, entries::Vector{MorphismEntry})

Write several morphisms to one `.wgm` file, each preceded by its metadata comment
line and separated by a bare `---` line. Order is preserved.
"""
save_morphisms(path::AbstractString, entries::Vector{MorphismEntry}) =
    save_entry_file(path, entries, morphismgraph_to_string)

"""
    load_morphisms(path::AbstractString) -> Vector{MorphismEntry}

Read a `.wgm` file written by `save_morphisms` back into its entries, in order.
"""
load_morphisms(path::AbstractString) =
    load_entry_file(path, MorphismEntry, morphismgraph_from_string)
