# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/IO.jl  —  the line-based text-format core shared by the graph formats
#
# There are two graph record formats in the package, and they are the SAME format
# with one line swapped:
#
#     .wg  / .wgm    WordGraph     / MorphismGraph          (diagram/GraphIO.jl,
#                                                            diagram/MorphismIO.jl)
#     .fwg / .fwgm   CircularGraph / CircularMorphismGraph  (circular/CircularGraphIO.jl)
#
# Everything except the NODE line is identical: the version tag, the `word` line,
# the `# nodes` / `# edges` sections, the port tokens, the optional `cuts` line,
# and the `---`-separated multi-record container with one free-form metadata
# comment per record. Only a node line differs, because a `Node` carries
# `colours` + `m` while a `CircularNode` carries an arbitrary `arms` vector.
#
# So a format is described by a `GraphTextFormat` value — tag, accepted read
# versions, the version to write, the node writer/parser, and the graph
# constructor — and the readers/writers below take that value. Adding a third
# graph type means writing two node functions, not another 110 lines.
#
# The `.cd`/`.cdb` formats (circular/CircularDecoratedIO.jl) are NOT built on this:
# they do not carry a graph record of their own, they WRAP a `.fwg` body with a
# label header, and so already reuse `circulargraph_to_string`.
#
# FILE LAYOUT (identical for both, `<tag>` being `wg` or `fwg`):
#
#     #<tag> <version>          format version
#     word 1212121              boundary word (compact; "eps" for ε)
#     # nodes                   one per line, in index order (1-based)
#     n <kind> ...                format-specific, see the node functions
#     # edges                   one per line: colour  portA  portB
#     e 1 L1 n2.3                 colour 1, leaf 1 — node 2 slot 3
#     cuts 0 3                  OPTIONAL, morphism cuts
#
# A port is either `L<k>` (boundary leaf k), `n<node>.<slot>` (node connector), or
# `C<colour>` (the free-circle marker, diagram/Graph.jl's `Circle`).

# ---- word <-> token ----------------------------------------------------------

_wg_word_token(w::CircularWord) = isempty(w) ? "eps" : w.str
_wg_token_word(tok::AbstractString) =
    tok == "eps" ? EMPTY : CircularWord([Int(c - '0') for c in tok])

# ---- port <-> token ----------------------------------------------------------

_port_to_token(p::Leaf)     = "L$(p.k)"
_port_to_token(p::NodePort) = "n$(p.node).$(p.slot)"
_port_to_token(p::Circle)   = "C$(p.colour)"     # free-circle marker (Graph.jl)

function _token_to_port(tok::AbstractString)
    if startswith(tok, "L")
        return Leaf(parse(Int, tok[2:end]))
    elseif startswith(tok, "n")
        node, slot = split(tok[2:end], ".")
        return NodePort(parse(Int, node), parse(Int, slot))
    elseif startswith(tok, "C")
        return Circle(parse(Int, tok[2:end]))
    end
    error("bad port token: $tok")
end

# ---- a format description ----------------------------------------------------

"""
    GraphTextFormat{G,N}

Describes one of the line-based graph record formats. `G` is the graph type read
back, `N` the node type. Fields:

  * `tag` — the magic word after `#` on the first line (`"wg"`, `"fwg"`), also
    used in error messages as the file extension;
  * `read_versions` — versions the reader accepts;
  * `write_version` — the version the writer stamps;
  * `node_to_line` — `N -> String`, the node line WITHOUT the trailing newline;
  * `line_to_node` — `String -> N`, its inverse;
  * `build` — `(word, Vector{N}, Vector{Edge}) -> G`.

The node functions are the only thing that differs between `.wg` and `.fwg`.
"""
struct GraphTextFormat{G,N}
    tag::String
    read_versions::Vector{Int}
    write_version::Int
    node_to_line::Any
    line_to_node::Any
    build::Any
end

_fmt_versions_text(fmt::GraphTextFormat) =
    length(fmt.read_versions) == 1 ? string(fmt.read_versions[1]) :
    join(fmt.read_versions, ", ", " or ")

# ---- graph <-> text ----------------------------------------------------------

"""
    graph_to_text(fmt::GraphTextFormat, g) -> String

Serialise `g` (anything with `word`/`nodes`/`edges`) in `fmt`. Inverse of
[`graph_from_text`](@ref) up to the optional `cuts` line, which a bare graph does
not carry.
"""
function graph_to_text(fmt::GraphTextFormat, g)
    io = IOBuffer()
    println(io, "#", fmt.tag, " ", fmt.write_version)
    println(io, "word ", _wg_word_token(g.word))
    println(io, "# nodes")
    for nd in g.nodes
        println(io, fmt.node_to_line(nd))
    end
    println(io, "# edges")
    for e in g.edges
        println(io, "e ", e.colour, " ", _port_to_token(e.a), " ", _port_to_token(e.b))
    end
    return String(take!(io))
end

"""
    graph_from_text(fmt::GraphTextFormat{G,N}, s) -> (G, Union{Nothing,Tuple{Int,Int}})

Parse text in `fmt`. Returns the graph and the `cuts` pair if the text carried a
`cuts` line, `nothing` otherwise — the graph readers drop it, the morphism
readers require it. `cuts` is parsed wherever it appears and may hold NEGATIVE
numbers (the `_EMPTY_BOTTOM_CUT` sentinel `-1`, morphism/MorphismGraph.jl).

Unknown lines are skipped rather than rejected, which is what lets `.cd`
(circular/CircularDecoratedIO.jl) hand its whole text down after stripping only
the lines it owns.
"""
function graph_from_text(fmt::GraphTextFormat{G,N}, s::AbstractString) where {G,N}
    word = EMPTY
    nodes = N[]
    edges = Edge[]
    cuts = nothing
    section = :header
    tagline = "#" * fmt.tag
    for raw in split(s, '\n')
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, tagline)
            v = parse(Int, split(line)[2])
            v in fmt.read_versions ||
                error("unsupported .$(fmt.tag) version $v (expected $(_fmt_versions_text(fmt)))")
        elseif startswith(line, "word ")
            word = _wg_token_word(split(line)[2])
        elseif line == "# nodes"
            section = :nodes
        elseif line == "# edges"
            section = :edges
        elseif startswith(line, "cuts ")
            parts = split(line)               # ["cuts", c1, c2] — c1/c2 may be negative
            cuts = (parse(Int, parts[2]), parse(Int, parts[3]))
        elseif section === :nodes && startswith(line, "n ")
            push!(nodes, fmt.line_to_node(line))
        elseif section === :edges && startswith(line, "e ")
            parts = split(line)              # ["e", colour, portA, portB]
            push!(edges, Edge(parse(Int, parts[2]),
                              _token_to_port(parts[3]), _token_to_port(parts[4])))
        end
    end
    return fmt.build(word, nodes, edges), cuts
end

# ---- morphism = graph + cuts -------------------------------------------------

"""
    morphism_to_text(fmt::GraphTextFormat, m) -> String

The graph text of `m.graph` plus the trailing `cuts <cut1> <cut2>` line.
"""
morphism_to_text(fmt::GraphTextFormat, m) =
    graph_to_text(fmt, m.graph) * "cuts $(m.cut1) $(m.cut2)\n"

"""
    morphism_from_text(fmt::GraphTextFormat, s, ctor, fname) -> ctor(...)

Parse `s` and hand `(graph, cut1, cut2)` to `ctor`. A `cuts` line is REQUIRED — a
text with none is an ERROR, not silently read as `(0, 0)`, since that would
misrepresent "no cuts recorded" as the very meaningful "whole boundary → ε" cut.
`fname` only names the caller in the error message.
"""
function morphism_from_text(fmt::GraphTextFormat, s::AbstractString, ctor, fname::AbstractString)
    graph, cuts = graph_from_text(fmt, s)
    cuts === nothing && error("$fname: no `cuts` line found — " *
                              "this text has no recorded morphism cuts")
    return ctor(graph, cuts[1], cuts[2])
end

# ---- the `---`-separated multi-record container ------------------------------
#
# One free-form metadata COMMENT line per record (preserved VERBATIM, never parsed
# for meaning — e.g. `# e=1,1,0 f=1,0,1 z=1 degree=2` or `# zamo=1,8`), then the
# record body, records separated by a bare `---` line.

# a metadata line always starts with "# " in the file, but callers may pass the
# comment with or without the leading marker — normalise on write, strip on read.
_wgm_meta_line(meta::AbstractString) =
    startswith(meta, "#") ? rstrip(meta) : "# " * meta
_wgm_strip_meta(line::AbstractString) =
    startswith(line, "# ") ? line[3:end] : (startswith(line, "#") ? line[2:end] : line)

"""
    save_entry_file(path, entries, to_text)

Write `entries` (each with `meta` and `morphism`) to `path`, one metadata comment
line per record, records separated by `---`. Order is preserved.
"""
function save_entry_file(path::AbstractString, entries, to_text)
    open(path, "w") do io
        for (i, entry) in enumerate(entries)
            println(io, _wgm_meta_line(entry.meta))
            print(io, to_text(entry.morphism))
            i < length(entries) && println(io, "---")
        end
    end
end

"""
    load_entry_file(path, EntryType, from_text) -> Vector{EntryType}

Read a container written by [`save_entry_file`](@ref) back into its entries, in
order. A stray empty trailing block is skipped.
"""
function load_entry_file(path::AbstractString, ::Type{E}, from_text) where {E}
    text = read(path, String)
    out = E[]
    for block in split(text, "\n---\n")
        lines = split(block, '\n')
        # the first non-blank line is the metadata comment; the rest is the body.
        first_idx = findfirst(l -> !isempty(strip(l)), lines)
        first_idx === nothing && continue          # skip a stray empty trailing block
        meta = _wgm_strip_meta(strip(lines[first_idx]))
        body = join(lines[(first_idx + 1):end], '\n')
        push!(out, E(meta, from_text(body)))
    end
    return out
end
