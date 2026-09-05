# CORE — the shared diagram/word layer under the circular pipeline.
# Checkpoint.jl
#
# Saving/loading a build to disk, so a long run can be time-bounded and RESUMED.
#
# Format choice: we follow the SPIRIT of Oscar's serialization (a versioned,
# human-readable, self-describing text file — Oscar writes JSON `.mrdi` with a
# version header) but implement it in a few dependency-free lines, because our
# data is tiny and simple (words over {1,2,3} + integer edge labels). The file is
# plain text, one record per line, so it diffs nicely in git and can be eyeballed.
#
# A CircularWord is written as its compact string ("1213", or "eps" for ε — we
# avoid the unicode ε on disk for portability). An edge is "v>w:label".
#
# FILE LAYOUT (a .cwg = "circular-word graph" file):
#
#     #cwg 1                         format version
#     maxlen 12                      the length cap
#     status timeout                 :done / :timeout / :maxsteps at save time
#     # nodes
#     N eps
#     N 11
#     ...
#     # frontier                     words still to process, one per line
#     F 1213
#     ...
#     # edges                        directed, with label in {-1,0,1}
#     E eps>11:1
#     E 11>eps:-1
#     ...

const CWG_VERSION = 1

# ---- word <-> on-disk string ----

_word_to_token(w::CircularWord) = isempty(w) ? "eps" : w.str

function _token_to_word(tok::AbstractString)
    tok == "eps" && return EMPTY
    return CircularWord([Int(c - '0') for c in tok])
end

"""
    save_build(path, st; status = :unknown)

Write the full build state `st` (nodes, frontier, edges) to `path` in the plain
text `.cwg` format. `status` records why the build stopped, for information.
"""
function save_build(path::AbstractString, st::BuildState; status::Symbol = :unknown)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "#cwg ", CWG_VERSION)
        println(io, "maxlen ", st.maxlen)
        println(io, "status ", status)

        println(io, "# nodes")
        for w in sort!(collect(keys(st.adj)))
            println(io, "N ", _word_to_token(w))
        end

        println(io, "# frontier")
        # flatten buckets low length -> high, matching processing order
        for bucket in st.frontier, w in bucket
            println(io, "F ", _word_to_token(w))
        end

        println(io, "# edges")
        for v in sort!(collect(keys(st.adj)))
            for (w, lab) in st.adj[v]
                println(io, "E ", _word_to_token(v), ">", _word_to_token(w), ":", lab)
            end
        end
    end
    return path
end

"""
    load_build(path) -> BuildState

Reconstruct a `BuildState` from a `.cwg` file written by `save_build`. The result
can be passed straight to `run_build!` to continue where it left off.
"""
function load_build(path::AbstractString)
    maxlen = nothing
    adj = Dict{CircularWord,Vector{Tuple{CircularWord,Int}}}()
    seen = Set{CircularWord}()
    frontier_words = CircularWord[]

    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "#cwg")
            v = parse(Int, split(line)[2])
            v == CWG_VERSION || error("unsupported .cwg version $v (expected $CWG_VERSION)")
        elseif startswith(line, "maxlen ")
            maxlen = parse(Int, split(line)[2])
        elseif startswith(line, "status ") || startswith(line, "#")
            continue
        elseif startswith(line, "N ")
            w = _token_to_word(split(line)[2])
            push!(seen, w)
            get!(adj, w, Tuple{CircularWord,Int}[])
        elseif startswith(line, "F ")
            push!(frontier_words, _token_to_word(split(line)[2]))
        elseif startswith(line, "E ")
            spec = split(line)[2]                 # "v>w:label"
            vw, labstr = split(spec, ":")
            vtok, wtok = split(vw, ">")
            v = _token_to_word(vtok); w = _token_to_word(wtok)
            push!(get!(adj, v, Tuple{CircularWord,Int}[]), (w, parse(Int, labstr)))
        end
    end
    maxlen === nothing && error("no maxlen in $path")

    st = BuildState(maxlen)
    st.adj = adj
    st.seen = seen
    for w in frontier_words
        push!(st.frontier[length(w) + 1], w)
    end
    return st
end

"""
    build_or_resume(; datafile, maxlen = 12, time_limit_s = Inf, save_every = true)
        -> (SequenceGraph, Symbol)

Convenience driver for long runs. If `datafile` exists, resume from it; otherwise
start fresh (seeded at ε). Run for at most `time_limit_s` seconds, then (if
`save_every`) checkpoint back to `datafile`. Returns the graph built up to that
point and
the stop status (`:done` / `:timeout` / `:maxsteps`).

Typical use: call repeatedly with a time budget; each call makes progress and
saves. When it returns `:done`, exploration (within `maxlen`) is complete.
"""
function build_or_resume(; datafile::AbstractString, maxlen::Int = 12,
                         time_limit_s::Real = Inf, save_every::Bool = true)
    st = isfile(datafile) ? load_build(datafile) : (let s = BuildState(maxlen); _seed_empty!(s); s end)
    status = run_build!(st; time_limit_s = time_limit_s)
    save_every && save_build(datafile, st; status = status)
    return SequenceGraph(st.maxlen, st.adj), status
end

# ---------------------------------------------------------------------------
# Saving a plain SequenceGraph (e.g. a down-graph) — no frontier, just nodes+edges
# ---------------------------------------------------------------------------

"""
    save_graph(path, g)

Write a `SequenceGraph` (nodes + directed labeled edges, no frontier) to `path`
in the same plain-text `.cwg` format used by `save_build` (with an empty frontier
and `status graph`). Round-trips via `load_graph`.
"""
function save_graph(path::AbstractString, g::SequenceGraph)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "#cwg ", CWG_VERSION)
        println(io, "maxlen ", g.maxlen)
        println(io, "status graph")
        println(io, "# nodes")
        for w in sort!(collect(keys(g.adj)))
            println(io, "N ", _word_to_token(w))
        end
        println(io, "# frontier")            # intentionally empty for a plain graph
        println(io, "# edges")
        for v in sort!(collect(keys(g.adj)))
            for (w, lab) in g.adj[v]
                println(io, "E ", _word_to_token(v), ">", _word_to_token(w), ":", lab)
            end
        end
    end
    return path
end

"""
    load_graph(path) -> SequenceGraph

Load a `SequenceGraph` written by `save_graph` (or the graph part of any `.cwg`).
The frontier, if present, is ignored.
"""
function load_graph(path::AbstractString)
    st = load_build(path)
    return SequenceGraph(st.maxlen, st.adj)
end

# ---------------------------------------------------------------------------
# Saving a plain WORD SET (e.g. the upward component) — nodes only, no edges
# ---------------------------------------------------------------------------
#
# `up_component` returns a Set{CircularWord} (membership is all we need), which can
# be a million words. We store just the nodes (one `N <token>` per line) plus the
# caps/status header — far smaller and faster than a full graph with edges.

"""
    save_wordset(path, words; maxlen, status = :unknown)

Write a set/collection of `CircularWord`s to `path` in the plain-text `.cws`
("circular-word set") format: a `#cws 1` header, `maxlen`, `status`, then one
`N <token>` line per word (sorted). Round-trips via `load_wordset`.
"""
function save_wordset(path::AbstractString, words; maxlen::Int, status::Symbol = :unknown)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "#cws ", CWG_VERSION)
        println(io, "maxlen ", maxlen)
        println(io, "status ", status)
        println(io, "# nodes")
        for w in sort!(collect(words))
            println(io, "N ", _word_to_token(w))
        end
    end
    return path
end

"""
    load_wordset(path) -> (words::Set{CircularWord}, maxlen::Int, status::Symbol)

Load a word set written by `save_wordset`. Returns the set together with the
`maxlen` and `status` recorded in the header (so a caller can check the cached
file was built to a large enough length before trusting membership).
"""
function load_wordset(path::AbstractString)
    words = Set{CircularWord}()
    maxlen = nothing
    status = :unknown
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "#cws")
            v = parse(Int, split(line)[2])
            v == CWG_VERSION || error("unsupported .cws version $v (expected $CWG_VERSION)")
        elseif startswith(line, "maxlen ")
            maxlen = parse(Int, split(line)[2])
        elseif startswith(line, "status ")
            status = Symbol(split(line)[2])
        elseif startswith(line, "N ")
            push!(words, _token_to_word(split(line)[2]))
        end
    end
    maxlen === nothing && error("no maxlen in $path")
    return words, maxlen, status
end

"""
    up_component_cached(; maxlen = 20, max_words = 1_000_000,
                        datafile = "data/up_component_len\$maxlen.cws",
                        rebuild = false) -> (words, status)

Build the upward component (`up_component`) OR load it from `datafile` if a cached
file built to at least `maxlen` exists (unless `rebuild = true`). After a fresh
build it saves to `datafile`, so subsequent runs are instant. Returns the word set
and the build status (`:done` / `:maxwords`, or `:cached` when loaded from disk).
"""
function up_component_cached(; maxlen::Int = 20, max_words::Int = 1_000_000,
                             datafile::AbstractString = "",
                             rebuild::Bool = false)
    isempty(datafile) && (datafile = "data/up_component_len$(maxlen).cws")
    if !rebuild && isfile(datafile)
        words, ml, st = load_wordset(datafile)
        if ml >= maxlen && st != :maxwords
            # cached set was built to at least this length and drained fully;
            # restrict to the requested length in case the file is bigger.
            ml > maxlen && (words = Set(w for w in words if length(w) <= maxlen))
            return words, :cached
        end
    end
    words, status = up_component(maxlen = maxlen, max_words = max_words)
    save_wordset(datafile, words; maxlen = maxlen, status = status)
    return words, status
end
