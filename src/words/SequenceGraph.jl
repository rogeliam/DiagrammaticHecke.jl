# CORE — the shared diagram/word layer under the circular pipeline.
# SequenceGraph.jl
#
# The reduction / reachability graph of `CircularWord`s under the rules in Rules.jl.
#
# GOAL. Understand which circular words can be reduced to the empty word ε, in
# type A3. We build this BACKWARDS: start from ε and repeatedly apply moves,
# discovering everything in ε's connected component. Every word we reach is, by
# construction, connected to ε — and since the moves are reversible, "connected
# to ε" == "reducible to ε".
#
# DIRECTED graph. Moves are reversible but their labels flip by direction: a
# reducing move (-1) reversed is an increasing move (+1). So we cannot use one
# undirected label. We store a directed graph: if v --n--> w is a move, we also
# know w --(-n)--> v. We record BOTH directed edges (the reverse is discovered
# when we process w, or added explicitly), each with its own label in {-1,0,+1}.
#
# LENGTH CAP. First pass: keep only words of length <= `maxlen` (default 12).
# When processing a word of length k we still COMPUTE increasing moves, but if a
# neighbor exceeds `maxlen` we skip that edge entirely (its far endpoint is out
# of scope). Reducing/neutral moves stay in scope. This yields the COMPLETE
# induced subgraph on words of length <= maxlen — with the caveat (to check
# later) that a length<=12 word reachable from ε only via a detour above 12 would
# be missed by this cap.
#
# FRONTIER ORDER. We process by increasing length ("take the lowest length whose
# bucket is nonempty"), as requested, using a bucket per length.

"""
    SequenceGraph

The explored graph. Fields:
  - `maxlen`   :: Int                              length cap used
  - `adj`      :: Dict{CircularWord, Vector{Tuple{CircularWord,Int}}}
                   adjacency: node -> list of (neighbor, label) directed edges.
                   Label ∈ {-1,0,+1} is the sign of the move node -> neighbor.
  - `discovered_len` :: Dict{CircularWord,Int}     length bucket a node was found in
All nodes in `adj` are reducible to ε (they are in ε's component, within the cap).
"""
struct SequenceGraph
    maxlen::Int
    adj::Dict{CircularWord,Vector{Tuple{CircularWord,Int}}}
end

nodes(g::SequenceGraph) = keys(g.adj)
edges_from(g::SequenceGraph, w::CircularWord) = get(g.adj, w, Tuple{CircularWord,Int}[])

# a compact one-line summary instead of dumping the whole adjacency dict
function Base.show(io::IO, g::SequenceGraph)
    nn = length(g.adj)
    ne = sum(length(v) for v in values(g.adj); init = 0)
    print(io, "SequenceGraph(maxlen=", g.maxlen, ", ", nn, " nodes, ", ne, " edges)")
end

"""
    edge_label(g, v, w) -> Union{Int,Nothing}

The label of the directed edge v -> w, or `nothing` if absent.
"""
function edge_label(g::SequenceGraph, v::CircularWord, w::CircularWord)
    for (nb, lab) in edges_from(g, v)
        nb == w && return lab
    end
    return nothing
end

# add a directed edge v -(label)-> w (dedup: keep first label seen)
function _add_edge!(adj, v::CircularWord, w::CircularWord, label::Int)
    lst = get!(adj, v, Tuple{CircularWord,Int}[])
    for (nb, _) in lst
        nb == w && return          # already present
    end
    push!(lst, (w, label))
end

"""
    build_reduction_graph(; maxlen = 12) -> SequenceGraph

Backward BFS from ε, exploring ε's connected component restricted to words of
length <= `maxlen`. Returns the directed, labeled graph. Every node is reducible
to ε.

Algorithm: length buckets `frontier[k]` hold discovered-but-unprocessed words of
length k. Repeatedly take the smallest nonempty k, pop a word, enumerate its
one-step neighbors; for each neighbor within the cap add the directed edge (with
its sign) and, if the neighbor is new, drop it into its length bucket.
"""
function build_reduction_graph(; maxlen::Int = 12)
    st = BuildState(maxlen)
    _seed_empty!(st)
    run_build!(st)                     # no time limit
    return SequenceGraph(st.maxlen, st.adj)
end

# ---------------------------------------------------------------------------
# Resumable / time-bounded build
# ---------------------------------------------------------------------------
#
# `BuildState` is the full mutable state of the BFS: the graph built up to that
# point, the set of
# discovered nodes, and the length-bucketed frontier of words still to process.
# It can be checkpointed (see Checkpoint.jl) and resumed exactly.

mutable struct BuildState
    maxlen::Int
    adj::Dict{CircularWord,Vector{Tuple{CircularWord,Int}}}
    seen::Set{CircularWord}
    frontier::Vector{Vector{CircularWord}}   # bucket k+1 holds words of length k
end

BuildState(maxlen::Int) = BuildState(
    maxlen,
    Dict{CircularWord,Vector{Tuple{CircularWord,Int}}}(),
    Set{CircularWord}(),
    [Vector{CircularWord}() for _ in 0:maxlen],
)

_bucket(st::BuildState, k::Int) = st.frontier[k + 1]

# seed the BFS with the empty word (only if not already present)
function _seed_empty!(st::BuildState)
    EMPTY in st.seen && return st
    push!(st.seen, EMPTY)
    st.adj[EMPTY] = Tuple{CircularWord,Int}[]
    push!(_bucket(st, 0), EMPTY)
    return st
end

"""
    run_build!(st; time_limit_s = Inf, max_steps = typemax(Int)) -> Symbol

Advance the BFS in `st` (processing lowest-length words first) until the frontier
is empty, or the wall-clock `time_limit_s` is exceeded, or `max_steps` words have
been processed. Returns why it stopped: `:done`, `:timeout`, or `:maxsteps`.

On `:timeout`/`:maxsteps` the state is fully consistent and can be checkpointed
and resumed later — the remaining frontier lists exactly which words still need
processing.
"""
function run_build!(st::BuildState; time_limit_s::Real = Inf,
                    max_steps::Integer = typemax(Int))
    t0 = time()
    steps = 0
    while true
        k = findfirst(!isempty, st.frontier)
        k === nothing && return :done
        steps >= max_steps && return :maxsteps
        (time() - t0) >= time_limit_s && return :timeout

        w = pop!(st.frontier[k])       # k is 1-based bucket index = length+1
        steps += 1
        for (nb, sgn) in neighbors(w)
            length(nb) > st.maxlen && continue    # out of scope: skip edge
            _add_edge!(st.adj, w, nb, sgn)
            _add_edge!(st.adj, nb, w, -sgn)       # reverse edge, flipped label
            if !(nb in st.seen)
                push!(st.seen, nb)
                get!(st.adj, nb, Tuple{CircularWord,Int}[])
                push!(_bucket(st, length(nb)), nb)
            end
        end
    end
end

"""
    frontier_remaining(st) -> Int

How many discovered-but-unprocessed words remain (across all length buckets).
Zero means the component (within the cap) is fully explored.
"""
frontier_remaining(st::BuildState) = sum(length, st.frontier; init = 0)

# ---------------------------------------------------------------------------
# Session cache + size guard
# ---------------------------------------------------------------------------
#
# The graph is a purely LOCAL computation (every move depends only on the word
# itself), so there is no global state to share BETWEEN runs — but within one
# session it is wasteful to rebuild the same graph twice. We memoize by `maxlen`
# in a module-global dict that lives only as long as the Julia session.
#
# Size guard: a runaway `maxlen` could build a graph that eats all RAM. Before
# accepting a freshly built graph into the cache we estimate its byte size; if it
# exceeds `size_limit_bytes` (default 10 MB) we DON'T cache it (and warn), so the
# cache can never grow without bound. The build itself still returns the graph —
# the cap governs caching, not computation. Pass a larger limit if you really
# want a big graph cached.

const _GRAPH_CACHE = Dict{Int,SequenceGraph}()
const DEFAULT_SIZE_LIMIT_BYTES = 10 * 1024 * 1024   # 10 MB

"""
    graph_bytes(g) -> Int

Rough in-memory size of a `SequenceGraph`: counts the words' letters and the edge
tuples. An estimate, not exact, but good enough for a soft cap.
"""
function graph_bytes(g::SequenceGraph)
    total = 0
    for (w, lst) in g.adj
        total += length(w) + 16              # word letters + small overhead
        total += length(lst) * 24            # (word ptr, Int) tuples, rough
    end
    return total
end

"""
    cached_graph(; maxlen = 12, size_limit_bytes = DEFAULT_SIZE_LIMIT_BYTES,
                 rebuild = false) -> SequenceGraph

Session-memoized `build_reduction_graph`. Returns the cached graph for `maxlen`
if present (unless `rebuild = true`); otherwise builds it, and caches it **only
if** its estimated size is within `size_limit_bytes` (warns and skips caching if
not). Clear the cache with `clear_graph_cache!()`.
"""
function cached_graph(; maxlen::Int = 12,
                      size_limit_bytes::Int = DEFAULT_SIZE_LIMIT_BYTES,
                      rebuild::Bool = false)
    if !rebuild && haskey(_GRAPH_CACHE, maxlen)
        return _GRAPH_CACHE[maxlen]
    end
    g = build_reduction_graph(maxlen = maxlen)
    b = graph_bytes(g)
    if b <= size_limit_bytes
        _GRAPH_CACHE[maxlen] = g
    else
        @warn "graph for maxlen=$maxlen is ~$(round(b/1024/1024, digits=1)) MB > \
               limit $(round(size_limit_bytes/1024/1024, digits=1)) MB; not cached"
    end
    return g
end

"Empty the session graph cache."
clear_graph_cache!() = (empty!(_GRAPH_CACHE); nothing)

"""
    reducible_words(g) -> Vector{CircularWord}

All words in the graph (== all words of length <= g.maxlen reducible to ε, modulo
the length-cap caveat), sorted by (length, lexicographic normal form).
"""
reducible_words(g::SequenceGraph) = sort!(collect(nodes(g)))

"""
    all_labels_nonincreasing(g) -> Bool

True iff every directed edge in `g` has label <= 0.
"""
all_labels_nonincreasing(g::SequenceGraph) =
    all(lab <= 0 for w in nodes(g) for (_, lab) in edges_from(g, w))

"""
    non_increasing_subgraph(g) -> SequenceGraph

The subgraph keeping only NON-INCREASING directed edges, i.e. labels in {-1, 0}
(reducing and neutral moves). Same node set. This is the graph in which we ask
"can we get down to ε using only reducing/neutral steps?".
"""
function non_increasing_subgraph(g::SequenceGraph)
    adj = Dict{CircularWord,Vector{Tuple{CircularWord,Int}}}()
    for w in nodes(g)
        adj[w] = Tuple{CircularWord,Int}[]   # keep isolated nodes too
    end
    for (v, lst) in g.adj
        for (w, lab) in lst
            if lab <= 0
                push!(adj[v], (w, lab))
            end
        end
    end
    return SequenceGraph(g.maxlen, adj)
end

# ---------------------------------------------------------------------------
# The "down graph": non-increasing moves under the EXTENDED rule set
# ---------------------------------------------------------------------------
#
# Second experiment (extended rules). Take a set of seed words. Build the graph
# whose edges are the NON-INCREASING moves of the extended rule set (all reducers
# BASE ∪ EXTRA applied downward, plus neutral braid/commutation), via
# `down_neighbors`. Then ask: from each seed, can we reach ε using only these
# edges? Because every edge has length change ≤ 0, the reachable set from any seed
# is finite (length can only drop or stay the same, and equal-length neutral
# classes are finite), so this terminates.

"""
    down_graph(seeds) -> SequenceGraph

Build the directed graph of non-increasing extended-rule moves (`down_neighbors`)
reachable from the given `seeds` (any iterable of `CircularWord`). Labels are −1
(reducing) or 0 (neutral). New words discovered along the way (always of length ≤
the seed that produced them) are included. `maxlen` of the result is the largest
seed length seen.
"""
function down_graph(seeds)
    adj = Dict{CircularWord,Vector{Tuple{CircularWord,Int}}}()
    stack = CircularWord[]
    maxlen = 0
    for s in seeds
        if !haskey(adj, s)
            adj[s] = Tuple{CircularWord,Int}[]
            push!(stack, s)
            maxlen = max(maxlen, length(s))
        end
    end
    while !isempty(stack)
        w = pop!(stack)
        outs = down_neighbors(w)
        adj[w] = outs
        for (nb, _) in outs
            if !haskey(adj, nb)
                adj[nb] = Tuple{CircularWord,Int}[]
                push!(stack, nb)
                maxlen = max(maxlen, length(nb))
            end
        end
    end
    return SequenceGraph(maxlen, adj)
end

"""
    reaches_empty_down(dg) -> Set{CircularWord}

Given a down-graph `dg` (from `down_graph`), the set of nodes that can reach ε
following directed edges. Computed by reverse reachability from ε: reverse all
edges, then BFS from ε.
"""
function reaches_empty_down(dg::SequenceGraph)
    # reverse adjacency
    rev = Dict{CircularWord,Vector{CircularWord}}()
    for w in nodes(dg); rev[w] = CircularWord[]; end
    for (v, lst) in dg.adj, (w, _) in lst
        haskey(rev, w) || (rev[w] = CircularWord[])
        push!(rev[w], v)
    end
    reach = Set{CircularWord}()
    EMPTY in nodes(dg) || return reach          # ε not even present
    q = CircularWord[EMPTY]; push!(reach, EMPTY)
    while !isempty(q)
        x = popfirst!(q)
        for pre in get(rev, x, CircularWord[])
            if !(pre in reach)
                push!(reach, pre); push!(q, pre)
            end
        end
    end
    return reach
end

"""
    reduces_to_empty_down(w; memo) -> Bool

Whether `w` can reach ε using only non-increasing extended moves (`down_neighbors`),
computed by a memoized search — each distinct word is visited ONCE and its answer
cached in `memo` (a `Dict{CircularWord,Bool}`, shared across calls for speed).

This avoids materialising one giant down-graph: we recurse through `down_neighbors`
with memoization. Because every move is non-increasing, the recursion is finite
(length can only drop or stay equal; equal-length neutral classes are finite and
guarded by an in-progress marker to avoid cycles).
"""
function reduces_to_empty_down(w::CircularWord;
                               memo::Dict{CircularWord,Bool} = Dict{CircularWord,Bool}())
    isempty(w) && return true
    haskey(memo, w) && return memo[w]

    # Explore the NEUTRAL-connected class of `w` (same-length words reachable by
    # 0-labelled moves). `w` reduces to ε iff SOME word in this class has a
    # strictly-reducing (−1) move to a word that itself reduces to ε. This handles
    # neutral cycles correctly (no premature false-memoization inside the class).
    class = Set{CircularWord}([w])
    stack = CircularWord[w]
    reducing_targets = CircularWord[]     # −1 neighbors leaving the class (shorter)
    while !isempty(stack)
        x = pop!(stack)
        for (nb, s) in down_neighbors(x)
            if s == 0
                if !(nb in class); push!(class, nb); push!(stack, nb); end
            else                          # s == -1 (down_neighbors excludes +1)
                push!(reducing_targets, nb)
            end
        end
    end

    ans = false
    for t in reducing_targets
        if reduces_to_empty_down(t; memo = memo)
            ans = true; break
        end
    end
    # memoize the whole neutral class with the same answer
    for c in class; memo[c] = ans; end
    return ans
end

"""
    down_path(w; maxsteps = 100_000) -> Union{Nothing, Vector{Tuple{CircularWord,Int}}}

A concrete non-increasing path from `w` down to ε, or `nothing` if `w` can't
reduce. Returned as a list of `(word, sign)` steps STARTING at `w` and ending at
ε, where `sign` is the label of the move that LEFT the previous word (the first
entry's sign is 0, a placeholder for the start). BFS over `down_neighbors`, so the
path is shortest in number of moves.
"""
function down_path(w::CircularWord; maxsteps::Int = 100_000)
    isempty(w) && return [(EMPTY, 0)]
    prev = Dict{CircularWord,Tuple{CircularWord,Int}}()   # node => (from, sign)
    seen = Set{CircularWord}([w])
    q = CircularWord[w]
    steps = 0
    found = false
    while !isempty(q) && steps < maxsteps
        x = popfirst!(q); steps += 1
        x == EMPTY && (found = true; break)
        for (nb, s) in down_neighbors(x)
            if !(nb in seen)
                push!(seen, nb); prev[nb] = (x, s); push!(q, nb)
                nb == EMPTY && (found = true)
            end
        end
        found && break
    end
    (EMPTY in seen) || return nothing
    # reconstruct ε <- … <- w
    path = Tuple{CircularWord,Int}[]
    cur = EMPTY
    while cur != w
        haskey(prev, cur) || return nothing
        p, s = prev[cur]
        pushfirst!(path, (cur, s))
        cur = p
    end
    pushfirst!(path, (w, 0))
    return path
end

"""
    down_reachability_report(seeds) -> (dg, reach, stuck)

Convenience: build the down-graph from `seeds`, compute which nodes reach ε
(`reach`) and which cannot (`stuck` = nodes not in `reach`). Returns the graph and
both sets. `stuck` being empty means EVERY seed reduces to ε using only
non-increasing extended-rule moves.
"""
function down_reachability_report(seeds)
    dg = down_graph(seeds)
    reach = reaches_empty_down(dg)
    stuck = Set(w for w in nodes(dg) if !(w in reach))
    return dg, reach, stuck
end

"""
    down_counts_by_length(seeds) -> Vector{NTuple{3,Int}}

Run the down-analysis on `seeds` and return, per length `k`, a tuple
`(k, n_reduce, n_stuck)`: how many words of length `k` in the down-graph can
reduce to ε and how many cannot. Only lengths that occur are returned, sorted.
"""
function down_counts_by_length(seeds)
    dg = down_graph(seeds)
    reach = reaches_empty_down(dg)
    red = Dict{Int,Int}(); stk = Dict{Int,Int}()
    for w in nodes(dg)
        k = length(w)
        if w in reach
            red[k] = get(red, k, 0) + 1
        else
            stk[k] = get(stk, k, 0) + 1
        end
    end
    ks = sort!(collect(union(keys(red), keys(stk))))
    return [(k, get(red, k, 0), get(stk, k, 0)) for k in ks]
end

# ---------------------------------------------------------------------------
# "Expand then reduce": increase the whole word, then come down to ε
# ---------------------------------------------------------------------------
#
# The property: a word `w` such that SOME full increasing expansion of `w` (every
# letter covered by an increasing move — see `increasing_expansions`) is reducible
# to ε using only non-increasing moves. The `1221` example:
#   1221 --(12->2121, 21->1212)--> 21211212 --down--> ε.
#
# QUICK CHECK: we reuse a precomputed oracle `reach` = the set of words that reach
# ε down-only (from `reaches_empty_down` on a down-graph). `w` has the property
# iff any expansion lands in `reach`.

"""
    expands_then_reduces(w, reach) -> Bool

True iff some full increasing expansion of `w` (all positions covered) is in the
set `reach` of down-reducible words. `reach` is typically
`reaches_empty_down(down_graph(seeds))` — the oracle computed once and reused, so
each word is a quick membership test.
"""
function expands_then_reduces(w::CircularWord, reach::AbstractSet{CircularWord})
    for ex in increasing_expansions(w)
        ex in reach && return true
    end
    return false
end

"""
    find_expand_then_reduce(words, reach) -> Vector{CircularWord}

Among `words`, those satisfying `expands_then_reduces(w, reach)`, sorted. `reach`
is the down-reducible oracle set (reuse the down-analysis result for speed).
"""
function find_expand_then_reduce(words, reach::AbstractSet{CircularWord})
    hits = [w for w in words if expands_then_reduces(w, reach)]
    return sort!(hits)
end

# ---------------------------------------------------------------------------
# The UPWARD component: ε, then only increasing moves (no ε-exit after step 1)
# ---------------------------------------------------------------------------
#
# Read FORWARD from ε as an upward walk:
#   - step 1 (from ε): an ε-exit — ε→11, ε→22, ε→121212 (drop ε→33, ε→232323 by
#     the 1↔3 symmetry). These are the SEEDS.
#   - every later step: a strictly-INCREASING move whose source is NONEMPTY
#     (`up_neighbors`, i.e. UP_RULES) — 1→11, 21→1212, 2→12121, 3→131, … . NO
#     decreasing moves, NO neutral moves, and NO further ε-exit.
#
# So this is the set of circular words reachable from ε by going strictly up the
# whole way with a single initial blow-up. Equivalently (reversing every edge) the
# words that reduce to ε using only these moves downward. We stop at length
# `maxlen` OR once `max_words` nodes have been collected, whichever comes first —
# the walk is length-monotone (every edge is +Δ), so a length cap makes it finite;
# the node cap guards against an explosion before the length cap bites.

"the non-ε seeds of the allowed ε-exits (11, 22, 121212; 33/232323 dropped by 1↔3 symmetry)"
const UP_SEEDS = CircularWord[CircularWord([1,1]), CircularWord([2,2]),
                              CircularWord([1,2,1,2,1,2])]

"""
    up_component(; maxlen = 20, max_words = 1_000_000, seeds = UP_SEEDS)
        -> (words::Set{CircularWord}, status::Symbol)

Enumerate the UPWARD component from ε: start at the `seeds` (the allowed ε-exits),
then repeatedly apply `up_neighbors` (strictly-increasing, nonempty-source moves —
no ε-exit, no decrease), collecting every word of length ≤ `maxlen`. BFS by
increasing length. Stops when the frontier is exhausted (`:done`) or when
`max_words` distinct words have been collected (`:maxwords`).

Returns the SET of words (not a graph — we only need membership/counts here, and a
full adjacency of ~10⁶ nodes is wasteful). Use `counts_by_length(collect(words))`
for the per-length breakdown.
"""
function up_component(; maxlen::Int = 20, max_words::Int = 1_000_000,
                      seeds::AbstractVector{CircularWord} = UP_SEEDS,
                      progress_every::Int = 0)
    seen = Set{CircularWord}()
    frontier = [Vector{CircularWord}() for _ in 0:maxlen]
    for s in seeds
        length(s) > maxlen && continue
        if !(s in seen)
            push!(seen, s)
            push!(frontier[length(s) + 1], s)
        end
    end
    next_report = progress_every
    while true
        if length(seen) >= max_words
            progress_every > 0 && @info "up_component: reached word cap" words=length(seen)
            return seen, :maxwords
        end
        k = findfirst(!isempty, frontier)
        k === nothing && return seen, :done
        w = pop!(frontier[k])
        for nb in up_neighbors(w)
            length(nb) > maxlen && continue
            if !(nb in seen)
                push!(seen, nb)
                push!(frontier[length(nb) + 1], nb)
                if progress_every > 0 && length(seen) >= next_report
                    println("  up_component: ", length(seen), " words so far ",
                            "(current length ", length(w), ")")
                    flush(stdout)
                    next_report += progress_every
                end
                length(seen) >= max_words && return seen, :maxwords
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Enumerating circular words + the "expandable all around" list
# ---------------------------------------------------------------------------

"""
    all_circular_words(maxlen) -> Vector{CircularWord}

Every circular word over {1,2,3} of length 1..`maxlen` (normal forms, deduped),
plus ε. Enumerated by generating all linear words and normalizing; deduplication
by normal form collapses rotations. (Fine for maxlen ≤ ~14; 3^n before dedup.)
"""
function all_circular_words(maxlen::Int)
    seen = Set{CircularWord}([EMPTY])
    for n in 1:maxlen
        for code in 0:(3^n - 1)          # base-3 enumeration of length-n words
            v = Vector{Int}(undef, n)
            c = code
            for i in 1:n
                v[i] = (c % 3) + 1; c ÷= 3
            end
            push!(seen, CircularWord(v))
        end
    end
    return collect(seen)
end

"""
    expandable_all_around(maxlen; require_all_letters = true)
        -> Vector{CircularWord}

The circular words of length ≤ `maxlen` that admit a full braid-expander tiling
(`has_full_expansion`), optionally also requiring that all three letters occur
(`uses_all_letters`). Sorted by (length, normal form).
"""
function expandable_all_around(maxlen::Int; require_all_letters::Bool = true)
    out = CircularWord[]
    for w in all_circular_words(maxlen)
        require_all_letters && !uses_all_letters(w) && continue
        has_full_expansion(w) || continue
        push!(out, w)
    end
    return sort!(out)
end

"""
    counts_by_length(words) -> Vector{Tuple{Int,Int}}

`(length, count)` pairs for a word list, sorted by length.
"""
function counts_by_length(words)
    d = Dict{Int,Int}()
    for w in words; d[length(w)] = get(d, length(w), 0) + 1; end
    return [(k, d[k]) for k in sort!(collect(keys(d)))]
end
