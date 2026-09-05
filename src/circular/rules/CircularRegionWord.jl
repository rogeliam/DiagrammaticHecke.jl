# circular/rules/CircularRegionWord.jl — the REGION-WORD criterion as the trigger
# of the 2parallel rule.
#
# Every region carries a region word. Where a region's word is not reduced, a rule
# must fire. The SHORTEST such place is what is looked for: region_i has reduced
# word w_i, and the next region has word w_{i+1} which is not reduced. In the
# simplest case w_i ends directly on that letter, which is 2parallel. In the
# second-simplest case the word reads `212` before and `1` after. The remaining
# cases are recorded, not acted on.
#
# THE REFORMULATION THAT CARRIES IT: "not reduced" is COXETER LENGTH
# (`is_reduced`, morphism/BraidMoves.jl). `w_i · s` is not reduced iff `s` is a
# right descent of `w_i`. No braid classes needed for the trigger itself.
#
# A TRANSITION is a pair (region i, letter s) with `w_i` reduced and `w_i · s`
# not reduced, where `s` is crossed when leaving region i:
#   * kind = :edge  — an `s`-edge on the boundary of region i leading to a
#                     region one step further from the mark;
#   * kind = :dot   — an `s`-dot INSIDE region i (its cap edge). Dot edges are
#                     transparent for the distances (CircularRegion.jl), so the dot
#                     colour is the only thing that "crosses" here — this is
#                     exactly the `121` + 2-dot pattern.
# The SHORTEST transition wins: smallest `dist[i]`, then region number, then
# edges before dots, then edge index (deterministic; documented here).
#
# THE CASES, all read off the reduced word `w_i`:
#   :direct    w_i ends on s               -> 2parallel at that edge (case 1)
#   :dihedral  ONE braid/commutation move of w_i ends on s (case 2)
#   :long      s is a descent, but only reachable through a longer chain
#              -> recorded only
#
# WHAT IS WIRED WHERE:
#   * Case 1 = condition 4 of 2parallel
#     (Circular2Parallel.jl, `find_circular_2parallel`): a candidate pair is
#     accepted only if it lives in region i, contains the transition edge, and
#     the partner edge leads towards the mark. The strict guard (condition 6)
#     and the connection test are UNCHANGED — so the termination measure of
#     Circular2Parallel.jl (file head) is untouched: `:regionword` only SELECTS among
#     the pairs `:d4` would allow (one side closer than the middle), it never
#     adds one. The default is `:d4`.
#   * Case 2 (dot) = `circular_regionword_dihedral_step`, a SEPARATE entry point
#     (not in the driver): it forces `apply_circular_d4` at an `s`-edge bounding
#     the dot's region.

"""
    CIRCULAR_REGIONWORD_TRIVALENT_START

Does the ε of the start region spread across a TRIVALENT start? Default
**`false`**.

⚠ **NOT A GOOD READING.** Forcing several regions to distance 0 is wrong:
distance 0 belongs to exactly ONE region, the marked one, and everything else is
higher. The test below also asks only for "two edges, same colour" and NOT for
`arm_count == 3`, so it fires at a 10-armed node as well; and the figure it
addresses — the last term of `id_1212` — is covered by
[`CIRCULAR_REGIONWORD_2K`](@ref). The switch stays for comparison runs.

**What it is for**: simplifying the last term of `id_1212` produces a trivalent
merge around the start region, then a 12-braid, then a 2-braid, and the criterion
must still hold there. If the distance-0 region is bounded by 2 edges of the same
colour (i.e. trivalent), its neighbours would also get ε as their word.

That term is the one that made the criterion non-terminating: its start region
is bounded by EXACTLY TWO edges of the SAME colour (the two lower legs of a
trivalent merge), and counting that colour once per leg produced the region word
`1212` all over again. The two legs are the same strand, though — crossing
either of them does not move you anywhere in the Coxeter group. So both
neighbours inherit the EMPTY word (and no parent edge: their path from the mark
is empty), and the BFS continues from there.

Set to `false` for the plain level BFS, for
comparison.
"""
const CIRCULAR_REGIONWORD_TRIVALENT_START = Ref(false)

"""
    _crw_trivalent_start_neighbours(g, adj, start) -> Vector{Int}

The neighbours that inherit the start region's ε
([`CIRCULAR_REGIONWORD_TRIVALENT_START`](@ref)): non-empty only if the boundary
of `start` consists of exactly two edges OF THE SAME COLOUR. Returns the two
regions on the other side (deduplicated — if both edges lead to the same region
it is listed once).
"""
function _crw_trivalent_start_neighbours(g::CircularGraph, adj, start::Int)
    CIRCULAR_REGIONWORD_TRIVALENT_START[] || return Int[]
    b = adj[start]
    length(b) == 2 || return Int[]
    e1, e2 = b[1][2], b[2][2]
    e1 == e2 && return Int[]
    g.edges[e1].colour == g.edges[e2].colour || return Int[]
    return unique(Int[b[1][1], b[2][1]])
end

"""
    circular_region_words(m::CircularMorphismGraph)
        -> (words, parent_edge)

The region words BY DEFINITION: a LEVEL BFS over the region adjacency — distance 0
is the marked region (empty word), every step crosses ONE edge and appends its
colour, each cell gets its word from exactly ONE parent cell plus edge, as in a
tree. Deterministic: FIFO BFS, neighbours in ascending edge index.

This is the pure EDGE distance (`circular_region_distances_edges_only`), NOT the
weighted one with node jumps (`circular_region_distances`). The distinction is
real: the two differ in 235 of 2130 end-term regions of the reference corpora.

Returns `words[R]` (`nothing` if unreachable / no mark) and `parent_edge[R]`
(the tree edge into `R`; `0` for the start region and unreachable regions).
"""
function circular_region_words(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    nreg = region_count(g)
    words = Vector{Union{Nothing, Vector{Int}}}(nothing, nreg)
    parent = zeros(Int, nreg)
    (n == 0 || _circular_left_mark(m) === nothing) && return (words, parent)
    regs = regions(g)
    gap = mod1(m.cut1, n)
    start = findfirst(R -> gap in R.gaps, regs)
    start === nothing && return (words, parent)

    adj, _ = circular_region_adjacency(g)
    words[start] = Int[]
    queue = [start]
    # the trivalent start (see `CIRCULAR_REGIONWORD_TRIVALENT_START`): both
    # neighbours inherit ε, with NO parent edge — their path from the mark is
    # empty, so `circular_path_edges` and the region word stay the same object.
    for R2 in _crw_trivalent_start_neighbours(g, adj, start)
        R2 == start && continue
        words[R2] = Int[]
        parent[R2] = 0
        push!(queue, R2)
    end
    while !isempty(queue)
        R = popfirst!(queue)
        for (R2, ei) in sort(adj[R]; by = x -> x[2])
            words[R2] === nothing || continue
            words[R2] = vcat(words[R]::Vector{Int}, g.edges[ei].colour)
            parent[R2] = ei
            push!(queue, R2)
        end
    end
    return (words, parent)
end

"""
    CIRCULAR_REGIONWORD_MVALENT

Does the region word measure at the M-VALENT NODE instead of edge by edge?
Default **`false`**: under this reading `id_1213213` term 5 fires, so there is
no termination measure yet (see below).

**What it is for.** On `L4`, the leaf whose region word the level BFS reads as
`22,2,2,eps,22`, walking edge by edge builds both distances and the region word
wrongly. At an m-valent vertex all other neighbours are only 2 steps away, so
regions 1 and 5 must carry the word `2` and the reading must be `2,2,2,eps,2`.

L4 has exactly ONE node with more than two arms, `v3` with arms `[2,2,2,2,2]`
(5-valent, `:gen2`), and ALL FIVE regions are its sectors. Walking around that
node edge by edge spells one letter per arm; around an m-valent node every
sector is ONE crossing away. That is the node adjacency (b) of
[`circular_region_distances`](@ref), which the plain level BFS sets aside.

**The reading** ([`circular_region_mvalent_words`](@ref)) is therefore
[`circular_region_tree_words`](@ref) PLUS the ε/fold collapse, the latter
restricted to a **trivalent** shared node:

1. edge steps cost 1 and append the edge colour;
2. node steps at `:gen2`/`:gen13` nodes connect ALL sector regions, cost =
   the distinct colours crossed on the cheaper arc, capped by
   [`node_cap`](@ref) — at a general-2 node exactly one letter, however many
   arms;
3. cost **0** across an edge bounding a FOLD POCKET: a region whose boundary
   is exactly two edges of the SAME colour sharing a node **with three arms**.
   Only there are the two legs one strand (strand + dot arm); two arms of a
   2k-valent node are not — that is precisely what L4 shows, and it is why a fold
   reading over ANY shared node is wrong.

This reading reaches `2,2,2,eps,2` on L4, silences every false alarm of the level
BFS (`id_1212` terms 1 and 4, `id_1213213` term 3), keeps every genuine trigger
(`id_11`, `id_22`, `id_1213213` term 2), and — unlike a plain node-jump reading —
keeps the α₁ terms clean, because the trivalent ε rule carries them. ⚠
`id_1213213` term 5 fires under it, so there is no termination measure yet.

Under this reading exactly ONE of the four degree-2 leaves of `21232/21232`
carries a non-reduced word: **L5**.

Read by [`circular_unreduced_transition`](@ref), which switches to these words
when this is `true`.
"""
const CIRCULAR_REGIONWORD_MVALENT = Ref(false)

"""
    _crw_fold_pockets(g, adj) -> Set{Int}

The edges of cost 0 in [`circular_region_mvalent_words`](@ref): both boundary
edges of every FOLD POCKET — a region bounded by exactly two edges of the same
colour that share a **trivalent** node (`arm_count == 3`).

The arm count matters: without it the two arms of an m-valent node would count as
one folded strand as well, and L4 would collapse to ε everywhere instead of
reading `2,2,2,eps,2`.
"""
function _crw_fold_pockets(g::CircularGraph, adj)
    out = Set{Int}()
    for R in 1:region_count(g)
        b = adj[R]
        length(b) == 2 || continue
        e1, e2 = b[1][2], b[2][2]
        e1 == e2 && continue
        E1, E2 = g.edges[e1], g.edges[e2]
        E1.colour == E2.colour || continue
        n1 = Int[p.node for p in (E1.a, E1.b) if p isa NodePort]
        n2 = Int[p.node for p in (E2.a, E2.b) if p isa NodePort]
        shared = intersect(n1, n2)
        any(v -> arm_count(g.nodes[v]) == 3, shared) || continue
        push!(out, e1); push!(out, e2)
    end
    return out
end

"""
    circular_region_mvalent_words(m::CircularMorphismGraph)
        -> (words, parent_edge)

The region words of the M-VALENT reading (see
[`CIRCULAR_REGIONWORD_MVALENT`](@ref)) — a function of its own next to
[`circular_region_words`](@ref), because the pure edge BFS stays as the comparison
quantity, and so does [`circular_region_tree_words`](@ref), which is this reading
WITHOUT the fold collapse).

Bellman-Ford over the region graph (weights ≠ 1, so no BFS), same step set as
`circular_region_tree_words` — edge steps weight 1, node steps at
`:gen2`/`:gen13` nodes weight = distinct colours of the cheaper arc, capped by
[`node_cap`](@ref) — plus weight **0** across the boundary edges of a trivalent
fold pocket ([`_crw_fold_pockets`](@ref)). Only a GENUINE improvement relaxes,
so the first witness stays and the tree is deterministic (as in
`circular_region_tree_words`).

⚠ `parent_edge[R]` is the crossed EDGE, and it is `0` when the step into `R`
was a NODE JUMP (a jump is no edge at all) as well as for the start region and
unreachable regions. So the `parent_edge` chain is not a full path here, and
the word is not the colour sequence of a chain of edges. What the cut curve of the
rex fusion is under this reading is not yet resolved, so `circular_path_edges` reads the
level BFS.

Returns `(words[R], parent_edge[R])` like `circular_region_words`.
"""
function circular_region_mvalent_words(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    nreg = region_count(g)
    words = Vector{Union{Nothing, Vector{Int}}}(nothing, nreg)
    parent = zeros(Int, nreg)
    (n == 0 || _circular_left_mark(m) === nothing) && return (words, parent)
    regs = regions(g)
    gap = mod1(m.cut1, n)
    start = findfirst(R -> gap in R.gaps, regs)
    start === nothing && return (words, parent)

    adj, edge_between = circular_region_adjacency(g)
    zerocost = _crw_fold_pockets(g, adj)

    # the directed steps `(target, letters, edge)`; `edge == 0` = a node jump.
    steps = [Tuple{Int, Vector{Int}, Int}[] for _ in 1:nreg]

    # (a) edge steps: one edge, one colour — none across a fold pocket wall.
    for R in 1:nreg, (R2, _) in adj[R]
        e = edge_between[(min(R, R2), max(R, R2))]
        push!(steps[R], (R2, e in zerocost ? Int[] : [g.edges[e].colour], e))
    end

    # (b) node steps at :gen2/:gen13 nodes — the m-valent part.
    sec = circular_node_sector_regions(g)
    for (v, nd) in enumerate(g.nodes)
        cap = node_cap(node_class(nd))
        cap === nothing && continue
        wired = [s for s in 1:arm_count(nd) if sec[v][s] != 0]
        k = length(wired)
        k < 2 && continue
        Rs   = [sec[v][s] for s in wired]
        cols = [arm_colour(nd, s) for s in wired]
        allcols = unique(cols)
        for i in 1:k, j in 1:k
            i == j && continue
            fwd = Int[]; t = i; while t != j; t = mod1(t + 1, k); push!(fwd, cols[t]); end
            bwd = Int[]; t = j; while t != i; t = mod1(t + 1, k); push!(bwd, cols[t]); end
            ufwd, ubwd = unique(fwd), unique(bwd)
            best = length(ufwd) <= length(ubwd) ? ufwd : ubwd
            length(best) > cap && (best = allcols)   # the cap: once per colour
            push!(steps[Rs[i]], (Rs[j], best, 0))
        end
    end

    dist = fill(typemax(Int), nreg)
    dist[start] = 0
    words[start] = Int[]
    for _ in 1:nreg, a in 1:nreg          # small graphs: Bellman-Ford is enough
        dist[a] == typemax(Int) && continue
        for (b, letters, e) in steps[a]
            d2 = dist[a] + length(letters)
            if d2 < dist[b]               # only a GENUINE improvement
                dist[b] = d2
                words[b] = vcat(words[a]::Vector{Int}, letters)
                parent[b] = e
            end
        end
    end
    return (words, parent)
end

# ---------------------------------------------------------------------------
# The 2k-SATURATED reading
# ---------------------------------------------------------------------------

"""
    _CRW_M_BRAID

The `m` at which the region word saturates around a `{1,2}`/`{2,3}` node: **3**,
fixed. Other Coxeter types would give a different value; only A₃ is in scope.
"""
const _CRW_M_BRAID = 3

"""
    _crw_node_cap(nd::CircularNode) -> Tuple{Int, Vector{Int}}

How far the region word may grow around ONE node, and the word it saturates to:
the length of the longest element of the parabolic subgroup generated by the
node's colours, written canonically (alternating, starting with the SMALLER
colour) so that both directions around the node give the same string.

| node colours | subgroup | cap | saturated word |
|---|---|---|---|
| `{c}` (`:mono`, single-colour `:mixed`, every trivalent) | `⟨c⟩` | 1 | `c` |
| `{1,3}` (the commuting pair) | `⟨1,3⟩` | 2 | `13` |
| `{1,2}` / `{2,3}` (`:braid`) | dihedral, `m = 3` | 3 | `121` / `232` |

This is the SAME statement for every node class: a node identifies everything
beyond its own longest element, so walking around it can never add more letters
than that. For the `2k`-armed node it reads "all other neighbours have `w·i·j·i`";
in the `cap = 1` instance (`circular_region_mvalent_words`) it reads "all sector
regions are ONE crossing away from each other".
"""
function _crw_node_cap(nd::CircularNode)
    cols = sort(unique(nd.arms))
    length(cols) == 1 && return (1, [cols[1]])
    m = Set(cols) == Set([1, 3]) ? 2 : _CRW_M_BRAID
    return (m, [isodd(k) ? cols[1] : cols[2] for k in 1:m])
end

"""
    CIRCULAR_REGIONWORD_2K

Saturate the region words around a `2k`-armed node, default **`true`**.

**The rule:** the first word touching the node is `w`.
Then the direct neighbours are `w·i` for `i` a matching colour. After two steps
it's `w·i·j`. All other neighbours then have `w·i·j·i`. So even with 20 arms
we get for example `ε, 1, 2, 12, 21`, and then `121` fifteen times.

So walking around a braid-like node never adds more than
[`_CRW_M_BRAID`](@ref) letters: beyond that the node identifies everything, and
the appended part is the longest element of the dihedral subgroup — written
canonically as the alternating word starting with the SMALLER colour, so that
both directions give the same string (`121`, never `212`).

**Why it matters**: the 10-armed node with two opposite dots and the diagram
"trivalent + 6-armed + trivalent" are the SAME diagram, yet under the plain level
BFS the first one carries the unreduced word `1212`. With the saturation it reads
`121`, and every one of its region words is reduced.

Read by [`circular_region_words_2k`](@ref) and by
[`circular_unreduced_transition`](@ref). Set to `false` for the plain level BFS.
"""
const CIRCULAR_REGIONWORD_2K = Ref(true)

"""
    CIRCULAR_REGIONWORD_DOT_POLICY

Let the dot policy ([`circular_dot_may_pass`](@ref)) gate the `:dot`
transitions of [`circular_unreduced_transition`](@ref), default **`true`**.

A dot sitting on a node it may not pass through cannot start a rule, so it must
not raise an alarm either: with more arms nothing happens, and that counts as
reduced. Without this the 10-armed diagram reports a `dot`/`:direct` transition at
its 1-dot although no rule can act there.
"""
const CIRCULAR_REGIONWORD_DOT_POLICY = CIRCULAR_DOT_POLICY   # ONE switch for both rule and criterion

"The `(node, slot)` pairs of `g` that carry a dot (a 1-armed neighbour)."
function _crw_dot_slots(g::CircularGraph)
    out = Set{Tuple{Int,Int}}()
    for e in g.edges
        (e.a isa NodePort && e.b isa NodePort) || continue
        arm_count(g.nodes[e.b.node]) == 1 && push!(out, (e.a.node, e.a.slot))
        arm_count(g.nodes[e.a.node]) == 1 && push!(out, (e.b.node, e.b.slot))
    end
    return out
end

"""
    _crw_walls(g, v, sec, dots) -> Vector{Tuple{Int,Int}}

The WALLS around node `v`, cyclically: `(colour, region behind that wall)`.
A wall is a wired, non-dot arm — dot arms separate no regions
(circular/CircularRegion.jl), so the two sectors next to a dot are one region and the
dot arm is skipped. `sec` is `circular_node_sector_regions(g)`.
"""
function _crw_walls(g::CircularGraph, v::Int, sec, dots::Set{Tuple{Int,Int}})
    nd = g.nodes[v]
    return [(arm_colour(nd, s), sec[v][s]) for s in 1:arm_count(nd)
            if sec[v][s] != 0 && !((v, s) in dots)]
end

"""
    _crw_sector_words(walls, R, w, cap, full) -> Vector{Tuple{Int,Vector{Int}}}

All sectors of ONE node, seen from the sector `R` that carries the word `w`:
walk the wall ring in both directions, append the crossed colours, and cap the
appended part at `cap` letters — beyond that it is `full`, the node's saturated
word ([`_crw_node_cap`](@ref)).
"""
function _crw_sector_words(walls::Vector{Tuple{Int,Int}}, R::Int, w::Vector{Int},
                           cap::Int, full::Vector{Int})
    p = length(walls)
    p == 0 && return Tuple{Int,Vector{Int}}[]
    regs = [walls[i][2] for i in 1:p]
    out = Dict{Int,Vector{Int}}()
    for i0 in (i for i in 1:p if regs[i] == R), dir in (+1, -1)
        cols = Int[]
        for step in 1:(p - 1)
            wi = dir == +1 ? mod1(i0 + step, p) : mod1(i0 - step + 1, p)
            push!(cols, walls[wi][1])
            R2 = dir == +1 ? regs[mod1(i0 + step, p)] : regs[mod1(i0 - step, p)]
            R2 == R && continue
            cand = vcat(w, length(cols) < cap ? cols : full)
            (!haskey(out, R2) || length(cand) < length(out[R2])) && (out[R2] = cand)
        end
    end
    return [(R2, wd) for (R2, wd) in out]
end

"""
    circular_region_words_2k(m::CircularMorphismGraph) -> Vector{Union{Nothing,Vector{Int}}}

The region words of the SATURATED reading ([`CIRCULAR_REGIONWORD_2K`](@ref)) —
an OWN function next to [`circular_region_words`](@ref), like the m-valent one.

Dijkstra by WORD LENGTH: a region is reached either by an ordinary edge
crossing (append the colour, cost 1) or through the sector ring of a node,
where the appended part is capped at that node's own bound
([`_crw_node_cap`](@ref) — 3 at a `{1,2}`/`{2,3}` node, 2 at a `{1,3}` one,
1 at a single-colour node). The shortest word wins; ties are broken by region
number, so the result is deterministic.

⚠ Returns the words ALONE — no `parent_edge`. A saturated step is not a single
edge crossing, so there is no tree edge to record; whoever needs the cut curve
keeps reading [`circular_region_words`](@ref) (as `circular_path_edges` does).
"""
function circular_region_words_2k(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    nreg = region_count(g)
    words = Vector{Union{Nothing, Vector{Int}}}(nothing, nreg)
    (n == 0 || _circular_left_mark(m) === nothing) && return words
    regs = regions(g)
    gap = mod1(m.cut1, n)
    start = findfirst(R -> gap in R.gaps, regs)
    start === nothing && return words

    adj, _ = circular_region_adjacency(g)
    sec  = circular_node_sector_regions(g)
    dots = _crw_dot_slots(g)
    # EVERY node saturates, each with its own cap (`_crw_node_cap`) — dots
    # (1 arm) have no sectors, so they drop out by themselves.
    nodeinfo = [(_crw_walls(g, v, sec, dots), _crw_node_cap(nd)...)
              for (v, nd) in enumerate(g.nodes) if arm_count(nd) >= 3]

    words[start] = Int[]
    done = Set{Int}()
    cand = Dict{Int, Vector{Int}}()        # survives the rounds — a Dijkstra queue
    while true
        for R in 1:nreg
            (words[R] === nothing || R in done) && continue
            push!(done, R)
            w = words[R]::Vector{Int}
            for (R2, ei) in sort(adj[R]; by = x -> x[2])
                words[R2] === nothing || continue
                c = vcat(w, g.edges[ei].colour)
                (!haskey(cand, R2) || length(c) < length(cand[R2])) && (cand[R2] = c)
            end
            for (walls, cap, full) in nodeinfo, (R2, wd) in _crw_sector_words(walls, R, w, cap, full)
                words[R2] === nothing || continue
                (!haskey(cand, R2) || length(wd) < length(cand[R2])) && (cand[R2] = wd)
            end
        end
        for R in collect(keys(cand))
            words[R] === nothing || delete!(cand, R)
        end
        isempty(cand) && break
        L = minimum(length(v) for v in values(cand))
        for R2 in sort(collect(keys(cand)))
            if length(cand[R2]) == L
                words[R2] = cand[R2]
                delete!(cand, R2)
            end
        end
    end
    return words
end


"""
    circular_unreduced_transition(m::CircularMorphismGraph) -> Union{Nothing, NamedTuple}

The SHORTEST transition from a region with reduced region word `w_i` to a
non-reduced continuation `w_i · s` (file head), on the words of
[`circular_region_words`](@ref) (the level-BFS tree). Returns
`(region, word, colour, kind, edge, prev_edge, to, case)`:

* `region`    — the region i (reduced word `word`),
* `colour`    — the letter `s` (a right descent of `word`),
* `kind`      — `:edge` (boundary `s`-edge, one level further, to region `to`)
                or `:dot` (`s`-dot node `to` inside region i; `edge` is its cap),
* `prev_edge` — the tree edge INTO region i (0 at the start region). For a
                `:direct` `:edge` transition this is the parallel partner: the
                pair of the 2parallel rule is exactly `(prev_edge, edge)`,
* `case`      — `:direct`, `:dihedral` or `:long` (see file head).

`nothing` if every continuation stays reduced, or without a mark. Shortest =
smallest word length, then region number, then edges before dots, then edge
index (deterministic).
"""
function circular_unreduced_transition(m::CircularMorphismGraph)
    g = m.graph
    nreg = region_count(g)
    nreg == 0 && return nothing
    # The m-valent reading, behind its switch.
    # Cost-0 (fold) crossings keep the word length, so the `length + 1` check
    # below skips them by itself — no transition is read across a fold leg. A
    # node jump is no edge, so it is never a transition either; the transitions
    # are read on the edges as before, only the WORDS come from the other
    # reading.
    # Reading order: the m-valent one wins if it is switched on, otherwise the
    # 2k-saturated one (the DEFAULT), otherwise the plain level BFS. `parent`
    # always comes from the level BFS — a saturated step is no single edge
    # crossing, so it has no tree edge (see `circular_region_words_2k`), and
    # `prev_edge` keeps its meaning.
    words, parent = if CIRCULAR_REGIONWORD_MVALENT[]
        circular_region_mvalent_words(m)
    elseif CIRCULAR_REGIONWORD_2K[]
        (circular_region_words_2k(m), circular_region_words(m)[2])
    else
        circular_region_words(m)
    end
    adj, _ = circular_region_adjacency(g)

    best = nothing
    bestkey = nothing
    function consider(T)
        k = (length(T.word), T.region, T.kind === :edge ? 0 : 1, T.edge)
        if best === nothing || k < bestkey
            best = T; bestkey = k
        end
    end

    for R in 1:nreg
        w = words[R]
        w === nothing && continue
        is_reduced(w) || continue
        for (R2, ei) in adj[R]
            w2 = words[R2]
            (w2 !== nothing && length(w2) == length(w) + 1) || continue
            s = g.edges[ei].colour
            is_reduced(vcat(w, s)) && continue
            consider((region = R, word = w, colour = s, kind = :edge, edge = ei,
                      prev_edge = parent[R], to = R2, case = _crw_case(w, s)))
        end
    end
    for (v, nd) in enumerate(g.nodes)
        arm_count(nd) == 1 || continue
        R = circular_region_of_dot(g, v)
        w = words[R]
        w === nothing && continue
        is_reduced(w) || continue
        s = arm_colour(nd, 1)
        is_reduced(vcat(w, s)) && continue
        at = _circular_edges_at_node(g, v)
        isempty(at) && continue
        # The DOT POLICY (`circular_dot_may_pass`): a dot that may not be pushed
        # through its node is no transition — with
        # more arms nothing happens, and that also counts as reduced.
        if CIRCULAR_REGIONWORD_DOT_POLICY[]
            e = g.edges[at[1][1]]
            other = e.a isa NodePort && e.a.node == v ? e.b : e.a
            other isa NodePort && !circular_dots_reducible(g, other.node) && continue
        end
        consider((region = R, word = w, colour = s, kind = :dot, edge = at[1][1],
                  prev_edge = parent[R], to = v, case = _crw_case(w, s)))
    end
    return best
end

"The case of a transition (file head): `:direct`, `:dihedral` or `:long`."
function _crw_case(w::Vector{Int}, s::Int)
    isempty(w) && return :long
    w[end] == s && return :direct
    for site in braid_move_sites(w)
        v = apply_braid_move(w, site...)
        (!isempty(v) && v[end] == s) && return :dihedral
    end
    return :long
end

"""
    circular_regionword_dihedral_step(fdm::CircularDecoratedMorphism)
        -> Union{Nothing, NamedTuple}

Case 2 (the `121` + 2-dot pattern): if the shortest transition is a `:dot`
transition of case `:dihedral`, apply the DOT-SLIDE rule
([`circular_dot_slide_step`](@ref)) — with the transition as the GATE.

⚠ Forcing `apply_circular_d4` at the transition edge does NOT work: the dot's
region carries no edge of the dot colour at all, and the D4 lives on the braid²
window AFTER the slide surgery. The move that performs exactly that (surgery + D4
+ rest window) is the dot-slide rule (circular/rules/CircularDotSlide.jl, switch
`CIRCULAR_DOT_SLIDE_ENABLED`, default OFF, no termination measure). So case 2 is a
criterion-gated dot slide, not a new surgery.

Returns `(combo, dot, transition)` or `nothing` (no `:dot`/`:dihedral`
transition, or the slide does not fire). A separate entry point: it is
not part of `reduce_to_circular_leave` while the dot-slide termination is
unwritten.
"""
function circular_regionword_dihedral_step(fdm::CircularDecoratedMorphism)
    T = circular_unreduced_transition(fdm.m)
    (T === nothing || T.case !== :dihedral || T.kind !== :dot) && return nothing
    isone(fdm.region_labels[T.region]) || return nothing
    r = circular_dot_slide_step(fdm)
    r === nothing && return nothing
    return (combo = r, dot = T.to, transition = T)
end

"""
    circular_path_edges(m::CircularMorphismGraph, R::Int)
        -> Union{Nothing, NamedTuple}

The region word of `R` AS A LIST OF EDGES — the tree path from the marked region
to `R`, read from the mark outwards.

The path is the `parent_edge` chain of [`circular_region_words`](@ref) (the level
BFS), walked upwards from `R` and reversed. Returns
`(edges, word, regions)`:

* `edges`   — the crossed edges `e_1..e_k`, mark first,
* `word`    — their colours, which IS `circular_region_words(m)[1][R]`,
* `regions` — the visited regions, `regions[1]` = the marked one, `regions[end] = R`.

`nothing` if `R` is unreachable (no mark, other component). For the marked
region itself the answer is the empty path.

These `k` edges are the CUT CURVE: cut each of them once and splice a morphism
`w → w` in between. Whether that is always planar is not established here — this
function supplies the curve, it does not judge it.
"""
function circular_path_edges(m::CircularMorphismGraph, R::Int)
    g = m.graph
    (1 <= R <= region_count(g)) || return nothing
    words, parent = circular_region_words(m)
    words[R] === nothing && return nothing
    adj, _ = circular_region_adjacency(g)

    edges = Int[]
    regs = Int[R]
    cur = R
    while parent[cur] != 0
        ei = parent[cur]
        push!(edges, ei)
        up = nothing
        for (R2, e2) in adj[cur]
            e2 == ei || continue
            w2 = words[R2]
            (w2 !== nothing && length(w2) == length(words[cur]::Vector{Int}) - 1) || continue
            up = R2
            break
        end
        up === nothing && return nothing        # broken tree — never seen, not guessed away
        push!(regs, up)
        cur = up
    end
    reverse!(edges)
    reverse!(regs)
    return (edges = edges, word = [g.edges[e].colour for e in edges], regions = regs)
end
