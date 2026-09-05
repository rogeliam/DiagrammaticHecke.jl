# circular/CircularRegion.jl — region adjacency and distances for CircularGraph.
#
# A REGION IS EXACTLY ONE CELL. The dart -> region assignment
# comes straight out of planarity: the tracer (diagram/Faces.jl) already runs
# the boundary arcs as genuine darts, so its face walk separates at every
# strand AND at every boundary arc. `dart_region[d] = cell_of_face[face_of_dart[d]]`,
# done.
#
# Rule for edges (unchanged): an edge separates two regions iff it is not
# incident to a dot (a 1-armed node).

"""
    _circular_dart_region(g::CircularGraph) -> Tuple{Dict{Int,Int}, Any}

Assigns each **inner** dart (a directed edge) the region on its left side.
Non-dot edges therefore separate two distinct regions.

Returns:
- `dart_region`: dart id -> region id
- `t`: the full tracer result (`_trace_generic`)
"""
function _circular_dart_region(g::CircularGraph)
    t = _trace_generic(g.word, g.nodes, g.edges, _circular_degree)

    # Region == cell, and the assignment falls out of the tracer
    # directly: it runs boundary arcs as genuine darts, so its face walk already
    # separates exactly where a region boundary should be — at every strand and
    # at every leaf-to-leaf connection. No gap labels, no backward search, no
    # special case for borderless regions.
    dart_region = Dict{Int, Int}()
    for (fi, cyc) in enumerate(t.cycles)
        fi == t.outer_face && continue
        c = t.cell_of_face[fi]
        for d in cyc
            dart_region[d] = c
        end
    end

    return dart_region, t
end

"""
    circular_node_sector_regions(g::CircularGraph) -> Vector{Vector{Int}}

For each node `v` and each wired slot `s`, the **region in the sector AFTER
`s`** (clockwise) — the area between slot `s` and the next wired slot. Unwired
slots get `0`.

The region analogue of `Cells.sector_cell` (diagram/Faces.jl), built by exactly
the same formula: *sector after slot `s` = the face left of the dart at the
NEXT slot* — slots are ordered clockwise (Faces.jl header), and the face walk
leaves a node it entered via slot `s` via slot `s+1`. Here `dart_region[…]`
from [`_circular_dart_region`](@ref) stands in for `cell_of_face[face_of_dart[…]]`.

Used by the D4 node variant (`find_circular_d4_node_match`, circular/CircularDecoratedRules.jl):
it needs to know which sector of a node points into the dot's region — that is
where the new arm gets inserted.
"""
function circular_node_sector_regions(g::CircularGraph)
    dart_region, t = _circular_dart_region(g)
    out = Vector{Vector{Int}}(undef, length(g.nodes))
    for (ni, nd) in enumerate(g.nodes)
        d = arm_count(nd)
        sr = zeros(Int, d)
        wired = [(s, t.slot_dart[(ni, s)]) for s in 1:d if haskey(t.slot_dart, (ni, s))]
        for (i, (s, _)) in enumerate(wired)
            d_next = wired[mod1(i + 1, length(wired))][2]
            sr[s] = get(dart_region, d_next, 0)
        end
        out[ni] = sr
    end
    return out
end

"""
    circular_region_adjacency(g::CircularGraph)
        -> Tuple{Vector{Vector{Tuple{Int,Int}}}, Dict{Tuple{Int,Int}, Int}}

Region adjacency: for each region `R`, `adj[R]` is a list of
`(neighbour region, edge index)` pairs sorted by edge index. The edge index
refers to `g.edges`; it is never a dot edge. `edge_between[(r1,r2)]` (with
`r1 ≤ r2`) returns the separating edge.
"""
function circular_region_adjacency(g::CircularGraph)
    dart_region, t = _circular_dart_region(g)
    nregs = region_count(g)
    adj = [Tuple{Int,Int}[] for _ in 1:nregs]
    edge_between = Dict{Tuple{Int,Int}, Int}()

    for (ei, e) in enumerate(g.edges)
        # Dot edges never separate regions.
        ((e.a isa NodePort && arm_count(g.nodes[e.a.node]) == 1) ||
         (e.b isa NodePort && arm_count(g.nodes[e.b.node]) == 1)) && continue
        haskey(t.port_dart, e.a) || continue

        d  = t.port_dart[e.a]      # e.a -> e.b
        dr = t.darts[d].rev        # e.b -> e.a
        R1 = get(dart_region, d, nothing)
        R2 = get(dart_region, dr, nothing)
        (R1 === nothing || R2 === nothing || R1 == R2) && continue

        pair = (min(R1, R2), max(R1, R2))
        haskey(edge_between, pair) && continue
        edge_between[pair] = ei
        push!(adj[R1], (R2, ei))
        push!(adj[R2], (R1, ei))
    end
    for lst in adj
        sort!(lst, by = x -> x[2])
    end
    return adj, edge_between
end

"""
    node_class(nd::CircularNode) -> Symbol

Classifies a node by the colours of its arms: `:small` (≤ 2 arms, i.e. a dot or
a pass-through edge), `:gen2` (all arms colour 2), `:gen13` (all arms in
`{1,3}`), `:braid` (mixed, i.e. a 12-/23-braid).

Used by [`circular_region_distances`](@ref): only at `:gen2`/`:gen13` nodes do the
regions around the node count as additionally adjacent.
"""
function node_class(nd::CircularNode)
    arm_count(nd) <= 2 && return :small
    cs = sort(unique(nd.arms))
    cs == [2]                && return :gen2
    all(c -> c in (1, 3), cs) && return :gen13
    return :braid
end

"""
    node_cap(cls::Symbol) -> Union{Int, Nothing}

Cap on the extra adjacency of a node class (see [`node_class`](@ref)): `1` at
general-2, `2` at general-13 nodes, `nothing` = no extra edges at all
(`:braid`, `:small`).

The cap is the number of colours that occur at such a node at all — going
around the node never has to cost more than once per colour.
"""
node_cap(cls::Symbol) = cls === :gen2 ? 1 : cls === :gen13 ? 2 : nothing

"""
    circular_region_distances(m::CircularMorphismGraph) -> Vector{Int}

Distance of each region from the black marking (the gap between the bottom and
top start): a touching region is 0, and so on. Unreachable regions get `-1`;
with no marking, all regions come back `-1`.

**THE BOUNDARY CIRCLE IS A WALL.**
The search never leaves the disc: it walks diagram edges (and the node
neighbourhoods below), never a boundary arc, so there is no passage
outer ↔ inner. The outer region is therefore NOT part of this vector — it has a
fixed distance of its own, see [`circular_outer_region_distance`](@ref). Boundary
regions of a component that the marking cannot reach through edges keep `-1`;
that is the wall, not a defect.

**DOT EDGES COST 0** (as elsewhere: crossing a non-dot edge raises the
distance by 1). A dot edge — an edge with a 1-armed node (a dot) at one end —
never separates two regions in the first place: the face walk enters the dot
and leaves along the reverse dart, so both sides lie in the same face and
hence in the same region ([`circular_region_adjacency`](@ref) skips them
explicitly). Crossing a dot edge therefore costs 0 *by construction*, and the
rule needs no weight of its own. A leaf–leaf strand is NOT a dot edge — its
ends are boundary leaves, not dots — so the minimal figure `1111` keeps its
distances `[0, 1, 2]`.

**Edges are not the only thing that counts.** Adjacencies are:

- **(a) edges**, as before, weight 1 (`circular_region_adjacency`);
- **(b) nodes**: at a `:gen2`/`:gen13` node ([`node_class`](@ref)) *all* sector
  regions of the node count as adjacent, with weight = the number of
  **distinct colours** crossed going around (the minimum over the two arcs),
  capped by [`node_cap`](@ref). Braid nodes are unaffected — there the colours
  genuinely separate.

Because (b) introduces weights ≠ 1, the search is Bellman-Ford, not BFS (the
region graphs are tiny, `O(n²)` is fine).

**ATTENTION, TWO DISTANCES.** This is the distance for rules that ask how far a
region is from the marking across the diagram — D4 and its node variant, the dot
slide, P3, and the region labels in the renderer. The rule that MOVES a
polynomial uses the pure edge distance
[`circular_region_distances_edges_only`](@ref) instead: `circular_fusion_step`
walks one edge at a time, and the node extra-adjacency (b) squeezes distances
together until the `d − 1` neighbour such a step needs is gone. That function's
docstring has the worked case.
"""
function circular_region_distances(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    nreg = region_count(g)
    (n == 0 || _circular_left_mark(m) === nothing) && return fill(-1, nreg)
    regs = regions(g)
    gap = mod1(m.cut1, n)
    start = findfirst(R -> gap in R.gaps, regs)
    start === nothing && return fill(-1, nreg)

    # Edge weights, taking the minimum over all reasons offered for the same pair.
    W = Dict{Tuple{Int,Int}, Int}()
    function add!(a, b, w)
        a == b && return
        k = (min(a, b), max(a, b))
        W[k] = min(get(W, k, typemax(Int)), w)
        return
    end

    # (a) the usual edge adjacencies
    adj, _ = circular_region_adjacency(g)
    for R in 1:nreg, (R2, _) in adj[R]
        add!(R, R2, 1)
    end

    # (b) the node adjacencies at general-2/general-13 nodes
    sec = circular_node_sector_regions(g)
    for (v, nd) in enumerate(g.nodes)
        cap = node_cap(node_class(nd))
        cap === nothing && continue
        wired = [s for s in 1:arm_count(nd) if sec[v][s] != 0]
        k = length(wired)
        k < 2 && continue
        Rs   = [sec[v][s] for s in wired]
        cols = [arm_colour(nd, s) for s in wired]
        for i in 1:k, j in 1:k
            i == j && continue
            # colours crossed by path i->j clockwise resp. counter-clockwise.
            fwd = Int[]; t = i; while t != j; t = mod1(t + 1, k); push!(fwd, cols[t]); end
            bwd = Int[]; t = j; while t != i; t = mod1(t + 1, k); push!(bwd, cols[t]); end
            add!(Rs[i], Rs[j], min(length(unique(fwd)), length(unique(bwd)), cap))
        end
    end

    # (c) NO passage across the boundary circle: the outer region is a wall, not
    # a node here. Dot edges are already absent from
    # `adj` — they separate nothing, so crossing them costs 0 by construction.

    nb = [Tuple{Int,Int}[] for _ in 1:nreg]
    for ((a, b), w) in W
        push!(nb[a], (b, w))
        push!(nb[b], (a, w))
    end

    dist = fill(typemax(Int), nreg)
    dist[start] = 0
    for _ in 1:nreg, a in 1:nreg          # small graphs: Bellman-Ford is enough
        dist[a] == typemax(Int) && continue
        for (b, w) in nb[a]
            dist[b] = min(dist[b], dist[a] + w)
        end
    end
    return [d == typemax(Int) ? -1 : d for d in dist]
end

"""
    circular_region_tree_words(m::CircularMorphismGraph) -> Vector{Union{Nothing, Vector{Int}}}

The **tree_word label** of each region: alongside the node-weighted
distance [`circular_region_distances`](@ref), every region gets the WORD of colours
that a shortest path from the marking region crosses to reach it — "walking up
a tree". The start region (distance 0) carries the empty word; each edge
crossing appends its colour.

Same adjacencies as [`circular_region_distances`](@ref): (a) edge steps append
the ONE edge colour; (b) node steps at `:gen2`/`:gen13` nodes append the
DISTINCT colours of the cheaper arc (in crossing order, each colour once) — at
a general-2 node this is exactly `[2]`, the distance grows by at most 1
([`node_cap`](@ref)). So always
`length(tree_word[R]) == circular_region_distances(m)[R]`.

When several shortest paths exist, the FIRST one found wins (Bellman-Ford only
relaxes on a genuine improvement). The region word is the TREE word: every cell
gets its word from exactly ONE predecessor cell + edge, and the fixed
construction order pins the tree down deterministically. In words: distances "0,
all neighbours 1, their neighbours 2 and so on", words "empty, all neighbours
i/j/k, their neighbours length 2 and so on" — level-by-level BFS, one letter
per level; the node jumps of the weighted version can deviate from that.
Unreachable regions (distance -1) and the no-marking case return `nothing`.

CONJECTURE (test in `test/circular.jl`): if a region carries a
NON-REDUCED tree_word ([`is_reduced`](@ref)), the diagram is not in normal
form — on the circular-leaf end terms of `reduce_to_circular_leave`, all tree_words are
reduced.
"""
function _circular_region_steps(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    nreg = region_count(g)
    (n == 0 || _circular_left_mark(m) === nothing) && return nothing
    regs = regions(g)
    gap = mod1(m.cut1, n)
    start = findfirst(R -> gap in R.gaps, regs)
    start === nothing && return nothing

    # Directed step edges (a -> b) with the appended colours; the weight is
    # ALWAYS `length(letters)` — exactly the weights from `circular_region_distances`,
    # just carrying the justification along.
    steps = [Tuple{Int, Vector{Int}}[] for _ in 1:nreg]

    # (a) edge adjacencies: one edge, one colour.
    adj, edge_between = circular_region_adjacency(g)
    for R in 1:nreg, (R2, _) in adj[R]
        e = edge_between[(min(R, R2), max(R, R2))]
        push!(steps[R], (R2, [g.edges[e].colour]))
    end

    # (b) node adjacencies at :gen2/:gen13 nodes (as in circular_region_distances,
    # but carrying the crossed colours as a word).
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
            push!(steps[Rs[i]], (Rs[j], best))
        end
    end

    return (steps = steps, start = start, nreg = nreg)
end

"""
    circular_region_tree_parents(m::CircularMorphismGraph)
        -> (words, parent, letters)

The region-word tree of [`circular_region_tree_words`](@ref) with the TREE
itself: `words[R]` is the word, `parent[R]` the region the winning step came
from (`0` at the root and where there is none), and `letters[R]` the letters
that step appended.

Why a parent REGION and not a parent EDGE: a step of this tree may cross a
`:gen2`/`:gen13` NODE instead of an edge, and then there is no edge to name. The
edges-only walk of [`circular_region_words`](@ref) can report a parent edge
because every one of its steps has one; this walk cannot, and a drawing that
insists on an edge shows a crossing that does not happen.
"""
function circular_region_tree_parents(m::CircularMorphismGraph)
    nreg = region_count(m.graph)
    words = Vector{Union{Nothing, Vector{Int}}}(nothing, nreg)
    parent = zeros(Int, nreg)
    letters_of = Vector{Union{Nothing, Vector{Int}}}(nothing, nreg)
    st = _circular_region_steps(m)
    st === nothing && return (words = words, parent = parent, letters = letters_of)
    steps, start = st.steps, st.start

    dist = fill(typemax(Int), nreg)
    dist[start] = 0
    words[start] = Int[]
    for _ in 1:nreg, a in 1:nreg          # small graphs: Bellman-Ford is enough
        dist[a] == typemax(Int) && continue
        for (b, letters) in steps[a]
            nd2 = dist[a] + length(letters)
            if nd2 < dist[b]              # only a GENUINE improvement: first witness stays
                dist[b] = nd2
                words[b] = vcat(words[a]::Vector{Int}, letters)
                parent[b] = a
                letters_of[b] = letters
            end
        end
    end
    return (words = words, parent = parent, letters = letters_of)
end

function circular_region_tree_words(m::CircularMorphismGraph)
    return circular_region_tree_parents(m).words
end

"""
    circular_region_distances_edges_only(m::CircularMorphismGraph) -> Vector{Int}

Distance of each region from the black marking, like
[`circular_region_distances`](@ref), but **ONLY genuine edges**
(`circular_region_adjacency`, weight 1) — without the node extra-adjacency at
`:gen2`/`:gen13` nodes. Plain BFS instead of Bellman-Ford.

**Same wall, same dot-edge rule** as [`circular_region_distances`](@ref): the outer
region is not a node here either, and dot edges cost 0 because they separate
nothing.

**What it's for — the FUSION**, the one rule that moves polynomials: it pushes
them through edges step by step, and that needs the distance that simply decreases by 1 at each step.
This version counts exactly the steps `circular_fusion_step` can take. Under the
node-weighted distance the `:gen13` star node has no `d − 1` neighbour at all
(`[0,1,1,2,2,2,1,1]` vs. `[0,1,2,3,4,3,2,1]`), so fusion cannot fire there.

**Every other caller wants the node-weighted [`circular_region_distances`](@ref)**
— D4, the dot slide, P3 and the renderer all ask about reachability across the
diagram, not about single edge steps.
"""
function circular_region_distances_edges_only(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    nreg = region_count(g)
    (n == 0 || _circular_left_mark(m) === nothing) && return fill(-1, nreg)
    regs = regions(g)
    gap = mod1(m.cut1, n)
    start = findfirst(R -> gap in R.gaps, regs)
    start === nothing && return fill(-1, nreg)

    # No node for the outer region: the boundary circle is a wall.
    adj, _ = circular_region_adjacency(g)

    dist = fill(-1, nreg)
    dist[start] = 0
    queue = [start]
    while !isempty(queue)
        R = popfirst!(queue)
        for (R2, _) in adj[R]
            dist[R2] == -1 || continue
            dist[R2] = dist[R] + 1
            push!(queue, R2)
        end
    end
    return dist
end

"""
    circular_boundary_edge(g::CircularGraph, r1::Int, r2::Int) -> Union{Int, Nothing}

The index in `g.edges` of the edge separating regions `r1` and `r2`, or
`nothing` if they are not adjacent. The edge is never a dot edge.
"""
function circular_boundary_edge(g::CircularGraph, r1::Int, r2::Int)
    _, edge_between = circular_region_adjacency(g)
    return get(edge_between, (min(r1, r2), max(r1, r2)), nothing)
end

"""
    circular_region_of_dot(g::CircularGraph, dot_node::Int) -> Int

The region the dot `dot_node` lies in. Dots are 1-armed and transparent; they
belong to the region their capping edge encloses.
"""
function circular_region_of_dot(g::CircularGraph, dot_node::Int)
    1 <= dot_node <= length(g.nodes) || error(
        "circular_region_of_dot: invalid node index $dot_node")
    arm_count(g.nodes[dot_node]) == 1 || error(
        "circular_region_of_dot: node $dot_node is not a dot")

    cap_ei = findfirst(e -> (e.a isa NodePort && e.a.node == dot_node) ||
                            (e.b isa NodePort && e.b.node == dot_node),
                       g.edges)
    cap_ei === nothing && error(
        "circular_region_of_dot: dot $dot_node has no edge")
    e = g.edges[cap_ei]
    p = (e.a isa NodePort && e.a.node == dot_node) ? e.b : e.a

    if p isa Leaf
        k = p.k
        n = length(g.word)
        regs = regions(g)
        # both gaps at leaf k belong to the same region
        R = findfirst(r -> k in r.gaps || mod1(k - 1, n) in r.gaps, regs)
        R === nothing && error(
            "circular_region_of_dot: leaf $k not found in any region")
        return R
    elseif p isa NodePort
        dart_region, t = _circular_dart_region(g)
        d = t.port_dart[p]   # dart p -> dot
        return dart_region[d]
    else
        error("circular_region_of_dot: unexpected port type $(typeof(p))")
    end
end

# ---------------------------------------------------------------------------
# THE DISTANCE-WORD CRITERION
#
#   "for dots, look at the distance word; if it can end on the dot's colour,
#    it can be reduced further. That only works if it ends with the
#    same colour, or a short dihedral."
#
# The word comes from `circular_region_tree_words` (distance 0 -> outward). A dot's
# DISTANCE WORD is that word read BACKWARDS — from the dot's cell to the
# marking. If braid moves can bring it to end on the dot's colour, the dot is
# reducible and the diagram is not in normal form.
# ---------------------------------------------------------------------------

"""
    circular_braid_class(w::Vector{Int}; limit = 20_000) -> Set{Vector{Int}}

All reduced words for the same A₃ element as `w` — i.e. everything reachable
from `w` by braid and commutation moves. Delegates to
[`reduced_words`](@ref) (`src/morphism/BraidMoves.jl`), the checked machinery.

A **non-reduced** `w` has no braid class in this sense (`reduced_words` demands
a reduced word) and gets `Set([w])` back — the singleton, no witnesses.
Non-reduced region words are the job of the SEPARATE region-word criterion
(E2), not of the dot criterion.

`limit` caps the returned set; A₃ classes are far smaller than the default.
"""
function circular_braid_class(w::Vector{Int}; limit::Int = 20_000)
    isempty(w) && return Set([w])
    is_reduced(w) || return Set([w])
    out = Set{Vector{Int}}()
    for v in reduced_words(w)
        length(out) < limit || break
        push!(out, v)
    end
    return out
end

"""
    circular_distance_word(m::CircularMorphismGraph, dot_node::Int) -> Union{Nothing, Vector{Int}}

The distance word of the dot node `dot_node`: the region word of its region
([`circular_region_words`](@ref)), read BACKWARDS — i.e. from the dot's cell
toward the marking. `nothing` if the region is unreachable or there is no
marking.

The path is NOT searched backwards: the word arises in the distance
computation from the marking outward and is only reversed
here.

Word source: the EDGE BFS `circular_region_words`, not the weighted
`circular_region_tree_words` (one notion of word in the whole
package). The two differ in 235/2130 end-term regions, but as GROUP
elements the two words agree, so the choice is free and falls on the one the
rules use.
"""
function circular_distance_word(m::CircularMorphismGraph, dot_node::Int)
    g = m.graph
    arm_count(g.nodes[dot_node]) == 1 || error(
        "circular_distance_word: node $dot_node is not a dot (arm_count != 1)")
    R = circular_region_of_dot(g, dot_node)
    w = circular_region_words(m)[1][R]
    w === nothing && return nothing
    return reverse(w)
end

"""
    circular_dot_reducible(m::CircularMorphismGraph, dot_node::Int; limit = 20_000) -> Bool

The distance-word criterion: `true` if the dot's distance word can be brought
to end on the **dot's colour** by braid moves — then the dot can be reduced
further, and the diagram is NOT in normal form.

An empty distance word (the dot sits in the marking region) and an
unreachable region give `false`.

Example: the 2-dot on the `1212` diagram sits in a region with
tree_word `121`; the braid class is `{121, 212}`, and `212` ends on 2 = the
dot's colour ⇒ `true`.
"""
function circular_dot_reducible(m::CircularMorphismGraph, dot_node::Int; limit::Int = 20_000)
    w = circular_distance_word(m, dot_node)
    (w === nothing || isempty(w)) && return false
    s = arm_colour(m.graph.nodes[dot_node], 1)
    return any(v -> !isempty(v) && v[end] == s, circular_braid_class(w; limit = limit))
end

"""
    _circular_dot_on_leaf(g::CircularGraph, v::Int) -> Bool

`true` if the dot `v` hangs on a BOUNDARY LEAF, i.e. its single edge ends on a
`Leaf` rather than on another node.
"""
function _circular_dot_on_leaf(g::CircularGraph, v::Int)
    for e in g.edges
        for (p, q) in ((e.a, e.b), (e.b, e.a))
            (p isa NodePort && p.node == v) || continue
            return q isa Leaf
        end
    end
    return false
end

"""
    circular_non_normal_dots(m::CircularMorphismGraph; limit = 20_000) -> Vector{Int}

The dots of `m` that witness it is NOT in normal form: those on a BOUNDARY LEAF
that also satisfy [`circular_dot_reducible`](@ref). Empty means no dot violates
the criterion.

ONLY LEAF DOTS COUNT. A dot sitting on another node — a braid node in
particular — is a normal form: the distance word measures a path through
regions to the dot, and a dot absorbed into a node is not something the
braid-class argument can move. Applying the criterion there reports diagrams
that are already done, which is what it used to do on the Zamolodchikov terms:
every dot flagged there hangs on the 12-armed braid node.
"""
function circular_non_normal_dots(m::CircularMorphismGraph; limit::Int = 20_000)
    out = Int[]
    for (v, nd) in enumerate(m.graph.nodes)
        arm_count(nd) == 1 || continue
        _circular_dot_on_leaf(m.graph, v) || continue
        circular_dot_reducible(m, v; limit = limit) && push!(out, v)
    end
    return out
end
