# circular/rules/CircularGen12Expand.jl — C13⁻¹: EXPAND a gbraid back into a cluster
# of two braids and two trivalent nodes.
#
# THE PREIMAGES ARE HARDCODED. Finding them by breadth-first search ("expand
# R9⁻¹, fully reduce, repeat", fixed point at depth 5) takes about 17 minutes,
# so the result is recorded here as a literal.
#
# THERE ARE FOUR, NOT TWELVE. The search returns twelve results, but under
# `circular_canonical_key` they collapse into **four classes of three**:
#
#     [1,4,7]   [2,9,11]   [3,8,10]   [5,6,12]
#
# Representatives used here: **1, 2, 3, 5** (in this order, `variant = 1…4`).
# Coarser still: up to leaf rotation there are **2** (trivalent nodes of one
# color or the other), up to rotation AND color swap there is **1**.
#
# ⚠ Twelve and four count different things: twelve `WordGraph` wirings (slots
# fixed and numbered) against four `CircularGraph`s (a node only carries its
# cyclic arm sequence) — exactly the distinction the circular type draws.
#
# NOT in `CIRCULAR_RULES`: this is an EXPANSION rule, it increases `circular_weight`
# (an 8-armed braid-like node becomes two 6-armed ones, component 1 rises
# 8 ↦ 12). It is a hand tool for notebooks and for searching upward — like
# `expand_circular_merge` in CircularInverseRules.jl.
#
# LITERAL FORMAT, one triple per preimage:
#   nodes  arm sequences of the four nodes (colors 1/2; the real colors come
#          from the gbraid, see below)
#   inner  (color, node a, slot a, node b, slot b) — the five inner edges
#   outer  (leaf, node, slot) — the eight outer ports, leaf `i` carries color
#          `i` odd ? 1 : 2 (boundary word `12121212`)
#
# ⚠ SCOPE: only the 8-armed case is covered by a literal. The 10/12/14-armed
# cases are handled by the general construction `gbraid_tower_*` below.

"""
    GBRAID_PREIMAGES

The four preimages of an 8-armed gbraid under C13 (`gen12_merge`), stored
as cluster literals `(nodes, inner, outer)`. See the file header for why there
are four and not twelve.
"""
const GBRAID_PREIMAGES = (
    # variant 1  (= preimage 1 of the search; class [1,4,7]) — trivalent nodes color 1
    ([[2,1,2,1,2,1],[1,1,1],[1,1,1],[1,2,1,2,1,2]],
     [(1,1,4,2,1),(1,1,6,3,3),(1,2,3,4,3),(1,3,1,4,1),(2,1,5,4,2)],
     [(1,2,2),(2,1,3),(3,1,2),(4,1,1),(5,3,2),(6,4,6),(7,4,5),(8,4,4)]),
    # variant 2  (= preimage 2; class [2,9,11]) — trivalent nodes color 2
    ([[1,2,1,2,1,2],[2,2,2],[1,2,1,2,1,2],[2,2,2]],
     [(1,1,5,3,1),(2,1,4,2,1),(2,1,6,4,1),(2,2,3,3,2),(2,3,6,4,2)],
     [(1,3,3),(2,2,2),(3,1,3),(4,1,2),(5,1,1),(6,4,3),(7,3,5),(8,3,4)]),
    # variant 3  (= preimage 3; class [3,8,10]) — trivalent nodes color 2, different placement
    ([[2,2,2],[2,1,2,1,2,1],[1,2,1,2,1,2],[2,2,2]],
     [(1,2,2,3,5),(2,1,1,3,4),(2,1,3,2,3),(2,2,1,4,3),(2,3,6,4,2)],
     [(1,2,6),(2,2,5),(3,2,4),(4,1,2),(5,3,3),(6,3,2),(7,3,1),(8,4,1)]),
    # variant 4  (= preimage 5; class [5,6,12]) — trivalent nodes color 1, different placement
    ([[1,1,1],[1,2,1,2,1,2],[2,1,2,1,2,1],[1,1,1]],
     [(1,1,1,3,4),(1,1,3,2,3),(1,2,1,4,3),(1,3,6,4,2),(2,2,2,3,5)],
     [(1,3,2),(2,3,1),(3,4,1),(4,2,6),(5,2,5),(6,2,4),(7,1,2),(8,3,3)]),
)

"""
    expand_gbraid(g::CircularGraph, v::Int, variant::Int; rot::Int = 0)
        -> Union{Nothing, CircularGraph}

**C13⁻¹.** Replaces the 8-armed gbraid node `v` of `g` with one of the
four preimage clusters (`variant ∈ 1:4`, see `GBRAID_PREIMAGES`) and
returns the resulting graph. `rot ∈ 0:7` rotates which arm of the gbraid
corresponds to the cluster's first leaf — for a fixed `variant`, the eight
values of `rot` run through the leaf rotations of the same preimage.

Colors are **not** guessed: the cluster is recorded in colors `1/2`, and
`1 ↦ arms(v)[rot+1]`, `2 ↦ arms(v)[rot+2]` translates it to the actual
gbraid's color pair (`{1,2}` as well as `{2,3}`). The boundary word and
`circular_degree` of `g` stay unchanged; `_fr_gen12_merge` maps the result back
onto `g` (round trip, checked in the test over all 4 × 8 combinations).

Returns `nothing` if `v` is not an 8-armed gbraid.

**Not in `CIRCULAR_RULES`** — this rule expands and increases `circular_weight`.
"""
function expand_gbraid(g::CircularGraph, v::Int, variant::Int; rot::Int = 0)
    (1 <= variant <= length(GBRAID_PREIMAGES)) ||
        throw(ArgumentError("expand_gbraid: variant must be in 1:$(length(GBRAID_PREIMAGES)), was $variant"))
    nd = g.nodes[v]
    (nd.kind === :braid && arm_count(nd) == 8) || return nothing
    return _circular_splice_cluster(g, v, GBRAID_PREIMAGES[variant]; rot = rot)
end

"""
    _circular_splice_cluster(g::CircularGraph, v::Int, lit; rot::Int = 0)
        -> Union{Nothing, CircularGraph}

Replaces node `v` (arm count `n`) with the cluster literal
`lit = (nodes, inner, outer)`, which has `n` outer leaves. Factored out of
`expand_gbraid` so that the 8-armed literal from `GBRAID_PREIMAGES` and the
constructed tower clusters (`gbraid_tower_literal`) go through the SAME
splice.

`nothing` if the color check fails.
"""
function _circular_splice_cluster(g::CircularGraph, v::Int, lit; rot::Int = 0)
    nd = g.nodes[v]
    n = arm_count(nd)
    a = collect(arms(nd))
    (nds, inner, outer) = lit
    length(outer) == n || return nothing

    # Leaf `i` of the cluster belongs at slot `mod1(2 - i + rot, n)` of the
    # node. ⚠ BACKWARDS, and that is not arbitrary: the house convention is
    # "leaves run backwards in slot order around the leaf ring"
    # (`check_wiring`, `:leaf_ring`). Replacing a node with a cluster walks
    # the boundary of the resulting hole in the opposite direction from the
    # diagram boundary the cluster is written against — the reversal is exactly
    # that switch. With `mod1(i + rot, n)` (forward) the result is NOT planar
    # (`euler = -4`, three `:leaf_ring` violations).
    slot_of_leaf = Dict(i => mod1(2 - i + rot, n) for i in 1:n)
    colour = Dict(1 => a[slot_of_leaf[1]], 2 => a[slot_of_leaf[2]])
    for i in 1:n                                  # color check, not guessed
        a[slot_of_leaf[i]] == colour[isodd(i) ? 1 : 2] || return nothing
    end

    # Nodes: all except `v`, then the four cluster nodes (colors translated).
    alt   = [u for u in 1:length(g.nodes) if u != v]
    idx   = Dict(u => t for (t, u) in enumerate(alt))
    basis = length(alt)
    nodes = vcat(CircularNode[g.nodes[u] for u in alt],
                 CircularNode[circular_node([colour[c] for c in arm]) for arm in nds])

    # outer port per leaf -> (new node, slot)
    port_of_leaf = Dict{Int,NodePort}()
    for (lf, n0, s0) in outer
        port_of_leaf[lf] = NodePort(basis + n0, s0)
    end
    port_of_slot = Dict(slot_of_leaf[lf] => p for (lf, p) in port_of_leaf)

    shift(p) = p isa Leaf ? p :
               (p.node == v ? port_of_slot[p.slot] : NodePort(idx[p.node], p.slot))

    edges = Edge[Edge(e.colour, shift(e.a), shift(e.b)) for e in g.edges]
    for (c, na, sa, nb, sb) in inner
        push!(edges, Edge(colour[c], NodePort(basis + na, sa), NodePort(basis + nb, sb)))
    end

    return CircularGraph(g.word, nodes, edges)
end

"""
    expand_gbraid(g::CircularGraph, variant::Int; rot::Int = 0)

Like above, but at the **first** gbraid node with eight arms. `nothing` if
there is none.
"""
function expand_gbraid(g::CircularGraph, variant::Int; rot::Int = 0)
    for (v, nd) in enumerate(g.nodes)
        (nd.kind === :braid && arm_count(nd) == 8) || continue
        return expand_gbraid(g, v, variant; rot = rot)
    end
    return nothing
end

# ---- THE TOWER PREIMAGE for arbitrary arm count --------------------------
#
# The fusion cluster's placement is FORCED: three cyclically adjacent slots at
# both nodes, wired in opposite direction, the middle channel direct, both outer
# channels each through a pure trivalent node in the color of the outer slots.
# `gbraid_tower_cluster` writes that cluster down directly.
#
# THE WIRING DERIVATION, so it stays checkable (A left, B right, t₁ top, t₂
# bottom; arms at the node run cw, leaves at the boundary run ccw):
#
#   * A has `2k` arms, channels at slots 1,2,3 (cw ⇒ 1 top-right, 2 right,
#     3 bottom-right). Colors `A[i] = i odd ? s : t`.
#   * B has `2l` arms, channels at 3,2,1 — OPPOSITE direction, the house
#     convention "a bigon closes its edges in opposite directions".
#     `A[1]=s=B[3]`, `A[2]=t=B[2]`, `A[3]=s=B[1]` ✓.
#   * Channel 2 (color `t`) runs DIRECTLY. Channels 1 and 3 (color `s`) each
#     go through a pure `[s,s,s]` trivalent node; both stems point outward.
#   * Boundary ccw starting at the stem of t₁: `t₁, A[2k], A[2k-1], …, A[4],
#     t₂, B[2l], B[2l-1], …, B[4]` — length `1 + (2k-3) + 1 + (2l-3) =
#     2(k+l-2)` ✓, and the sequence alternates `s,t,s,t,…`, starting with `s`
#     as the bare node does. (That slots run DESCENDING per node is literally
#     `check_wiring`'s `:leaf_ring` condition.)
#
# Over all choices of `qB`, wiring direction, and both trivalent rotations, only
# the opposite-direction wirings survive, and they ALL give the same
# `circular_canonical_key`, so the construction below is unique. For `(k,l) ∈
# {(3,3),(4,3),(3,4),(4,4),(5,3),(3,5)}`: `euler = 2`, `check_wiring` empty,
# `contract_circular_cluster` returns exactly the bare `2(k+l-2)`-armed node, and
# the degree is preserved (`-2 ↦ -2`, `-4 ↦ -4`, `-6 ↦ -6`).

"""
    gbraid_tower_decompositions(m::Int) -> Vector{Tuple{Int,Int}}

All decompositions `m = k + l - 2` with `k, l ≥ 3` — the possible expansions
of a `2m`-armed gbraid into two braid-like nodes `2k`, `2l`.

    m = 4 ⇒ [(3,3)]          m = 5 ⇒ [(3,4),(4,3)]
    m = 6 ⇒ [(3,5),(4,4),(5,3)]
"""
gbraid_tower_decompositions(m::Int) = [(k, m + 2 - k) for k in 3:(m - 1)]

"""
    gbraid_tower_literal(k::Int, l::Int) -> (nodes, inner, outer)

The fusion cluster `(2k, 2l)` as a cluster literal in the shape of
[`GBRAID_PREIMAGES`](@ref) — colors `1` (that of the outer channels and
trivalent nodes) and `2`. Node `1 = A` (`2k` arms), `2 = B` (`2l` arms),
`3 = t₁`, `4 = t₂`. See the comment block above for the derivation.

Throws for `k < 3` or `l < 3`.
"""
function gbraid_tower_literal(k::Int, l::Int)
    (k >= 3 && l >= 3) ||
        throw(ArgumentError("gbraid_tower_literal: k, l must be ≥ 3, were $k, $l"))
    nA, nB = 2k, 2l
    A = [isodd(i) ? 1 : 2 for i in 1:nA]
    B = [isodd(j) ? 1 : 2 for j in 1:nB]
    nodes = [A, B, [1, 1, 1], [1, 1, 1]]

    # bslots[i] = B-slot of the i-th channel; opposite direction from A-slot i.
    bslots = [3, 2, 1]
    inner = [(2, 1, 2, 2, bslots[2]),          # channel 2: A-B direct
             (1, 1, 1, 3, 3), (1, 3, 2, 2, bslots[1]),   # channel 1 via t₁
             (1, 1, 3, 4, 2), (1, 4, 3, 2, bslots[3])]   # channel 3 via t₂

    outer = Tuple{Int,Int,Int}[(1, 3, 1)]                # leaf 1 = stem of t₁
    for i in nA:-1:4;  push!(outer, (length(outer) + 1, 1, i)); end
    push!(outer, (length(outer) + 1, 4, 1))              # stem of t₂
    for j in 1:(nB - 3)
        push!(outer, (length(outer) + 1, 2, mod1(bslots[3] - j, nB)))
    end
    @assert length(outer) == nA + nB - 4
    return (nodes, inner, outer)
end

"""
    gbraid_tower_cluster(k::Int, l::Int; s = 1, t = 2) -> CircularGraph

The fusion cluster `(2k, 2l)` as a STANDALONE diagram with `2(k+l-2)`
boundary leaves (color pair `(s,t)`). A check/notebook tool only:
`contract_circular_cluster(g, [1,2,3,4])` recovers the bare `2(k+l-2)`-armed
gbraid from it.
"""
function gbraid_tower_cluster(k::Int, l::Int; s::Int = 1, t::Int = 2)
    (nds, inner, outer) = gbraid_tower_literal(k, l)
    colour = Dict(1 => s, 2 => t)
    nodes = CircularNode[circular_node([colour[c] for c in arm]) for arm in nds]
    edges = Edge[Edge(colour[c], NodePort(na, sa), NodePort(nb, sb))
                 for (c, na, sa, nb, sb) in inner]
    word = Int[]
    for (lf, nn, ss) in outer
        push!(edges, Edge(arm_colour(nodes[nn], ss), Leaf(lf), NodePort(nn, ss)))
        push!(word, arm_colour(nodes[nn], ss))
    end
    return CircularGraph(CircularWord(word), nodes, edges)
end

"""
    expand_gbraid_tower(g::CircularGraph, v::Int, k::Int, l::Int; rot::Int = 0)
        -> Union{Nothing, CircularGraph}

**C13⁻¹ for arbitrary arm count.** Replaces the `2(k+l-2)`-armed gbraid
`v` with the tower cluster `(2k, 2l)` (see `gbraid_tower_literal`). `rot`
rotates which arm of the node corresponds to the cluster's first leaf, same
as in [`expand_gbraid`](@ref).

`nothing` if `v` is not a gbraid with `2(k+l-2)` arms, or the color check
fails. **Not in `CIRCULAR_RULES`** — this rule expands and increases
`circular_weight`.
"""
function expand_gbraid_tower(g::CircularGraph, v::Int, k::Int, l::Int; rot::Int = 0)
    nd = g.nodes[v]
    (nd.kind === :braid && arm_count(nd) == 2 * (k + l - 2)) || return nothing
    return _circular_splice_cluster(g, v, gbraid_tower_literal(k, l); rot = rot)
end

"""
    expand_gbraid_variants(g::CircularGraph, v::Int) -> Vector{CircularGraph}

**All** planar expansions of the gbraid `v`, deduplicated by
`circular_canonical_key`: for eight arms, the four `GBRAID_PREIMAGES`
(eight rotations each); for `2m` arms, the tower clusters over all
decompositions `m = k + l - 2` (`2m` rotations each). Analogous to
`expand_circular_braid_relation_variants`.

Empty vector if `v` is not a gbraid.
"""
function expand_gbraid_variants(g::CircularGraph, v::Int)
    nd = g.nodes[v]
    # The `2k`-nodes with `2k ≥ 8` can
    # be expanded — the 6-armed braid is itself a generator.
    (nd.kind === :braid && arm_count(nd) >= 8) || return CircularGraph[]
    n = arm_count(nd); m = n ÷ 2
    out = CircularGraph[]; seen = Set{Any}()
    kand = CircularGraph[]
    if n == 8
        for variant in 1:length(GBRAID_PREIMAGES), r in 0:(n - 1)
            h = expand_gbraid(g, v, variant; rot = r)
            h === nothing || push!(kand, h)
        end
    end
    for (k, l) in gbraid_tower_decompositions(m), r in 0:(n - 1)
        h = expand_gbraid_tower(g, v, k, l; rot = r)
        h === nothing || push!(kand, h)
    end
    for h in kand
        kk = circular_canonical_key(h)
        kk in seen && continue
        push!(seen, kk); push!(out, h)
    end
    return out
end
