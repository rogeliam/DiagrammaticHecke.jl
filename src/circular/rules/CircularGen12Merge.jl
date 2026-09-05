# circular/rules/CircularGen12Merge.jl — C13 `gen12_merge`: two braid-like nodes
# (`2k` + `2l` arms) + 2 trivalent nodes ↦ ONE gbraid with `2(k+l-2)` arms.
#
# ARM-COUNT GENERIC: the contraction itself
# (`contract_circular_cluster` doesn't know about `kind` at all) works for any
# `2k`/`2l`; only the coarse filter (`_circular_braidlike`) and the final
# check (`arm_count == expected`) are specific to the braid-like shape. The
# case `:braid` + `:braid` ↦ gbraid₈ (k = l = 3) is one instance of the
# general family.
#
# THE MATH behind it: the `2k`-armed node is `π_k = id_{B_{w₀}}` between the
# extreme shifts, degree `6 - 2k`, and that degree space is ONE-DIMENSIONAL.
# Arm and degree balance:
#
#     arms:   2k + 2l + 3 + 3 - 2·5 = 2(k+l-2)                     ✓
#     degree: (6-2k) + (6-2l) + (-1) + (-1) = 6 - 2(k+l-2)         ✓
#
# The degree balance requires a gbraid's `circular_degree = 6 - #arms`; with
# `circular_degree = 2 - #arms/2` the degree would not balance and the map
# would not be a relation. The composite is `≠ 0`
# (`reduce_to_circular_leave`: one term each, coefficient `+1`), and `(4,3)`
# is the double tower `π₅`.
#
# WHEN IT FIRES. Two braid-like nodes (per
# `_circular_braidlike`: `:braid` with 6 arms, or the general `2k >= 8`-armed
# one) with the same color pair are connected by three channels;
# if BOTH OUTER channels each carry a trivalent role (`t = 1 + 1`), the
# cluster `{b₁, b₂, t₁, t₂}` contracts to a single gbraid. The arm balance
# is forced by the match condition (5 inner edges): `2k + 2l + 3 + 3 - 2·5 =
# 2(k+l-2)`.
#
# THE PLACEMENT IS FORCED: three cyclically ADJACENT slots at both nodes, wired in opposite
# directions; the MIDDLE channel runs directly, both OUTER channels each go
# through a pure trivalent node in the color of the outer slots (so both
# trivalent nodes share a color). Any other placement (trivalent at the
# middle channel, gap in the bundle) has ZERO planar embeddings — the free
# stem would land in an enclosed region. The coarse filter below therefore
# does not need to check for this placement at all: anything not placed this
# way simply fails to contract to a valid gbraid.
#
# WHEN IT DOES NOT FIRE. If both trivalent roles sit
# in ONE channel (a channel node with `nc ≥ 4`), the contraction does run and
# produces eight arms — but NOT alternating (`[3,3,2,3,2,2,3,2]`), so not a
# valid `CircularNode`. Deleting the
# foreign-colored letters from the input diagram's boundary word leaves
# exactly this sequence — the BOUNDARY itself doesn't alternate, so no node
# with alternating arms can carry it. This case stays with
# C8 (`_fr_braid_relation_gen`), where `CircularBraidChannels.jl` already handles
# it. Here it is not matched.
#
# CHANNEL NORMALIZATION, the same pattern as `_fr_braid_relation_gen`: try
# the bare rule first; otherwise pull a pure trivalent node out of the
# general-13 node with `_circular_pull_trivalent`, at BOTH channels in turn — that
# is the only difference from C8. The rule preserves the boundary word and
# `circular_degree`, keeps `euler = 2`, passes `check_wiring`, and gives a
# result unique up to `circular_canonical_key`. Without the
# channel normalization, such diagrams would be
# FIXED POINTS of `CIRCULAR_RULES`.
#
# TERMINATION. Component 1 of `circular_weight` is the ARM SUM over braid-like
# nodes (`_circular_braidlike`), and the rule strictly decreases it:
# `2k + 2l ↦ 2(k+l-2)`, i.e. `-4` in every case — even when both partners are
# gbraids. (For `k = l = 3` specifically, that is `12 ↦ 8`; the weight goes
# `(2,4,22) ↦ (0,3,16)` with a gen-13 node, `(2,4,18) ↦ (0,1,8)` in the
# classical case — arm-sum component `12 ↦ 8` in both.) So the assert in
# `CircularDriver.jl` is satisfied without touching `circular_weight` itself (unlike
# C14).
#
# PLACE IN `CIRCULAR_RULES`: AFTER C7/C8. Placed before them, it would swallow
# sites where F2-2/F2-3 is actually the right rule.
#
# ⚠ The gbraid is locked out of the generic merge machinery
# (`_circular_braidlike`), and `merge_at_edge` throws at `:braid` nodes. So this
# file builds the contraction itself — that's exactly what
# `_circular_raw_merge_at` is for.
#
# FORWARD REFERENCE to `expand_circular_merge` (via `_circular_pull_trivalent`,
# CircularBraidChannels.jl): same situation as there, see that file's header.

"""
    _circular_raw_merge_at(g::CircularGraph, ei::Int) -> (CircularGraph, Int)

Like `merge_at_edge`, but WITHOUT the `:braid` lock and without the
`circular_node` check. Meant only as an intermediate step of the contraction; the
final result is validated with `circular_node` in `contract_circular_cluster`. The
splice formula is literally `_merge_at` (circular/CircularMerge.jl): the merged node's
arm sequence is "`i`'s arms from `sᵢ+1` on, then `j`'s arms from `sⱼ+1` on" —
purely planar, doesn't know about `kind` at all.
"""
function _circular_raw_merge_at(g::CircularGraph, ei::Int)
    e = g.edges[ei]
    (e.a isa NodePort && e.b isa NodePort) || error("_circular_raw_merge_at: not an inner edge")
    i, si = e.a.node, e.a.slot
    j, sj = e.b.node, e.b.slot
    i == j && error("_circular_raw_merge_at: self-connection")
    ndi, ndj = g.nodes[i], g.nodes[j]
    di, dj = arm_count(ndi), arm_count(ndj)

    new_arms = Int[]
    for off in 1:(di - 1); push!(new_arms, arm_colour(ndi, si + off)); end
    for off in 1:(dj - 1); push!(new_arms, arm_colour(ndj, sj + off)); end

    remap_i = Dict{Int,Int}(); remap_j = Dict{Int,Int}()
    for off in 1:(di - 1); remap_i[mod1(si + off, di)] = off; end
    for off in 1:(dj - 1); remap_j[mod1(sj + off, dj)] = off + (di - 1); end

    new_index = idx -> idx > j ? idx - 1 : idx
    merged = i > j ? i - 1 : i
    remap(p) = !(p isa NodePort) ? p :
        p.node == i ? NodePort(merged, remap_i[p.slot]) :
        p.node == j ? NodePort(merged, remap_j[p.slot]) :
        NodePort(new_index(p.node), p.slot)

    new_nodes = CircularNode[]
    for (idx, nd) in enumerate(g.nodes)
        idx == j && continue
        push!(new_nodes, idx == i ? _unchecked_circularnode(:mixed, new_arms) : nd)
    end
    new_edges = Edge[Edge(ee.colour, remap(ee.a), remap(ee.b))
                     for (k, ee) in enumerate(g.edges) if k != ei]
    return CircularGraph(g.word, new_nodes, new_edges), merged
end

"""
    _circular_drop_adjacent_selfloop(g::CircularGraph, v::Int) -> (CircularGraph, Bool)

Deletes a self-connection at `v` at cyclically ADJACENT slots, together with
those two slots. This is the two-connection convention of `merge_at_edges`
(documented in the header of circular/CircularMerge.jl): after the first contraction,
every further connection within the cluster remains as such a self-loop.

The tension with C1 (`_fr_needle` says: self-connection ⇒ 0) is not a new
convention here — it already stands in the header of `CircularMerge.jl`. A
NON-adjacent self-connection, by contrast, aborts the contraction (the
cluster was not a disk after all).
"""
function _circular_drop_adjacent_selfloop(g::CircularGraph, v::Int)
    nd = g.nodes[v]; d = arm_count(nd)
    for (ei, e) in enumerate(g.edges)
        (e.a isa NodePort && e.b isa NodePort && e.a.node == v && e.b.node == v) || continue
        s1, s2 = e.a.slot, e.b.slot
        s, snext = mod1(s1 + 1, d) == s2 ? (s1, s2) :
                   mod1(s2 + 1, d) == s1 ? (s2, s1) : (0, 0)
        s == 0 && continue
        keep = [mod1(snext + off, d) for off in 1:(d - 2)]
        remap = Dict(keep[k] => k for k in eachindex(keep))
        new_nodes = copy(g.nodes)
        new_nodes[v] = _unchecked_circularnode(:mixed, [arm_colour(nd, snext + off) for off in 1:(d - 2)])
        rm(p) = p isa NodePort && p.node == v ? NodePort(v, remap[p.slot]) : p
        new_edges = Edge[Edge(ee.colour, rm(ee.a), rm(ee.b))
                         for (k, ee) in enumerate(g.edges) if k != ei]
        return CircularGraph(g.word, new_nodes, new_edges), true
    end
    return g, false
end

"""
    contract_circular_cluster(g::CircularGraph, cluster::Vector{Int}) -> nothing | (CircularGraph, Int)

Contracts the node set `cluster` to ONE node and returns the new graph
together with that node's index. `nothing` if the result is not a valid
`CircularNode` arm sequence (the cluster was then not a disk), or if a
non-adjacent self-connection remains.
"""
function contract_circular_cluster(g::CircularGraph, cluster::Vector{Int})
    cur = g; live = collect(cluster)
    while true
        hit = nothing
        for (ei, e) in enumerate(cur.edges)
            (e.a isa NodePort && e.b isa NodePort) || continue
            (e.a.node in live && e.b.node in live && e.a.node != e.b.node) || continue
            hit = ei; break
        end
        hit === nothing && break
        e = cur.edges[hit]; i, j = e.a.node, e.b.node
        cur, merged = _circular_raw_merge_at(cur, hit)
        live = sort(unique([(k == i || k == j) ? merged : (k > j ? k - 1 : k) for k in live]))
    end
    length(live) == 1 || return nothing
    v = live[1]
    while true
        cur, ok = _circular_drop_adjacent_selfloop(cur, v)
        ok || break
    end
    for e in cur.edges
        (e.a isa NodePort && e.b isa NodePort && e.a.node == v && e.b.node == v) && return nothing
    end
    newnd = try circular_node(collect(arms(cur.nodes[v]))) catch; return nothing end
    nodes2 = copy(cur.nodes); nodes2[v] = newnd
    return CircularGraph(cur.word, nodes2, copy(cur.edges)), v
end

"Pure trivalent node: three arms, one color."
_circular_is_pure_trivalent_node(nd::CircularNode) =
    arm_count(nd) == 3 && length(unique(arms(nd))) == 1

"""
    _circular_gen12_clusters(g::CircularGraph) -> Vector{Vector{Int}}

The COARSE FILTER for C13: two **braid-like** nodes (`_circular_braidlike`:
`:braid` with 6 arms, or the general `2k ≥ 8`-armed one) sharing a color pair,
two pure trivalent nodes with colors from that pair, **exactly 5** edges
between nodes of the cluster, and **none** of them a self-connection.

The return value is `[b₁, b₂, t₁, t₂]` with `b₁ < b₂` in the first two
slots — `_circular_gen12_contract` reads the target arm count off of those.

The coarse filter alone is NOT enough — it produces false positives. The
actual decision is made in `_fr_gen12_merge` by computing the contraction.
"""
function _circular_gen12_clusters(g::CircularGraph)
    braids = [i for (i, nd) in enumerate(g.nodes) if _circular_braidlike(nd)]
    trivs = [i for (i, nd) in enumerate(g.nodes) if _circular_is_pure_trivalent_node(nd)]
    out = Vector{Int}[]
    for b1 in braids, b2 in braids
        b1 < b2 || continue
        paar = Set(arms(g.nodes[b1])[1:2])
        paar == Set(arms(g.nodes[b2])[1:2]) || continue
        for t1 in trivs, t2 in trivs
            t1 < t2 || continue
            (arms(g.nodes[t1])[1] in paar && arms(g.nodes[t2])[1] in paar) || continue
            cl = [b1, b2, t1, t2]
            nint = 0; selfloop = false
            for e in g.edges
                (e.a isa NodePort && e.b isa NodePort) || continue
                if e.a.node in cl && e.b.node in cl
                    e.a.node == e.b.node && (selfloop = true)
                    nint += 1
                end
            end
            (nint == 5 && !selfloop) && push!(out, cl)
        end
    end
    return out
end

# The bare rule: run the coarse filter + contraction, require a gbraid
# with the TARGET arm count `2k + 2l - 4`. One term, coefficient 1.
function _circular_gen12_contract(g::CircularGraph)
    for cl in _circular_gen12_clusters(g)
        expected = arm_count(g.nodes[cl[1]]) + arm_count(g.nodes[cl[2]]) - 4
        r = contract_circular_cluster(g, cl)
        r === nothing && continue
        (g2, v) = r
        (g2.nodes[v].kind === :braid && arm_count(g2.nodes[v]) == expected) || continue
        return CircularComboR(g2)
    end
    return nothing
end

"""
    _fr_gen12_merge(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C13 `gen12_merge`.** Two braid-like nodes with `2k` resp. `2l` arms, whose
two outer channels each carry a trivalent role, merge into ONE gbraid
with `2(k+l-2)` arms. For `k = l = 3` this is the base case (2 braids ↦
gbraid₈).

Tries the bare rule first (both roles are already pure trivalent nodes),
otherwise channel normalization: pull a trivalent node out of the
general-13 node at BOTH channels (`_circular_pull_trivalent` = `C5⁻¹` here) and
try again. See the file header for the required placement, and for the
`nc ≥ 4` case, which does not match here.

⚠ The **normalization stage** only runs over pairs of plain `:braid`
nodes (`_circular_braid_pairs` in `CircularBraidChannels.jl` filters on
`kind === :braid`) — only the BARE contraction handles general `2k`/`2l`
arm counts. A `2k`/`2l` pair whose channels carry
general-13 nodes therefore stays a fixed point.
"""
function _fr_gen12_merge(g::CircularGraph)
    r = _circular_gen12_contract(g)
    r === nothing || return r
    for (a, b) in _circular_braid_pairs(g)
        g1 = _circular_pull_trivalent(g, a, b)
        g1 === nothing && continue
        r1 = _circular_gen12_contract(g1)
        r1 === nothing || return r1
        # second channel: node indices have shifted, so search again.
        for (a2, b2) in _circular_braid_pairs(g1)
            g2 = _circular_pull_trivalent(g1, a2, b2)
            g2 === nothing && continue
            r2 = _circular_gen12_contract(g2)
            r2 === nothing || return r2
        end
    end
    return nothing
end
