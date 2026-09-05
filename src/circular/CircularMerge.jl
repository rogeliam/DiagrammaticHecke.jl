# circular/CircularMerge.jl — the merge machinery: merge_at_edge, merge_nodes,
# merge_at_edges, _merge_at, _merge_remap_port, merge_all.
#
# ---- merging --------------------------------------------------------------------
#
# This implements the OPERATION; WHEN to merge is decided by the reduction
# rules elsewhere. `merge_nodes` is exactly the step that turns "two trivalents
# + one connecting edge" into ONE node; because `CircularNode` equality only sees
# the cyclic arm sequence, the two bracketings of `111 → 1` become the same
# node in the process.
#
# "gbraid" (general braid) below means a `:braid` node with `2k >= 8` arms;
# the matchers check `kind === :braid && arm_count == 8`.

"""
    merge_at_edge(g::CircularGraph, ei::Int) -> CircularGraph

Contracts the CONCRETE edge `g.edges[ei]` — both ends must be `NodePort`s at
DIFFERENT, non-`:braid` nodes `i` (slot `si`) and `j` (slot `sj`). The new arm
sequence is `i`'s arms starting at `si+1`, clockwise (length `deg_i − 1`),
followed by `j`'s arms starting at `sj+1` (length `deg_j − 1`) — the unique
planar splice: at the contracted edge the two node boundaries meet, and the
clockwise order around the new node runs through the rest of `i` first, then
the rest of `j`.

Unlike [`merge_nodes`](@ref), this function does NOT require `i` and `j` to be
connected by only ONE edge — multiply connected nodes (e.g. a bigon edge
between two `[1,1,1]`) are explicitly allowed: a SURVIVING edge that ends up
with BOTH ends at the new node after contraction automatically becomes a
SELF-CONNECTION (`_merge_remap_port` already handles this case unchanged) —
exactly the mechanism that subsumes R11 (bigon → needle) without a rule branch
of its own.

Throws if `ei` is not a valid edge between two different non-`:braid` nodes.
"""
function merge_at_edge(g::CircularGraph, ei::Int)
    1 <= ei <= length(g.edges) || throw(ArgumentError("merge_at_edge: edge index $ei outside 1:$(length(g.edges))"))
    e = g.edges[ei]
    (e.a isa NodePort && e.b isa NodePort) || throw(ArgumentError(
        "merge_at_edge: edge $ei does not have a NodePort end on both sides"))
    i, si = e.a.node, e.a.slot
    j, sj = e.b.node, e.b.slot
    i == j && throw(ArgumentError("merge_at_edge: edge $ei is already a self-connection (i == j)"))
    ndi, ndj = g.nodes[i], g.nodes[j]
    (_circular_braidlike(ndi) || _circular_braidlike(ndj)) &&
        throw(ArgumentError("merge_at_edge: braid-like nodes cannot be merged"))
    return _merge_at(g, i, si, j, sj, ei)
end

"""
    merge_nodes(g::CircularGraph, i::Int, j::Int) -> CircularGraph

Contracts the ONE edge between nodes `i` and `j` into a single node — a thin
wrapper around [`merge_at_edge`](@ref) that first looks up the (unique)
connecting edge. Throws if `i` and `j` are not connected by EXACTLY ONE edge:
with multiple edges the ambiguity is REAL (contracting edge A vs. edge B gives
different wirings), and a function taking only `(i,j)` cannot express that —
use `merge_at_edge` with the concrete edge index in that case.
"""
function merge_nodes(g::CircularGraph, i::Int, j::Int)
    i == j && throw(ArgumentError("merge_nodes: i and j must be different"))
    ndi, ndj = g.nodes[i], g.nodes[j]
    (_circular_braidlike(ndi) || _circular_braidlike(ndj)) &&
        throw(ArgumentError("merge_nodes: braid-like nodes cannot be merged"))

    conn = Int[]                          # edge indices between i and j
    for (ei, e) in enumerate(g.edges)
        a, b = e.a, e.b
        if a isa NodePort && b isa NodePort
            if (a.node == i && b.node == j) || (a.node == j && b.node == i)
                push!(conn, ei)
            end
        end
    end
    length(conn) == 1 || throw(ArgumentError(
        "merge_nodes: $i and $j are not connected by exactly one edge " *
        "(have $(length(conn))) — use merge_at_edge for multiple edges"))
    return merge_at_edge(g, conn[1])
end

# ---- merging with EXACTLY TWO adjacent connections ----------------------------
#
# THE RULE: "For exactly 2 connections we proceed
# similarly. Both get deleted. So `1....(k,k+1)....n` connected to
# `1.....(l,l+1)....m` becomes
# `1....k-1, l+2.....n.....l-1, k+2....n`."
#
# READING (machine-verified): read cyclically, this is the exact analogue of
# the one-connection case
# (`_merge_at` above), except TWO adjacent slots are consumed on EACH side
# instead of one: **`A`'s arms from `k+2` forward (`n−2` of them), then `B`'s
# arms from `l+2` forward (`m−2` of them)**. This reading matches exactly what
# TWO successive `_merge_at` contractions would give, PROVIDED `_fr_needle` is
# disabled in between — the self-connection arising after the first
# contraction always sits at two cyclically adjacent slots (`gap == d−1`),
# which is exactly the condition the adjacency requirement below enforces.
#
# COUNTER-ROTATION is mandatory, not optional (machine-checked): the pair on
# `A` sits at `(k, k+1)`, the CONNECTED pair on `B` must sit at `(l+1, l)`
# (i.e. `k↔l+1` AND `k+1↔l`) — not `k↔l, k+1↔l+1` (same-rotation). With
# DIFFERENT colours on the two connections, colour fidelity often forces this
# already; with the SAME colours both assignments are colour-valid, but only
# the counter-rotating one yields an ADJACENT self-connection after the first
# contraction (the other stays topologically NON-adjacent, and the boundary
# would fall apart into two separate arcs — the same trap as the
# `double_leaf(12321,…)` case, but here produced by wrong orientation rather
# than genuine non-adjacency). This function therefore
# checks BOTH assignments for colour fidelity and requires the
# counter-rotating one to actually hold.

"""
    merge_at_edges(g::CircularGraph, ei::Int, ej::Int) -> CircularGraph

Contracts TWO CONCRETE edges `g.edges[ei]` and `g.edges[ej]` at once — the
two-connection case, analogue of [`merge_at_edge`](@ref) for the
one-connection case. Preconditions (violation throws `ArgumentError`, NOT
a silently different result):

- both edges have `NodePort` ends at the **same** two nodes `i ≠ j`
  (different, non-`:braid`);
- the two slots of `i` (`si`, `sj_i`, one from each edge) are cyclically
  ADJACENT (`sj_i == si+1` or vice versa), likewise on `j`;
- the assignment is COUNTER-ROTATING: if `A`'s first pair slot `si` connects
  to `B`'s slot `t`, `A`'s second pair slot `si+1` must connect to `B`'s OTHER
  pair slot (`k↔l+1`, `k+1↔l`) — the same-rotation assignment is REJECTED (it
  produces a topologically inconsistent boundary, see the file header
  comment).

Result: ONE node with arm sequence `A from k+2 forward (n−2 of them), then B
from l+2 forward (m−2 of them)` — both connecting edges disappear, no
self-connection remains. All remaining arms keep their previous wiring
(slot remap analogous to `_merge_at`/`_merge_remap_port`).
"""
function merge_at_edges(g::CircularGraph, ei::Int, ej::Int)
    ei != ej || throw(ArgumentError("merge_at_edges: ei and ej must be different"))
    (1 <= ei <= length(g.edges) && 1 <= ej <= length(g.edges)) ||
        throw(ArgumentError("merge_at_edges: edge index outside 1:$(length(g.edges))"))
    e1, e2 = g.edges[ei], g.edges[ej]
    (e1.a isa NodePort && e1.b isa NodePort) || throw(ArgumentError(
        "merge_at_edges: edge $ei does not have a NodePort end on both sides"))
    (e2.a isa NodePort && e2.b isa NodePort) || throw(ArgumentError(
        "merge_at_edges: edge $ej does not have a NodePort end on both sides"))

    # Both edges must connect the same two nodes.
    nodes1 = Set((e1.a.node, e1.b.node))
    nodes2 = Set((e2.a.node, e2.b.node))
    nodes1 == nodes2 || throw(ArgumentError(
        "merge_at_edges: the two edges do not connect the same two nodes"))
    length(nodes1) == 2 || throw(ArgumentError(
        "merge_at_edges: the edges are already self-connections (i == j)"))
    i, j = e1.a.node, e1.b.node
    ndi, ndj = g.nodes[i], g.nodes[j]
    (_circular_braidlike(ndi) || _circular_braidlike(ndj)) &&
        throw(ArgumentError("merge_at_edges: braid-like nodes cannot be merged"))
    di, dj = arm_count(ndi), arm_count(ndj)

    # Extract the slot of i/j for each edge.
    si1, sj1 = e1.a.node == i ? (e1.a.slot, e1.b.slot) : (e1.b.slot, e1.a.slot)
    si2, sj2 = e2.a.node == i ? (e2.a.slot, e2.b.slot) : (e2.b.slot, e2.a.slot)

    # Adjacency on i.
    si, si_next = if mod1(si1 + 1, di) == si2
        si1, si2
    elseif mod1(si2 + 1, di) == si1
        si2, si1
    else
        throw(ArgumentError("merge_at_edges: slots $si1/$si2 of node $i are not cyclically adjacent"))
    end
    # the two far slots on j correspondingly — which belongs to si, which to si_next.
    t_of_si = si == si1 ? sj1 : sj2
    t_of_si_next = si_next == si1 ? sj1 : sj2

    # Counter-rotation: the j-slot belonging to si_next must be the (cyclic)
    # PREDECESSOR of the j-slot belonging to si — i.e. sj = t_of_si,
    # sj_next = t_of_si_next with mod1(sj_next + 1, dj) == sj.
    mod1(t_of_si_next + 1, dj) == t_of_si ||
        throw(ArgumentError(
            "merge_at_edges: the connection is not counter-rotating " *
            "(expected k↔l+1, k+1↔l) — same-rotation assignment is not supported"))
    sj, sj_next = t_of_si_next, t_of_si    # sj = l (predecessor), sj_next = l+1

    # Colour fidelity (in addition to the structural check; is_wired normally
    # already guarantees this, but this function must not rely on that alone
    # — see the file header, "better to throw than to guess").
    arm_colour(ndi, si) == arm_colour(ndj, sj_next) ||
        throw(ArgumentError("merge_at_edges: colour fidelity violated at slot $i.$si / $j.$sj_next"))
    arm_colour(ndi, si_next) == arm_colour(ndj, sj) ||
        throw(ArgumentError("merge_at_edges: colour fidelity violated at slot $i.$si_next / $j.$sj"))

    # New arm sequence: A from si_next+1 (= si+2) forward (di-2 of them), then B
    # from sj_next+1 (= sj+2, since sj_next = mod1(sj+...)) — IMPORTANT: "from
    # l+2" means from the slot AFTER the pair (l, l+1) = (sj, sj_next); the
    # slot after the pair is mod1(sj_next+1, dj), which by construction equals
    # l+2 when sj=l, sj_next=l+1.
    new_arms = Int[]
    for off in 1:(di - 2)
        push!(new_arms, arm_colour(ndi, si_next + off))
    end
    for off in 1:(dj - 2)
        push!(new_arms, arm_colour(ndj, sj_next + off))
    end

    # NO back door: if the arm sequence violates the CircularNode invariant,
    # `circular_node` throws — and it should ("try/catch out"). An
    # `_unchecked_circularnode(:mixed, …)` here would turn a reported invariant
    # violation into a silent invalid graph (e.g. `[1,1,3]`, two colours with
    # only ONE 3). Anyone who lands here has found a missing rule — not an
    # overly strict constructor.
    new_node = circular_node(new_arms)

    # Slot remap: old slot (≠ the pair slots) → slot in the new node.
    remap_i = Dict{Int,Int}()
    for off in 1:(di - 2)
        remap_i[mod1(si_next + off, di)] = off
    end
    remap_j = Dict{Int,Int}()
    for off in 1:(dj - 2)
        remap_j[mod1(sj_next + off, dj)] = off + (di - 2)
    end

    new_index = idx -> idx > j ? idx - 1 : idx
    merged_index = i > j ? i - 1 : i

    new_nodes = CircularNode[]
    for (idx, nd) in enumerate(g.nodes)
        idx == j && continue
        push!(new_nodes, idx == i ? new_node : nd)
    end

    new_edges = Edge[]
    for (idx, e) in enumerate(g.edges)
        (idx == ei || idx == ej) && continue      # both consumed edges disappear
        a = _merge_remap_port(e.a, i, j, remap_i, remap_j, merged_index, new_index)
        b = _merge_remap_port(e.b, i, j, remap_i, remap_j, merged_index, new_index)
        push!(new_edges, Edge(e.colour, a, b))
    end

    return CircularGraph(g.word, new_nodes, new_edges)
end

# The shared body of merge_at_edge/merge_nodes: merges `i` (slot `si`) and `j`
# (slot `sj`) via the edge `conn_edge`.
function _merge_at(g::CircularGraph, i::Int, si::Int, j::Int, sj::Int, conn_edge::Int)
    ndi, ndj = g.nodes[i], g.nodes[j]
    di, dj = arm_count(ndi), arm_count(ndj)

    new_arms = Int[]
    for off in 1:(di - 1)
        push!(new_arms, arm_colour(ndi, si + off))
    end
    for off in 1:(dj - 1)
        push!(new_arms, arm_colour(ndj, sj + off))
    end

    # Try to construct a new CircularNode; if invalid (e.g. mixed colours but too few
    # occurrences), fall back to splicing logic to remove the created small node.
    new_node = nothing
    try
        new_node = circular_node(new_arms)
    catch err
        # Splice fallback: gather the far ports (other ends of the incident slots)
        far_ports = Port[]
        far_colours = Int[]
        # helper to find edge incident to a given node slot (excluding conn_edge)
        function find_far(ei_skip, node_idx, slot_target)
            for (ei, e) in enumerate(g.edges)
                ei == ei_skip && continue
                if e.a isa NodePort && e.a.node == node_idx && e.a.slot == slot_target
                    return (e.b, e.colour)
                elseif e.b isa NodePort && e.b.node == node_idx && e.b.slot == slot_target
                    return (e.a, e.colour)
                end
            end
            return (nothing, nothing)
        end
        for off in 1:(di - 1)
            slot = mod1(si + off, di)
            p, col = find_far(conn_edge, i, slot)
            p === nothing && continue
            push!(far_ports, p)
            push!(far_colours, col)
        end
        for off in 1:(dj - 1)
            slot = mod1(sj + off, dj)
            p, col = find_far(conn_edge, j, slot)
            p === nothing && continue
            push!(far_ports, p)
            push!(far_colours, col)
        end

        nfar = length(far_ports)
        # Build new node list with BOTH i and j removed
        function new_index_remove_two(idx)
            if idx == i || idx == j
                error("invalid mapping for removed node")
            end
            # assume i < j for simplicity of indexing adjustments
            if i < j
                return idx > j ? idx - 2 : (idx > i ? idx - 1 : idx)
            else
                return idx > i ? idx - 2 : (idx > j ? idx - 1 : idx)
            end
        end

        new_nodes = CircularNode[]
        for (idx, nd) in enumerate(g.nodes)
            if idx == i || idx == j
                continue
            end
            push!(new_nodes, nd)
        end

        # Rebuild edges: keep edges that are not incident to i or j and not the conn_edge
        new_edges = Edge[]
        for (ei, e) in enumerate(g.edges)
            ei == conn_edge && continue
            # skip edges incident to i or j
            a_i = (e.a isa NodePort && (e.a.node == i || e.a.node == j))
            b_i = (e.b isa NodePort && (e.b.node == i || e.b.node == j))
            if a_i || b_i
                continue
            end
            # remap node indices
            a = e.a
            b = e.b
            if a isa NodePort
                a = NodePort(new_index_remove_two(a.node), a.slot)
            end
            if b isa NodePort
                b = NodePort(new_index_remove_two(b.node), b.slot)
            end
            push!(new_edges, Edge(e.colour, a, b))
        end

        # If exactly two far ports remain, splice them into a single edge.
        if nfar == 2
            p1, p2 = far_ports[1], far_ports[2]
            col = far_colours[1]
            # remap NodePorts in p1/p2
            if p1 isa NodePort
                p1 = NodePort(new_index_remove_two(p1.node), p1.slot)
            end
            if p2 isa NodePort
                p2 = NodePort(new_index_remove_two(p2.node), p2.slot)
            end
            push!(new_edges, Edge(col, p1, p2))
            return CircularGraph(g.word, new_nodes, new_edges)
        end

        # If no far ports remain, just return the graph with i and j removed.
        if nfar == 0
            return CircularGraph(g.word, new_nodes, new_edges)
        end

        # nfar >= 3: there is nothing to splice here, and there is NO unchecked
        # back door: `_unchecked_circularnode(:mixed, new_arms)` would silently
        # build nodes `circular_node` rightly rejects. `circular_node`'s error
        # propagates instead; anyone who lands here needs a rule BEFORE `merge`
        # (e.g. C2 `_fr_dot_walks_on`), not a softer constructor.
        rethrow(err)
    end

    # old slot (≠ the contracted one) → slot in the new node
    remap_i = Dict{Int,Int}()
    for off in 1:(di - 1)
        remap_i[mod1(si + off, di)] = off
    end
    remap_j = Dict{Int,Int}()
    for off in 1:(dj - 1)
        remap_j[mod1(sj + off, dj)] = off + (di - 1)
    end

    # j disappears; the new node sits at i's position, indices > j shift down.
    new_index = idx -> idx > j ? idx - 1 : idx
    merged_index = i > j ? i - 1 : i

    new_nodes = CircularNode[]
    for (idx, nd) in enumerate(g.nodes)
        idx == j && continue
        push!(new_nodes, idx == i ? new_node : nd)
    end

    new_edges = Edge[]
    for (ei, e) in enumerate(g.edges)
        ei == conn_edge && continue      # the contracted edge disappears
        a = _merge_remap_port(e.a, i, j, remap_i, remap_j, merged_index, new_index)
        b = _merge_remap_port(e.b, i, j, remap_i, remap_j, merged_index, new_index)
        push!(new_edges, Edge(e.colour, a, b))
    end

    return CircularGraph(g.word, new_nodes, new_edges)
end

# Remaps a single port onto the new node indices/slots.
function _merge_remap_port(p::Port, i::Int, j::Int,
                           remap_i::Dict{Int,Int}, remap_j::Dict{Int,Int},
                           merged_index::Int, new_index)
    p isa NodePort || return p
    if p.node == i
        return NodePort(merged_index, remap_i[p.slot])
    elseif p.node == j
        return NodePort(merged_index, remap_j[p.slot])
    else
        return NodePort(new_index(p.node), p.slot)
    end
end

"""
    merge_all(g::CircularGraph) -> CircularGraph

Fixed point of [`merge_at_edge`](@ref) over all edges between two DIFFERENT
non-`:braid` nodes. EDGE-DRIVEN, not via `nshared == 1`: iterating by shared
edge count instead would let a bigon edge (two trivalents, doubly connected)
survive, because the first contraction turns the second edge into a
self-connection, which no `(i,j)` pair would find afterwards. Terminates
because every merge lowers the node count by one.
"""
function merge_all(g::CircularGraph)
    cur = g
    while true
        merged = false
        for (ei, e) in enumerate(cur.edges)
            e.a isa NodePort && e.b isa NodePort || continue
            e.a.node == e.b.node && continue        # already a self-connection
            _circular_braidlike(cur.nodes[e.a.node]) && continue
            _circular_braidlike(cur.nodes[e.b.node]) && continue
            cur = merge_at_edge(cur, ei)
            merged = true
            break
        end
        merged || return cur
    end
end
