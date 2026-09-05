# circular/rules/CircularBraidRules.jl — C6–C8: the two-colour rules.
#
# _fr_dot_into_braid (C6), _fr_braid_back (C7), _fr_braid_relation (C8) —
# the two-colour relations R4/R5/R9 on circular nodes.
#
# ---- C6–C8: the two-colour rules ---------------------------------------
#
# The rules are stated on circular nodes. Substitutions used throughout:
#   kind === :trivalent    → arms == [c,c,c]
#   :braid && m==3         → :braid
#   :braid && m==2         → :mixed with 4 alternating 1/3 arms
#   _slot_colour           → arm_colour
#   opposite_slot          → circular_opposite_slot
#   Node(...)              → circular_node([...])
#
# ⚠ In the circular world a crossing `[1,3,1,3]` is NOT distinguishable from a
# merged 4-armed 1/3 cluster, so R12 matches more shapes than the plain
# original does.

"""
    _circular_is_pure_trivalent(nd::CircularNode) -> Bool

Pure trivalent: **three arms, one colour** — nothing more.

The check is on `arm_count` and colour uniformity ALONE, not on `nd.kind ===
:mixed`: `kind` is derived in `circular_node` from the COLOURS (all arms `2`
⇒ `:mono`, all arms in `{1,3}` ⇒ `:mixed`), so a `kind` check would really
mean "colour 1 or 3" and silently exclude a **colour-2 trivalent**
`circular_node([2,2,2])` — exactly the common case in A₃ with 12-braids,
where C8 `braid_relation` must still recognize diagrams with two
doubly-connected 12-braids plus a colour-2 trivalent as matching, not as
normal form.

This is also redundant with `arm_count` alone: `:braid` has 6 arms, gbraid
has ≥ 8, so neither can be 3-armed. `expand_circular_braid_relation_variants`
(CircularInverseRules.jl) checks the arm condition directly, and
`CircularGen12Merge.jl` carries its own copy `_circular_is_pure_trivalent_node`.
"""
_circular_is_pure_trivalent(nd::CircularNode) = arm_count(nd) == 3 && length(unique(arms(nd))) == 1

# ---- C6: _fr_dot_into_braid (R2-2/R4) -----------------------------------------
"""
    _fr_dot_into_braid(g::CircularGraph) -> Union{Nothing, CircularComboR}

R4: a dot caps
one slot of a `:braid` node (m=3, 6 arms, colour pair `{s,t}`). Rewrite to a
SUM of two diagrams: term A (two dots of the minority colour + a pure
trivalent of the majority colour) + term B (one through-strand per colour +
a dot of the majority colour on the remaining leg). `nothing` if no match or
`m != 3`.
"""
function _fr_dot_into_braid(g::CircularGraph)
    for (di, dnd) in enumerate(g.nodes)
        arm_count(dnd) == 1 || continue
        at = _circular_edges_at_node(g, di)
        length(at) == 1 || continue
        (ei, _, other) = at[1]
        other isa NodePort || continue
        other.node == di && continue
        bn = g.nodes[other.node]
        bn.kind === :braid || continue
        arm_count(bn) == 6 || continue     # six arms explicitly (kind alone allows 2k arms too)
        arm_colour(bn, other.slot) == arms(dnd)[1] || continue
        cap = other.slot
        M = 6
        minc = arm_colour(bn, cap)
        majc = minc == arms(bn)[1] ? arms(bn)[2] : arms(bn)[1]
        # Count SLOTS, not edges: a
        # SELF-LOOP at the braid occupies two slots but `_circular_edges_at_node`
        # returns only ONE entry for it — a raw edge count would come out as 5
        # instead of 6, aborting the rule with the WRONG reason ("braid not
        # fully occupied").
        slotport = _circular_slot_ports(g, other.node)
        length(slotport) == M || continue

        # A SELF-LOOP AT THE BRAID is left UNMATCHED, guarded explicitly here.
        # The rule deletes the braid node; under a self-loop the "opposite port"
        # of a slot points back at the node itself, so that arm pair has no
        # right-hand side (terms A/B would wire onto a dead node).
        any(p -> p isa NodePort && p.node == other.node, values(slotport)) && continue
        ring = [mod1(cap + k, M) for k in 1:(M - 1)]
        rleaves = [(arm_colour(bn, k), slotport[k]) for k in ring]
        maj_leaves = [p for (col, p) in rleaves if col == majc]
        min_leaves = [p for (col, p) in rleaves if col == minc]
        (length(maj_leaves) == 3 && length(min_leaves) == 2) || continue
        adjpair = nothing
        for i in 1:length(rleaves)
            col_i, p_i = rleaves[i]
            col_j, p_j = rleaves[mod1(i + 1, length(rleaves))]
            if col_i == majc && col_j == majc
                adjpair = (p_i, p_j); break
            end
        end
        adjpair === nothing && continue

        dropbase = Set{Int}([e for (e, _, _) in _circular_edges_at_node(g, other.node)])
        push!(dropbase, ei)
        dead = Set([other.node, di])

        # Term A: 2 dots(minc) + a pure trivalent(majc)
        gA = CircularGraph(g.word,
                     vcat(g.nodes, [circular_node([majc, majc, majc]), circular_node([minc]), circular_node([minc])]),
                     g.edges)
        triA = length(gA.nodes) - 2; dA1 = length(gA.nodes) - 1; dA2 = length(gA.nodes)
        keepA = Edge[e for (j, e) in enumerate(gA.edges) if !(j in dropbase)]
        for (i, p) in enumerate(maj_leaves); push!(keepA, Edge(majc, p, NodePort(triA, i))); end
        push!(keepA, Edge(minc, min_leaves[1], NodePort(dA1, 1)))
        push!(keepA, Edge(minc, min_leaves[2], NodePort(dA2, 1)))
        termA = _circular_delete_nodes(gA, dead, keepA)

        # Term B: connect the two ADJACENT majc leaves, connect the two minc
        # leaves, dot(majc) on the remaining majc leaf.
        adj1, adj2 = adjpair
        lone = first(p for p in maj_leaves if p != adj1 && p != adj2)
        gB = CircularGraph(g.word, vcat(g.nodes, [circular_node([majc])]), g.edges)
        dB = length(gB.nodes)
        keepB = Edge[e for (j, e) in enumerate(gB.edges) if !(j in dropbase)]
        push!(keepB, Edge(majc, adj1, adj2))
        push!(keepB, Edge(minc, min_leaves[1], min_leaves[2]))
        push!(keepB, Edge(majc, lone, NodePort(dB, 1)))
        termB = _circular_delete_nodes(gB, dead, keepB)

        return CircularComboR(termA) + CircularComboR(termB)
    end
    return nothing
end

# ---- C7: _fr_braid_back (R2-3/R5) --------------------------------------------
"""
    _fr_braid_back(g::CircularGraph) -> Union{Nothing, CircularComboR}

R5: two
`:braid` nodes (m=3, same colour pair) connected back-to-back by exactly
`m=3` internal edges. Rewrite to a SUM (term `id`: three through-strands;
term `pinch`: two trivalents + two dots). `nothing` if no match.
"""

# Order of the two triC arms of a braid when it is replaced by a trivalent
# (arm 3 = the merged inner sector). The order must follow the original
# counter-clockwise ordering: arm 1 lies directly BEHIND the inner sector
# (CCW), arm 2 directly BEFORE it. Plain sorting by slot number is wrong when
# the inner sector doesn't sit at the slot-number boundary.
function _circular_triC_trivalent_order(node, internal_slots::Set{Int}, triC::Int)
    d = _circular_degree(node)
    outer_slots = Set(s for s in 1:d if s ∉ internal_slots)
    tri = [s for s in outer_slots if arm_colour(node, s) == triC]
    length(tri) == 2 || return tri
    start = first(s for s in internal_slots if mod1(s - 1, d) ∉ internal_slots)
    arm2 = mod1(start - 1, d)          # directly before the inner sector
    arm1 = only(s for s in tri if s != arm2)
    return [arm1, arm2]
end

function _fr_braid_back(g::CircularGraph)
    # `arm_count == 6` EXPLICITLY: `kind === :braid` alone covers both the
    # 6-armed ordinary braid and the `2k >= 8`-armed general one, and C7 is a
    # six-arm rule (m = 3), so the arm count is part of the guard.
    braids = [ni for (ni, nd) in enumerate(g.nodes) if nd.kind === :braid && arm_count(nd) == 6]
    for a in braids, b in braids
        a < b || continue
        r = _fr_braid_back_at(g, a, b)
        r === nothing || return r
    end
    return nothing
end

"""
    _fr_braid_back_at(g::CircularGraph, a::Int, b::Int) -> Union{Nothing, CircularComboR}

C7 at ONE given braid pair instead of first-match: same body as
`_fr_braid_back`, just without the pair search — same pattern as
`_fr_braid_relation_at` (C8, below). `nothing` if the pattern doesn't sit at
exactly this pair. This split exists for the region matcher
(circular/rules/CircularRegionRules.jl), which finds candidates `(a, b)` geometrically
via the bigon region between the braids and then calls this body — no second
copy of the id/pinch surgery.
"""
function _fr_braid_back_at(g::CircularGraph, a::Int, b::Int)
    nn = length(g.nodes)
    (1 <= a <= nn && 1 <= b <= nn && a != b) || return nothing
    na, nb = g.nodes[a], g.nodes[b]
    (na.kind === :braid && nb.kind === :braid) || return nothing
    # Six arms explicitly — see `_fr_braid_back` above.
    (arm_count(na) == 6 && arm_count(nb) == 6) || return nothing
    begin
        Set(arms(na)[1:2]) == Set(arms(nb)[1:2]) || return nothing
        m = 3
        internal = [(ei, pa, pb) for (ei, pa, pb) in _circular_edges_at_node(g, a)
                    if pb isa NodePort && pb.node == b]
        length(internal) == m || return nothing
        s, t = arms(na)[1], arms(na)[2]
        outa = sort([(ei, pa, pb) for (ei, pa, pb) in _circular_edges_at_node(g, a)
                     if !(pb isa NodePort && pb.node == b)], by = x -> x[2].slot)
        outb = sort([(ei, pb2, po) for (ei, pb2, po) in _circular_edges_at_node(g, b)
                     if !(po isa NodePort && po.node == a)], by = x -> x[2].slot)
        length(outa) == m && length(outb) == m || return nothing

        internal_bslot = Dict{Int,Int}()
        for (_, pa, pb) in internal
            internal_bslot[pa.slot] = pb.slot
        end
        aslots = [pa.slot for (_, pa, _) in outa]
        bslot_for_aslot = Dict{Int,Int}()
        for k in aslots
            a_int_slot = circular_opposite_slot(na, k)
            haskey(internal_bslot, a_int_slot) || continue
            b_int_slot = internal_bslot[a_int_slot]
            bslot_for_aslot[k] = circular_opposite_slot(nb, b_int_slot)
        end
        length(bslot_for_aslot) == m || return nothing
        outa_port = Dict(pa.slot => po for (_, pa, po) in outa)
        outb_port = Dict(pb2.slot => po for (_, pb2, po) in outb)

        # ---- Term id
        deadid = Set([a, b])
        keepid = Edge[e for (j, e) in enumerate(g.edges)
                      if !(j in Set(vcat([x[1] for x in internal],
                                         [x[1] for x in outa], [x[1] for x in outb])))]
        for k in aslots
            poutA = outa_port[k]
            poutB = outb_port[bslot_for_aslot[k]]
            push!(keepid, Edge(arm_colour(na, k), poutA, poutB))
        end
        termid = _circular_delete_nodes(g, deadid, keepid)

        # ---- Term pinch
        colours_a = [arm_colour(na, pa.slot) for (_, pa, _) in outa]
        triC = count(==(s), colours_a) == 2 ? s : t
        dotC = triC == s ? t : s
        internal_A_slots = Set(pa.slot for (_, pa, _) in internal)
        internal_B_slots = Set(pb.slot for (_, _, pb) in internal)
        triA_aslots = _circular_triC_trivalent_order(na, internal_A_slots, triC)
        triB_aslots = _circular_triC_trivalent_order(nb, internal_B_slots, triC)
        dotA_aslots = [k for k in aslots if arm_colour(na, k) == dotC]
        (length(triA_aslots) == 2 && length(triB_aslots) == 2 &&
         length(dotA_aslots) == 1) || return nothing
        legsA_tri = [outa_port[k] for k in triA_aslots]
        legsA_dot = [outa_port[k] for k in dotA_aslots]
        # Order of the partners at B likewise follows B's geometry, but
        # relative to the inner sector (not raw slot number) — this rules out
        # the case where the inner sector sits mid-range in the slot numbers
        # and sorting by slot number would flip the CCW order.
        legsB_tri_nested = [outb_port[s] for s in triB_aslots]
        legsB_dot_nested = [outb_port[bslot_for_aslot[k]] for k in dotA_aslots]
        gp = CircularGraph(g.word,
                     vcat(g.nodes,
                          [circular_node([triC, triC, triC]), circular_node([triC, triC, triC]),
                           circular_node([dotC]), circular_node([dotC])]),
                     g.edges)
        triAn = length(gp.nodes) - 3; triBn = length(gp.nodes) - 2
        dotAn = length(gp.nodes) - 1; dotBn = length(gp.nodes)
        deadp = Set([a, b])
        dropp = Set(vcat([x[1] for x in internal], [x[1] for x in outa], [x[1] for x in outb]))
        keepp = Edge[e for (j, e) in enumerate(gp.edges) if !(j in dropp)]
        push!(keepp, Edge(triC, legsA_tri[1], NodePort(triAn, 1)))
        push!(keepp, Edge(triC, legsA_tri[2], NodePort(triAn, 2)))
        push!(keepp, Edge(triC, legsB_tri_nested[1], NodePort(triBn, 1)))
        push!(keepp, Edge(triC, legsB_tri_nested[2], NodePort(triBn, 2)))
        push!(keepp, Edge(triC, NodePort(triAn, 3), NodePort(triBn, 3)))
        push!(keepp, Edge(dotC, legsA_dot[1], NodePort(dotAn, 1)))
        push!(keepp, Edge(dotC, legsB_dot_nested[1], NodePort(dotBn, 1)))
        termpinch = _circular_delete_nodes(gp, deadp, keepp)

        return CircularComboR(termid) + CircularComboR(termpinch)
    end
    return nothing
end

# ---- C8: _fr_braid_relation (R2-4/R9) ----------------------------------------

# Roles of the six slots of a path-A braid, circular version of `_r9_braid_roles`.
function _fr_braid_roles(g::CircularGraph, X::Int, Y::Int, tv::Int)
    slot = Dict{Int,Port}()
    for (_, pa, po) in _circular_edges_at_node(g, X)
        pa isa NodePort || return nothing
        haskey(slot, pa.slot) && return nothing
        slot[pa.slot] = po
    end
    length(slot) == 6 || return nothing
    tvs = [k for k in 1:6 if slot[k] isa NodePort && slot[k].node == tv]
    ys  = [k for k in 1:6 if slot[k] isa NodePort && slot[k].node == Y]
    (length(tvs) == 1 && length(ys) == 2) || return nothing
    a = tvs[1]
    for d in (1, -1)
        (mod1(a + d, 6) in ys && mod1(a + 2d, 6) in ys) || continue
        return (d, mod1(a - d, 6), mod1(a - 2d, 6), mod1(a - 3d, 6),
                mod1(a + d, 6), mod1(a + 2d, 6), slot)
    end
    return nothing
end

"""
    _fr_braid_relation(g::CircularGraph) -> Union{Nothing, CircularComboR}

R9: two `:braid`
nodes (6-armed, same colour pair) + a pure trivalent(q), wired in the path-A
pattern (trivalent with one leg into each braid, one EXTERNAL connection; the
two braids connected by exactly two direct edges). Rewrite (path B) to ONE
`:braid`(q,p) + ONE pure trivalent(p), node count 3 ↦ 2. `nothing` if no match.

The result uses LOOSE LEAVES, matching the plain version: roles are read off
the wiring STRUCTURE, not off boundary positions. The third node must be a
PURE trivalent — the rewrite is stated for that shape, not for an arbitrary
circular node.
"""
function _fr_braid_relation(g::CircularGraph)
    # `arm_count == 6` explicitly, same reason as in C7.
    braids = [ni for (ni, nd) in enumerate(g.nodes) if nd.kind === :braid && arm_count(nd) == 6]
    trivs  = [ni for (ni, nd) in enumerate(g.nodes) if _circular_is_pure_trivalent(nd)]
    for b1 in braids, b2 in braids, tv in trivs
        b1 >= b2 && continue
        r = _fr_braid_relation_at(g, b1, b2, tv)
        r === nothing || return r
    end
    return nothing
end

"""
    _fr_braid_relation_at(g::CircularGraph, b1::Int, b2::Int, tv::Int)
        -> Union{Nothing, CircularComboR}

C8 at ONE given site instead of first-match: same body as `_fr_braid_relation`,
just without the site search. `nothing` if the pattern doesn't sit at exactly
this triple.

This split exists for the round-trip check of the inverse direction
(`expand_circular_braid_relation_variants`, circular/rules/CircularInverseRules.jl):
that check needs a SPECIFIC site rather than first-match, since first-match
would discard valid preimages whenever the diagram has a SECOND C8 site. The
order of `b1`/`b2` doesn't matter — which braid plays role S and which plays
B is read off the wiring itself.
"""
function _fr_braid_relation_at(g::CircularGraph, b1::Int, b2::Int, tv::Int)
    m = length(g.nodes)
    (1 <= b1 <= m && 1 <= b2 <= m && 1 <= tv <= m) || return nothing
    (b1 != b2 && tv != b1 && tv != b2) || return nothing
    n1, n2, nt = g.nodes[b1], g.nodes[b2], g.nodes[tv]
    (n1.kind === :braid && n2.kind === :braid) || return nothing
    _circular_is_pure_trivalent(nt) || return nothing
    (arm_count(n1) == 6 && arm_count(n2) == 6) || return nothing
    Set(arms(n1)[1:2]) == Set(arms(n2)[1:2]) || return nothing
    q = arms(nt)[1]
    q in arms(n1) || return nothing
    p = q == arms(n1)[1] ? arms(n1)[2] : arms(n1)[1]
    tv_edges = _circular_edges_at_node(g, tv)
    length(tv_edges) == 3 || return nothing
    to_1 = [po for (_, _, po) in tv_edges if po isa NodePort && po.node == b1]
    to_2 = [po for (_, _, po) in tv_edges if po isa NodePort && po.node == b2]
    t_ext = [po for (_, _, po) in tv_edges
             if !(po isa NodePort && (po.node == b1 || po.node == b2))]
    (length(to_1) == 1 && length(to_2) == 1 && length(t_ext) == 1) || return nothing

    r1 = _fr_braid_roles(g, b1, b2, tv)
    r2 = _fr_braid_roles(g, b2, b1, tv)
    (r1 === nothing || r2 === nothing) && return nothing
    bS, bB, rS, rB = r1[1] == 1 && r2[1] == -1 ? (b1, b2, r1, r2) :
                     r2[1] == 1 && r1[1] == -1 ? (b2, b1, r2, r1) : (0, 0, r1, r2)
    bS == 0 && return nothing
    sS, sB = rS[7], rB[7]
    (sS[rS[5]] isa NodePort && sS[rS[5]].node == bB && sS[rS[5]].slot == rB[5]) || return nothing
    (sS[rS[6]] isa NodePort && sS[rS[6]].node == bB && sS[rS[6]].slot == rB[6]) || return nothing
    arm_colour(g.nodes[bS], rS[5]) == p || return nothing
    arm_colour(g.nodes[bB], rB[5]) == p || return nothing

    S_far, S_mid, S_near = sS[rS[4]], sS[rS[3]], sS[rS[2]]
    B_far, B_mid, B_near = sB[rB[4]], sB[rB[3]], sB[rB[2]]
    e_T = t_ext[1]

    gplus = CircularGraph(g.word,
                     vcat(g.nodes, [circular_node([q, p, q, p, q, p]),
                                    circular_node([p, p, p])]),
                     g.edges)
    NB, NT = length(gplus.nodes) - 1, length(gplus.nodes)
    dead = Set([bS, bB, tv])
    drop = Set{Int}()
    for ni in (bS, bB, tv), (ei, _, _) in _circular_edges_at_node(g, ni)
        push!(drop, ei)
    end
    keep = Edge[e for (j, e) in enumerate(gplus.edges) if !(j in drop)]
    push!(keep, Edge(q, S_mid,  NodePort(NB, 1)))
    push!(keep, Edge(p, S_near, NodePort(NB, 2)))
    push!(keep, Edge(q, e_T,    NodePort(NB, 3)))
    push!(keep, Edge(p, B_near, NodePort(NB, 4)))
    push!(keep, Edge(q, B_mid,  NodePort(NB, 5)))
    push!(keep, Edge(p, NodePort(NB, 6), NodePort(NT, 2)))
    push!(keep, Edge(p, S_far,  NodePort(NT, 1)))
    push!(keep, Edge(p, B_far,  NodePort(NT, 3)))
    return CircularComboR(_circular_delete_nodes(gplus, dead, keep))
end
