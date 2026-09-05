# circular/rules/CircularRules.jl — the Circular rules (the C-series).
#
# Uniform signature:
#   _fr_<name>(g::CircularGraph) -> Union{Nothing, CircularComboR}
# `nothing` = no match; an empty `CircularComboR` = the diagram is 0.
#
# ORDER IS SUBSTANTIVE, not cosmetic:
#   1. _fr_needle        → 0; cheapest, must run first
#   2. _fr_dot_walks_on  → MUST fire before merge (ordering trap, see there)
#   3. _fr_barbell       → MUST fire before merge (degree-1 + degree-1 would
#      otherwise hit merge, which is defined on CircularCombo and needs a
#      CircularDecorated result — see CircularDecorated.jl, Part D)
#   4. _fr_bead          → degree-2 normalization ("bead")
#   5. _fr_mono_double   → a SAME-coloured double connection ↦ 0. MUST fire
#      before C12: C12 only matches DIFFERENT colours,
#      so it would otherwise never see the same-coloured case — the two
#      rules are complementary at the same spot.
#   6. _fr_two_adjacent_merge (C12) → MUST fire before C10 (ordering trap, see
#      there). Only DIFFERENT-coloured double connections (see
#      _fr_mono_double).
#   7. _fr_commutation_merge (C10) → the adjacent-or-separated 1-3 pair
#   8. _fr_merge         → the workhorse
#   9. the three two-colour rules (C6–C8)
#
# The MERGE family — C5 `_fr_merge`, the merge lock for mixed 1/3 nodes, Cxx
# `_fr_mono_double`, C12 `_fr_two_adjacent_merge` and C10
# `_fr_commutation_merge` — lives in circular/rules/CircularMergeRules.jl,
# included directly after this file. What stays here are the helpers and the
# local rules C1–C4 and C17.

# ---- Helpers -----------------------------------------------------------------

"The edges of `g` at node `ni` (as (edge index, the port at ni, the other port))."
function _circular_edges_at_node(g::CircularGraph, ni::Int)
    hits = Tuple{Int,Port,Port}[]
    for (ei, e) in enumerate(g.edges)
        if e.a isa NodePort && e.a.node == ni
            push!(hits, (ei, e.a, e.b))
        elseif e.b isa NodePort && e.b.node == ni
            push!(hits, (ei, e.b, e.a))
        end
    end
    return hits
end

"""
    _circular_slot_ports(g::CircularGraph, ni::Int) -> Dict{Int,Port}

The OCCUPIED SLOTS of node `ni`, as `slot => opposite port`.

Why not `_circular_edges_at_node`: its `elseif` above yields exactly ONE entry per
edge. A SELF-LOOP (both ends at `ni`) occupies TWO slots but would appear only
once there — counting occupied slots from it undercounts. Here each end is
recorded separately, whether it is the `a`- or `b`-end, and whether both ends
hang off the same node.

For a self-loop `slot p ↔ slot q`, `d[p] = NodePort(ni, q)` and
`d[q] = NodePort(ni, p)` — the opposite port points back to `ni` itself.
"""
function _circular_slot_ports(g::CircularGraph, ni::Int)
    d = Dict{Int,Port}()
    for e in g.edges
        e.a isa NodePort && e.a.node == ni && (d[e.a.slot] = e.b)
        e.b isa NodePort && e.b.node == ni && (d[e.b.slot] = e.a)
    end
    return d
end

"""
The connecting edges between the DIFFERENT nodes `ni` and `nj`, as
`slot(ni) -> (edge index, slot(nj))`. Self-connections (`ni == nj`) and edges
with a leaf/circle end are excluded.

Used by `_fr_two_adjacent_merge` (C12), `_circular_adjacent_conn_pairs` (C10), and
`_fr_merge` (C5, for counting).
"""
function _circular_conn_slots(g::CircularGraph, ni::Int, nj::Int)
    at = Dict{Int,Tuple{Int,Int}}()
    ni == nj && return at
    for (ei, e) in enumerate(g.edges)
        a, b = e.a, e.b
        a isa NodePort && b isa NodePort || continue
        if a.node == ni && b.node == nj
            at[a.slot] = (ei, b.slot)
        elseif b.node == ni && a.node == nj
            at[b.slot] = (ei, a.slot)
        end
    end
    return at
end

"The number of connecting edges between the different nodes `ni` and `nj`."
_circular_conn_count(g::CircularGraph, ni::Int, nj::Int) = length(_circular_conn_slots(g, ni, nj))

"""
Does `nd` carry BOTH colours of the 1/3 world? That's the criterion for the
merge lock in `_fr_merge` (C5): a node carrying 1 and 3
simultaneously is not monochromatic and must not be merged across a SINGLE
connection with another such node.

BY THE TYPE: this is equivalent to "a pure 1/3 node that carries both
colours". `circular_node` (CircularGraph.jl:125-134) never mixes colour 2 with 1/3 — an
arm sequence containing both 1 AND 3 necessarily lies ENTIRELY in `{1,3}`. A
`[1,3,2,3]` node cannot exist. The phrase "carries both colours" is checked
directly here, matching the stated condition, independent of that type
invariant.
"""
_circular_is_bicoloured(nd::CircularNode) = (1 in arms(nd)) && (3 in arms(nd))

# Rebuilds `g` with the nodes in `dead` removed. Every surviving NodePort is
# rewritten to the new node numbering; `keep_edges` (still using OLD node
# indices) is rewritten too. Leaves/circles are left unchanged.
function _circular_delete_nodes(g::CircularGraph, dead::Set{Int}, keep_edges::Vector{Edge})
    remap = Dict{Int,Int}()
    newnodes = CircularNode[]
    for (ni, nd) in enumerate(g.nodes)
        ni in dead && continue
        push!(newnodes, nd)
        remap[ni] = length(newnodes)
    end
    fixport(p::Leaf) = p
    fixport(p::Circle) = p
    fixport(p::NodePort) = NodePort(remap[p.node], p.slot)
    newedges = [Edge(e.colour, fixport(e.a), fixport(e.b)) for e in keep_edges]
    return CircularGraph(g.word, newnodes, newedges)
end

# ---- C1: _fr_needle — self-connection ↦ 0 ------------------------------------
#
# Replaces `_r7_free_circle` (free circle) AND `_r10_trivalent_self_loop`
# (trivalent self-loop) with ONE rule: an edge with both ports on the SAME
# node (any slots), node NOT `:braid` ⇒ empty `CircularComboR`. A free `Circle`
# port in the same body means the same thing ("a closed monochromatic loop
# enclosing nothing").
#
# A self-connection at a `:braid` node, including the 2k-armed gbraid, IS 0
# too — C1 is responsible for ALL self-connections, braid-like or not.
"""
    _fr_needle(g::CircularGraph) -> Union{Nothing, CircularComboR}

Needle rule: a free circle (`Edge(c, Circle(c), Circle(c))`) OR an edge with
both ends at the SAME node ⇒ the diagram is 0 (empty `CircularComboR`). `nothing`
if neither holds. Applies to `:braid` nodes too, including the 2k-armed
gbraid (see above).
"""
function _fr_needle(g::CircularGraph)
    for e in g.edges
        if e.a isa Circle || e.b isa Circle
            return CircularComboR()
        end
        if e.a isa NodePort && e.b isa NodePort && e.a.node == e.b.node
            return CircularComboR()
        end
    end
    return nothing
end

# ---- C21: _fr_pitchfork — a dot between two same-coloured arms ↦ 0 ------------
#
# The same statement as C1, for the general figure: C1 needs a genuine self-edge
# (`e.a.node == e.b.node`), while here the two same-coloured neighbour arms may
# run to ANOTHER node and meet there. Whether they close on their own node (then
# it is the self-loop, and `:needle` comes first in `CIRCULAR_RULES` anyway) or
# on a foreign one makes no difference to the conclusion: the term is 0.
#
# The figure this covers is a dot capping one arm of a many-armed node while its
# two neighbours, of one colour, both land on a single `:mono` node — the shape
# that survives rex fusion on `121321/121321` and that C1 does not see.
#
# The smallest instance is the reduction `121 → 212 → 22 → ε` drawn out: a
# 6-armed braid, one 2-arm capped by a dot, its two neighbouring 1-arms meeting at
# one trivalent (`test/circularregionrules.jl`).
"""
    _fr_pitchfork(g::CircularGraph) -> Union{Nothing, CircularComboR}

The general pitchfork rule: at a node `v` with at least 3 arms a **dot** caps
the arm `s`, and the two neighbouring arms `s−1`/`s+1` have the **same colour**
and end at the **same** node `w` (at two different ports) ⇒ the term is **0**
(the empty `CircularComboR`). Otherwise `nothing`.

Slots are read cyclically (`mod1`) and the arm count is arbitrary — the rule has
no special case for `:braid`/`:gbraid`.
"""
function _fr_pitchfork(g::CircularGraph)
    ends = Dict{Tuple{Int,Int},Any}()
    for e in g.edges, (p, q) in ((e.a, e.b), (e.b, e.a))
        p isa NodePort || continue
        ends[(p.node, p.slot)] = q
    end
    for (v, nd) in enumerate(g.nodes)
        d = arm_count(nd)
        d >= 3 || continue
        for s in 1:d
            dot = get(ends, (v, s), nothing)
            dot isa NodePort || continue
            arm_count(g.nodes[dot.node]) == 1 || continue         # the dot
            sm, sp = mod1(s - 1, d), mod1(s + 1, d)
            a = get(ends, (v, sm), nothing)
            b = get(ends, (v, sp), nothing)
            (a isa NodePort && b isa NodePort) || continue
            a.node == b.node || continue                          # they meet
            (a.node, a.slot) != (b.node, b.slot) || continue       # at two ports
            nd.arms[sm] == nd.arms[sp] || continue                 # same colour
            return CircularComboR()
        end
    end
    return nothing
end

# ---- C2: _fr_dot_walks_on — case distinction ---------------------------------
#
# Match: node `nd` (`:mixed`/`:mono`), edge from slot `s` to a degree-1 node
# of colour `c = arms(nd)[s]`, AND `count(==(c), arms(nd)) == 2`. Let `s'` be
# the OTHER slot of this colour. Result (one term, coefficient 1): new node =
# `arms(nd)` without both slots `s` and `s'`, read cyclically; a FRESH dot
# `circular_node([c])` is attached to the far end of the old `s'` edge.
#
# THE ORDERING TRAP (see the order list above): `[1,3,1,3]` + dot gives
# `[3,1,3]` under merge, but the spec requires a PURE 3-node for the two-arm
# case. `_fr_dot_walks_on` MUST therefore fire before `_fr_merge`.
#
# THE `[1,1,3,3]` CASE ("a dot on 1133 gives a through line and the
# dot keeps moving"): the same rule, except the two slots of colour c are
# ADJACENT rather than opposite. Result: the two 1-arms become a continuous
# 1-line, the 3-dot slides to the far arm of the other 3-slot.
"""
    _fr_dot_walks_on(g::CircularGraph) -> Union{Nothing, CircularComboR}

A dot (a degree-1 node of colour `c`) caps slot `s` of a node `nd`, and `c`
occurs in `arms(nd)` exactly TWICE (slots `s` and `s'`). The dot "keeps
moving": the result is a new node without slots `s`/`s'` (read cyclically —
the block after `s` up to before `s'`, then the block after `s'` up to before
`s`), with a fresh dot of colour `c` at the far end of the old `s'` edge.
`arm_count(nd) >= 3` is assumed (with exactly 2 arms of colour `c` nothing
would remain — the `circular_node` constructor would throw anyway). `nothing` if
no such dot exists.

SCOPE: the case where the OTHER `c`-edge (at slot `s'`) is ITSELF capped by a
dot (two dots of the same colour on a node with exactly 2 arms of that
colour) is NOT matched here — `_circular_far_end` would then
return a `NodePort` to another degree-1 node, which would not be
structurally wrong, but the guard below only checks that slot `s` goes to a
dot; slot `s'`'s far end can be anything, INCLUDING a second dot of the same
colour. In that special case a new dot would arise whose edge leads to
another dot — geometrically plausible but NOT separately tested.
"""
function _fr_dot_walks_on(g::CircularGraph)
    for (di, dnd) in enumerate(g.nodes)
        arm_count(dnd) == 1 || continue
        c = arms(dnd)[1]
        at = _circular_edges_at_node(g, di)
        length(at) == 1 || continue                 # a genuine degree-1 dot
        (cap_ei, _dotport, other) = at[1]
        other isa NodePort || continue
        other.node == di && continue                # no self-loop here (C1's job)
        nd = g.nodes[other.node]
        _circular_braidlike(nd) && continue               # braid-like nodes (incl. gbraid) can't be merged/capped
        arm_colour(nd, other.slot) == c || continue   # capping edge must carry colour c
        s = other.slot
        count(==(c), arms(nd)) == 2 || continue       # ONLY the two-arm case belongs here
        deg = arm_count(nd)
        # find the other slot of the same colour c (s' ≠ s).
        sprime = findfirst(k -> k != s && arm_colour(nd, k) == c, 1:deg)
        sprime === nothing && continue

        at_nd = _circular_edges_at_node(g, other.node)
        sp_entry = findfirst(e -> e[2].slot == sprime, at_nd)
        sp_entry === nothing && continue              # s' not wired — no match
        (sp_ei, _sp_port, sp_far) = at_nd[sp_entry]

        # the remaining slots in order: the block after s up to before s',
        # then the block after s' up to before s (cyclic, both boundary slots
        # s/s' excluded). Both the arm sequence AND the slot remap come from
        # THESE — the two blocks may have any length, including zero.
        rem_slots = Int[]
        k = mod1(s + 1, deg)
        while k != sprime
            push!(rem_slots, k)
            k = mod1(k + 1, deg)
        end
        k = mod1(sprime + 1, deg)
        while k != s
            push!(rem_slots, k)
            k = mod1(k + 1, deg)
        end
        new_arms = Int[arm_colour(nd, k) for k in rem_slots]
        isempty(new_arms) && continue                 # nothing left — circular_node would throw

        if length(new_arms) == 2
            # Exactly 2 arms remain — this WOULD be a bead, but beads are NOT
            # valid `CircularNode`s ("no degree-2 node" invariant). The
            # remainder collapses IMMEDIATELY instead (as `_fr_bead` would):
            # the two remaining neighbours are connected directly, no
            # intermediate node.
            #
            # The restover slots come from `rem_slots`, the same source as
            # the arm sequence — this correctly handles both the layout with
            # one restover arm per side (e.g. `[1,3,1,3]`) and the layout
            # with one side empty and both restover arms on the other (e.g.
            # `[1,1,3,3]` with the dot on a 3-arm).
            #
            # Both restover arms NECESSARILY carry the same colour (the
            # other one): the node has four arms, two of colour c, and a
            # `:mixed`/`:mono` node only knows two colours. `arm_colour(nd,
            # slot1)` is thus the colour of the through line.
            slot1, slot2 = rem_slots[1], rem_slots[2]
            far1_at = findfirst(e -> e[1] != cap_ei && e[1] != sp_ei &&
                                     (e[2] isa NodePort && e[2].node == other.node && e[2].slot == slot1),
                                 at_nd)
            far2_at = findfirst(e -> e[1] != cap_ei && e[1] != sp_ei &&
                                     (e[2] isa NodePort && e[2].node == other.node && e[2].slot == slot2),
                                 at_nd)
            (far1_at === nothing || far2_at === nothing) && continue
            (_, _p1, far1) = at_nd[far1_at]
            (_, _p2, far2) = at_nd[far2_at]

            gplus = CircularGraph(g.word, vcat(g.nodes, [circular_node([c])]), g.edges)
            freshdot = length(gplus.nodes)
            drop = Set([e for (e, _, _) in at_nd])
            keep = Edge[e for (j, e) in enumerate(gplus.edges) if !(j in drop)]
            push!(keep, Edge(arm_colour(nd, slot1), far1, far2))
            push!(keep, Edge(c, sp_far, NodePort(freshdot, 1)))
            dead = Set([di, other.node])
            return CircularComboR(_circular_delete_nodes(gplus, dead, keep))
        end

        new_node = circular_node(new_arms)

        # Slot remap: old slots (≠ s, s') → new slots in the order above.
        remap = Dict{Int,Int}(oldslot => pos for (pos, oldslot) in enumerate(rem_slots))

        gplus = CircularGraph(g.word, vcat(g.nodes, [new_node, circular_node([c])]), g.edges)
        newid = length(gplus.nodes) - 1
        freshdot = length(gplus.nodes)
        # ALL edges at other.node are dropped (not just cap/sp) — they are
        # rewired individually below; otherwise stale NodePorts on the
        # removed node `other.node` would remain.
        drop = Set([e for (e, _, _) in at_nd])
        keep = Edge[e for (j, e) in enumerate(gplus.edges) if !(j in drop)]
        for (oldslot, newslot) in remap
            far_at = findfirst(e -> e[1] != cap_ei && e[1] != sp_ei &&
                                    ((e[2] isa NodePort && e[2].node == other.node && e[2].slot == oldslot)),
                                at_nd)
            far_at === nothing && continue
            (_, _p, far) = at_nd[far_at]
            push!(keep, Edge(arm_colour(nd, oldslot), NodePort(newid, newslot), far))
        end
        # the fresh dot goes to the far end of the s' edge
        push!(keep, Edge(c, sp_far, NodePort(freshdot, 1)))
        dead = Set([di, other.node])
        return CircularComboR(_circular_delete_nodes(gplus, dead, keep))
    end
    return nothing
end

# ---- C3: _fr_barbell — factor α_i ---------------------------------------------
#
# An edge between two degree-1 nodes of the SAME colour `i` ⇒ the edge and
# both nodes disappear, the cell gets decoration α_i. MUST fire before
# `_fr_merge` (degree-1 + degree-1 would otherwise throw at `merge_at_edge`,
# which by definition works on CircularCombo and needs a CircularDecorated result — see
# CircularDecorated.jl, Part D).
"""
    _fr_barbell(fd::CircularDecorated) -> Union{Nothing, CircularComboR}

Two degree-1 nodes of the same colour `i`, joined by exactly one edge of
colour `i` (a "barbell"). Both nodes and the edge disappear; the **region**
the barbell lies in gets the decoration α_i. Template: `apply_barbell`. If
the barbell lies in a region WITHOUT a boundary gap (free circle, lens),
`α_i` goes into the `outer_label` of the result instead of the coefficient.
`nothing` if no barbell matches.

**Regions, not faces.** A barbell consists of two dots, and dot edges never
separate regions by definition — so a barbell always lies entirely in ONE
region, which `circular_region_of_dot` returns directly. No separate
precondition is required.
"""
function _fr_barbell(fd::CircularDecorated)
    g = fd.graph
    for (ei, e) in enumerate(g.edges)
        e.a isa NodePort && e.b isa NodePort || continue
        e.a.node == e.b.node && continue
        na, nb = g.nodes[e.a.node], g.nodes[e.b.node]
        arm_count(na) == 1 && arm_count(nb) == 1 || continue
        i = arms(na)[1]
        i == arms(nb)[1] == e.colour || continue
        R = circular_region_of_dot(g, e.a.node)
        dead = Set([e.a.node, e.b.node])
        keep = [ee for (j, ee) in enumerate(g.edges) if j != ei]
        g2 = _circular_delete_nodes(g, dead, keep)
        gaps = regions(g)[R].gaps
        rn = regions(g2)
        nid = isempty(gaps) ? nothing : findfirst(S -> gaps[1] in S.gaps, rn)
        if nid !== nothing
            labels = _circular_transfer_labels(g, g2, fd.region_labels;
                overrides = Dict(nid => fd.region_labels[R] * alpha(i)))
            return CircularComboR(CircularDecorated(g2, labels, fd.outer_label))
        end
        # Region without a boundary gap (free circle, lens): can't be
        # assigned via the boundary — α_i goes into the OUTER region instead.
        # A label on this region goes the same way: `outer`
        # absorbs it instead of throwing. The Ref is seeded with the CURRENT
        # outer label so it isn't lost.
        outer = Ref(fd.outer_label * alpha(i))
        labels = _circular_transfer_labels(g, g2, fd.region_labels; outer = outer)
        return CircularComboR(CircularDecorated(g2, labels, outer[]))
    end
    return nothing
end

# ---- C4: _fr_bead — degree-2 normalization ("bead") ---------------------------
#
# A node with `arm_count == 2` (e.g. `[1,1]`, `circular_degree == 0`) disappears,
# its two edges get spliced into one. Without this rule, reduced diagrams
# would carry such beads.
"""
    _fr_bead(g::CircularGraph) -> Union{Nothing, CircularComboR}

A node with exactly TWO arms ("bead") disappears; its two edges are spliced
into one continuous edge of the same colour. `nothing` if no such node
exists.
"""
function _fr_bead(g::CircularGraph)
    for (ni, nd) in enumerate(g.nodes)
        arm_count(nd) == 2 || continue
        at = _circular_edges_at_node(g, ni)
        length(at) == 2 || continue
        (e1, _, p1) = at[1]
        (e2, _, p2) = at[2]
        # A self-loop at a bead is C1's (_fr_needle's) business, not this one.
        (p1 isa NodePort && p1.node == ni) && continue
        (p2 isa NodePort && p2.node == ni) && continue
        c = arm_colour(nd, 1)
        drop = Set([e1, e2])
        keep = Edge[e for (j, e) in enumerate(g.edges) if !(j in drop)]
        push!(keep, Edge(c, p1, p2))
        return CircularComboR(_circular_delete_nodes(g, Set([ni]), keep))
    end
    return nothing
end

# ---- shared union–find ---------------------------------------------------------

"""
    _circular_union_find(n::Int) -> (root, unite!)

A union–find over `1:n`, returned as the pair of closures the rules use:
`root(x)` gives the representative (path halving), `unite!(a, b)` joins two
classes and always keeps the SMALLER index as the representative — the callers
rely on that, they read the roots back as stable piece labels.

Shared by `circular_edge_linkage_components` (over edges),
`_circular_cluster2_pieces` and `_fcs_braid_links` (both over nodes).
"""
function _circular_union_find(n::Int)
    parent = collect(1:n)
    function root(x::Int)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    function unite!(a::Int, b::Int)
        ra, rb = root(a), root(b)
        ra == rb && return nothing
        parent[max(ra, rb)] = min(ra, rb)
        return nothing
    end
    return root, unite!
end
