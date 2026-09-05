# circular/rules/CircularD4Node.jl — the D4 NODE variant.
#
# Defines apply_circular_d4_node/find_circular_d4_node_match plus their
# _circular_d4_node_*-helpers. Needs CircularDecorated/CircularGraph
# (circular/CircularDecoratedRules.jl) — included right after it.
#
# ---- D4 NODE VARIANT -----------------------------------------------------
#
# WHY IT'S NEEDED. The existing D4 rule above looks, next to the dot, for an
# EDGE of the same colour that borders a region of smaller distance. Because
# of the merge lock in C5, 1/3-nodes remain standing where such an edge would
# otherwise form — e.g. for `id₁ ⊗ id₃ ⊗ dot₁`, the edge variant does not
# fire at all.
#
# THE RULE:
#
#   For dot 1 or dot 3: is there a circular 13-node in the region such that it
#   has 2 strands of the same colour bordering regions of smaller distance? If
#   so, a D4 variant should apply here. To do this, fix 2 strand-arms of the
#   node that share a colour and such that the region lying between them has
#   minimal distance.
#
#   The current diagram is in the dot_A case. dot_B would mean: insert an arm
#   at the 13-node, in the same region as the boundary dot, and connect the
#   node to the boundary. At the same time, one of the other two arms is
#   removed and turned into a dot at its neighbour. dot_C works the same way
#   with the other arm.
#
#   And for the trivalent part, one simply connects and writes a `-alpha_s`
#   into the region bounded by both arms.
#
# I.e. the same relation `dot_A = dot_B + dot_C - alpha_s · (just connected)` as
# above, except the role of the "connecting edge" is played here by a PAIR of
# same-coloured ARMS of one node.
#
# "THE REGION BETWEEN THE TWO ARMS": for a 1-dot with a 13-node, if a 1
# borders the same region as the boundary dot, D4 applies there directly. If
# not, the node is `1133`, bordered by 3s, and the region lies between the
# two 1-arms instead. Since `3·alpha_1 = alpha_1` it genuinely doesn't matter
# which region the alpha goes into. Only important that it is BETWEEN the
# 1-arms — and as far away from the dot as possible.
#
# This sharpens the division of labour between the two variants:
#   * an edge of the DOT COLOUR borders the dot region ⇒ that's the EDGE
#     variant's case (`find_circular_d4_match`);
#   * the dot region is instead bounded there by arms of the OTHER colour
#     (`1133`, `33131`, …) ⇒ this NODE variant, and `alpha_s` goes into the arc
#     BETWEEN the two same-coloured arms, i.e. the arc that does NOT contain
#     the dot region.
#
# This is exactly why the alternating crossing `[1,3,1,3]` does not match here:
# every sector there borders both a 1 and a 3, so the dot region already has an
# edge of the dot colour on its boundary — that's the edge variant's case.
#
# WHICH region of the far arc gets `alpha_s` is mathematically IRRELEVANT
# (`3·alpha_1 = alpha_1`): the one with smallest distance to the mark is
# chosen, because that region is exactly the one that must also carry the
# "smaller distance than the dot region" condition.
#
# WHERE THE NEW ARM IS INSERTED: the sector of `v` that faces the dot's region
# (`circular_node_sector_regions`, circular/CircularRegion.jl). If it isn't unique, the match
# is SKIPPED rather than guessed — same discipline as `_circular_d4_cyclic_ends`
# above, for the same reason: a wrong slot position produces a non-planar
# node, and that doesn't show up as an exception but as a silently collapsed
# cell structure.

"""
    find_circular_d4_node_match(m::CircularMorphismGraph)
        -> Union{Nothing, NamedTuple}

Looks for a match of the **D4 node variant** (see file header): a circular dot of
colour `s` in region `R`, a node `v` on the boundary of `R`, and two arms `p`,
`q` of `v` of colour `s` (same as the dot). The two arms split
the boundary of `v` into two arcs; `R` lies in one of them, and in the OTHER —
the arc "between the `s`-arms", away from the dot — there must be a region
`R_between` of **strictly smaller** distance than `R`.

ADDITIONALLY: in that far arc there must be **no further arm of colour `s`** —
`p` and `q` must be neighbours in the dot colour. Without this condition,
several equally good pairs would be indistinguishable except by plain loop
order, and the fresh dots could land on the wrong strands.

Returns `(dot_node, v, p, q, insert_after, R, R_between, dist_between)` or
`nothing`. `insert_after` is the slot of `v` AFTER which the new arm is
inserted (the sector `insert_after → insert_after+1` is the dot region `R`).

Dots are searched, as in [`find_circular_d4_match`](@ref), in DESCENDING region
distance so that `Σ_dots dist` decreases; within one dot, `R_between` is chosen
with **minimal** distance. Which region of the far arc it ends up being is
mathematically irrelevant (`3·alpha_1 = alpha_1`) — the minimal one is the one
that must carry the distance condition anyway.

SKIPPED (not thrown) when the situation isn't unique: several sectors of `v`
face `R`, or `R` lies in BOTH arcs. Never guessed.
"""
function find_circular_d4_node_match(m::CircularMorphismGraph)
    g = m.graph
    dist = circular_region_distances(m)
    all(==(-1), dist) && return nothing              # no mark
    sector = circular_node_sector_regions(g)

    # region of each dot, largest distance first (as in find_circular_d4_match)
    dots = Tuple{Int,Int}[]
    for (ni, nd) in enumerate(g.nodes)
        arm_count(nd) == 1 || continue
        push!(dots, (ni, circular_region_of_dot(g, ni)))
    end
    sort!(dots; by = t -> (-dist[t[2]], t[1]))

    for (ni, R) in dots
        dist[R] > 0 || continue                      # nothing to gain at the mark
        s = arm_colour(g.nodes[ni], 1)
        (s == 1 || s == 3) || continue               # "for dot 1 or dot 3"
        best = nothing
        for (v, nd) in enumerate(g.nodes)
            v == ni && continue
            _circular_braidlike(nd) && continue
            d = arm_count(nd)
            d >= 3 || continue
            # the sector facing the dot region — must be UNIQUE
            facing = [j for j in 1:d if sector[v][j] == R]
            length(facing) == 1 || continue
            insert_after = facing[1]
            slots = [j for j in 1:d if arm_colour(nd, j) == s]
            for a in 1:length(slots), b in (a + 1):length(slots)
                p, q = slots[a], slots[b]
                far = _circular_d4_node_far_arc(sector[v], d, p, q, R)
                far === nothing && continue
                # NEIGHBOURHOOD IN THE DOT COLOUR: in the far arc there must be
                # NO further `s`-arm. Otherwise the pair doesn't bound the
                # region the `alpha_s` is written into at all, and the fresh
                # dots land on the wrong strands — several pairs can tie at the
                # same `dist_between`, and loop order alone would pick an
                # arbitrary one. Combined with distance minimisation, this
                # condition leaves exactly one candidate pair.
                # The relation `dot_A = dot_B + dot_C − trivalent(α)` lives on
                # THREE strands of ONE colour; a further `s`-arm in between
                # breaks that picture.
                inner_slots = _circular_d4_node_far_inner_slots(sector[v], d, p, q, R)
                inner_slots === nothing && continue
                any(j -> arm_colour(nd, j) == s, inner_slots) && continue
                for Rb in far
                    Rb >= 1 || continue
                    dist[Rb] >= 0 && dist[Rb] < dist[R] || continue
                    cand = (dot_node = ni, v = v, p = p, q = q,
                            insert_after = insert_after, R = R,
                            R_between = Rb, dist_between = dist[Rb])
                    if best === nothing || cand.dist_between < best.dist_between
                        best = cand
                    end
                end
            end
        end
        best === nothing || return best
    end
    return nothing
end

"""
The sector regions in the arc between the arms `p` and `q` of `v` that does
NOT contain the dot region `R_dot` ("between the s-arms, away from the dot").
`sec` is `circular_node_sector_regions(g)[v]`, `d` the arm count.

`nothing` if `R_dot` lies in both arcs or in neither — then "the dot's side"
isn't determined, and this never guesses.
"""
function _circular_d4_node_far_arc(sec::Vector{Int}, d::Int, p::Int, q::Int, R_dot::Int)
    # sec[j] is the sector AFTER slot j; the arc from p forward to q consists
    # of the sectors after p, p+1, …, q-1.
    arc(from, to) = [sec[mod1(from + k, d)] for k in 0:(mod(to - from, d) - 1)]
    a1, a2 = arc(p, q), arc(q, p)
    in1, in2 = R_dot in a1, R_dot in a2
    in1 == in2 && return nothing        # in both or in neither — not unique
    return in1 ? a2 : a1
end

"""
The SLOTS strictly inside the same far arc that [`_circular_d4_node_far_arc`](@ref)
returns as regions (`p`/`q` themselves don't count). `nothing` under the same
condition as there.

Used for the neighbourhood condition in `find_circular_d4_node_match`: between `p`
and `q`, on the side away from the dot, there must be no FURTHER arm of the
dot colour — see the reasoning there.
"""
function _circular_d4_node_far_inner_slots(sec::Vector{Int}, d::Int, p::Int, q::Int,
                                      R_dot::Int)
    arc(from, to)   = [sec[mod1(from + k, d)] for k in 0:(mod(to - from, d) - 1)]
    inner(from, to) = [mod1(from + k, d)      for k in 1:(mod(to - from, d) - 1)]
    a1, a2 = arc(p, q), arc(q, p)
    in1, in2 = R_dot in a1, R_dot in a2
    in1 == in2 && return nothing
    return in1 ? inner(q, p) : inner(p, q)
end

"""
    apply_circular_d4_node(fd::CircularDecorated, dot_node, v, p, q, insert_after)
        -> Union{Nothing, CircularComboR}

The **D4 node variant** (see file header). The dot `dot_node` (colour `s`) and
its capping edge disappear in ALL three terms; at node `v` a new arm of colour
`s` is inserted AFTER slot `insert_after` and connected to the dot's strand
end, in each term:

1. **dot_B** (`+1`): additionally arm `p` is dropped; its far end gets a fresh
   circular dot (`circular_node([s])`).
2. **dot_C** (`+1`): the same with `q`.
3. **Trivalent term** (`-1`): only connects, nothing is removed; `alpha(s)` is
   placed in a region of the arc BETWEEN `p` and `q` that does not contain the
   dot region ([`_circular_d4_node_between_region`](@ref)). Which region there is
   mathematically irrelevant (`3·alpha_1 = alpha_1`) — the choice is only made
   deterministic so two runs give the same graph. It therefore need not match
   the matcher's `R_between`, which is chosen by smallest DISTANCE (there it
   carries the match condition).

`nothing` if the arguments don't fit (no degree-1 dot, `p`/`q` not of the dot
colour, one of the pair arms leads back to `v` itself).
"""
function apply_circular_d4_node(fd::CircularDecorated, dot_node::Int, v::Int,
                           p::Int, q::Int, insert_after::Int)
    g = fd.graph
    1 <= dot_node <= length(g.nodes) || return nothing
    1 <= v <= length(g.nodes) || return nothing
    dot_node == v && return nothing
    ndot = g.nodes[dot_node]
    arm_count(ndot) == 1 || return nothing
    at = _circular_edges_at_node(g, dot_node)
    length(at) == 1 || return nothing
    (cap_ei, _, dot_strand_end) = at[1]
    s = arm_colour(ndot, 1)

    nd = g.nodes[v]
    d = arm_count(nd)
    (1 <= p <= d && 1 <= q <= d && p != q) || return nothing
    (arm_colour(nd, p) == s && arm_colour(nd, q) == s) || return nothing
    (1 <= insert_after <= d) || return nothing
    # the dot's strand end must not itself hang off `v` (that would create a
    # self-connection, which is `needle`'s business, not D4's)
    dot_strand_end isa NodePort && dot_strand_end.node == v && return nothing

    # far ends of the two pair arms (excluding the capping edge)
    far_p = _circular_far_end_of(g, NodePort(v, p), Set([cap_ei]))
    far_q = _circular_far_end_of(g, NodePort(v, q), Set([cap_ei]))
    (far_p === nothing || far_q === nothing) && return nothing
    (far_p[1] isa NodePort && far_p[1].node == v) && return nothing
    (far_q[1] isa NodePort && far_q[1].node == v) && return nothing

    # ---- the new node: an arm of colour `s` AFTER slot `insert_after` --------
    # `drop` is the slot that disappears in this term (0 = none). Returns:
    # the new node and the slot map old -> new (without `drop`), plus the new
    # slot of the inserted arm.
    function rebuilt_node(drop::Int)
        old_slots = Int[]
        for j in 1:d
            j == drop && continue
            push!(old_slots, j)
            j == insert_after && push!(old_slots, 0)   # 0 = the NEW arm
        end
        # special case: if exactly `insert_after` is dropped, the new arm
        # sits in its place — the sector facing the dot region is the same.
        0 in old_slots || insert_after == drop ||
            error("apply_circular_d4_node: the new arm was not inserted — unexpected")
        if !(0 in old_slots)
            pos = findfirst(==(mod1(insert_after + 1, d)), old_slots)
            pos === nothing ? push!(old_slots, 0) : insert!(old_slots, pos, 0)
        end
        new_arms = [j == 0 ? s : arm_colour(nd, j) for j in old_slots]
        newnd = try
            circular_node(new_arms)
        catch
            _unchecked_circularnode(:mixed, new_arms)
        end
        remap = Dict{Int,Int}(j => k for (k, j) in enumerate(old_slots) if j != 0)
        newslot = findfirst(==(0), old_slots)
        return newnd, remap, newslot
    end

    # Builds the graph of one term: `drop` = the slot dropped at `v`
    # (0 = none), `far` = that arm's far end (its edge is replaced by an edge
    # to a fresh dot).
    #
    # ORDER MATTERS HERE: the EXISTING edges are rewritten to the new slot
    # numbering FIRST, and only THEN is the new edge added. Reversing the
    # order breaks when an existing slot numerically coincides with
    # `newslot` (e.g. `[3,1,3,1,1]`, `insert_after = 3` ⇒ `newslot = 4`,
    # and there is already a slot 4) — the new edge would then become
    # indistinguishable from the existing ones, and the face walk would run past
    # the arm count.
    function build(drop::Int, far)
        newnd, remap, newslot = rebuilt_node(drop)
        nodes2 = CircularNode[i == v ? newnd : g.nodes[i] for i in 1:length(g.nodes)]
        dropped_ei = drop == 0 ? 0 : far[3]

        fixp(x::Leaf) = x
        fixp(x::Circle) = x
        function fixp(x::NodePort)
            x.node == v || return x
            haskey(remap, x.slot) || error(
                "apply_circular_d4_node: slot $(x.slot) at node $v has no image " *
                "— unexpected")
            return NodePort(v, remap[x.slot])
        end

        keep = Edge[]
        for (ei, e) in enumerate(g.edges)
            (ei == cap_ei || ei == dropped_ei) && continue
            push!(keep, Edge(e.colour, fixp(e.a), fixp(e.b)))
        end
        # only now the new edge: the inserted arm to the dot's strand end
        push!(keep, Edge(s, fixp(dot_strand_end), NodePort(v, newslot)))
        if drop != 0
            push!(nodes2, circular_node([s]))
            push!(keep, Edge(far[2], fixp(far[1]), NodePort(length(nodes2), 1)))
        end
        return _circular_delete_nodes(CircularGraph(g.word, nodes2, keep), Set([dot_node]), keep)
    end

    # ---- term 1: just connect, alpha_s into the region between p and q -------
    g1 = build(0, nothing)
    R_between = _circular_d4_node_between_region(g, v, p, q, circular_region_of_dot(g, dot_node))
    between_gaps = regions(g)[R_between].gaps
    isempty(between_gaps) && error(
        "apply_circular_d4_node: the region between the arms has no boundary gaps " *
        "— label transfer impossible (same limitation as _circular_transfer_labels)")
    rn1 = regions(g1)
    between_new = findfirst(S -> between_gaps[1] in S.gaps, rn1)
    between_new === nothing && error(
        "apply_circular_d4_node: gap $(between_gaps[1]) lies in no region of the new graph")
    labels1 = _circular_transfer_labels(g, g1, fd.region_labels;
                                   overrides = Dict(between_new => alpha(s)))
    # carry along the outer label.
    term1 = -one(SoergelPoly) *
            CircularCombo{SoergelPoly}(CircularDecorated(g1, labels1, fd.outer_label))

    # ---- terms 2/3: additionally swap arm p resp. q for a dot ----------------
    R_dot = circular_region_of_dot(g, dot_node)
    isone(fd.region_labels[R_dot]) || error(
        "apply_circular_d4_node: the dot's own region carries a non-trivial label " *
        "— unexpected (the dot should be a pure bimodule cap)")
    function dotted(drop::Int, far)
        g2 = build(drop, far)
        return CircularCombo{SoergelPoly}(
            CircularDecorated(g2, _circular_transfer_labels(g, g2, fd.region_labels),
                         fd.outer_label))
    end
    return term1 + dotted(p, far_p) + dotted(q, far_q)
end

"""
The region the trivalent term's `alpha_s` is written into: a region in the arc
BETWEEN the arms `p` and `q` of `v` that does not contain the dot region
`R_dot` — only important that it lies between the two `s`-arms, and as far
away from the dot as possible.

WHICH one it is is mathematically irrelevant (`3·alpha_1 = alpha_1`);
the one with smallest region id is chosen, so the choice is reproducible.
Throws if the arc isn't determined or stays empty — never guessed.
"""
function _circular_d4_node_between_region(g::CircularGraph, v::Int, p::Int, q::Int, R_dot::Int)
    sector = circular_node_sector_regions(g)
    d = arm_count(g.nodes[v])
    far = _circular_d4_node_far_arc(sector[v], d, p, q, R_dot)
    far === nothing && error(
        "_circular_d4_node_between_region: the dot region $R_dot lies at node $v in both " *
        "or in neither of the arcs between arms $p/$q — the dot's side is " *
        "not determined")
    cands = sort(unique(filter(>=(1), far)))
    isempty(cands) && error(
        "_circular_d4_node_between_region: the arc between arms $p/$q at node $v " *
        "contains no region")
    return cands[1]
end
