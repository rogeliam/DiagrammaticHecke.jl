# circular/rules/CircularMergeRules.jl — the MERGE family of the circular rules.
#
# Included directly after CircularRules.jl, so the helpers
# (`_circular_edges_at_node`, `_circular_conn_slots`, `_circular_delete_nodes`, …)
# are already defined.
#
# What lives here: C5 `_fr_merge` (the workhorse), the merge lock for mixed 1/3
# nodes, Cxx `_fr_mono_double`, C12 `_fr_two_adjacent_merge` and C10
# `_fr_commutation_merge`.

# ---- C5: _fr_merge — the workhorse --------------------------------------------
#
# First edge with both ports `NodePort`, different nodes, neither `:braid`,
# NOT both degree 1 (that's barbell, C3) ⇒ `merge_at_edge`. In one rule this
# covers: merging, the general branch of the dot rule, and — together with C1
# (needle on the resulting self-connection) — R11.
# ---- Cxx: _fr_mono_double — EXACTLY TWO adjacent, SAME-coloured connections ↦ 0
#
# THE RULE: "if the double connection has different colours,
# i.e. 1+3, then it gets merged. With 2 equal colours it's needle and 0."
#
# This is the EXACT SAME spot as C12 below (two non-braid-like nodes, EXACTLY
# two connecting edges, cyclically adjacent on BOTH nodes, wired in opposite
# orientation) — only the colour condition is inverted: here both connecting
# edges carry the SAME colour, in C12 below DIFFERENT colours.
#
# ⚠ NOT A VALID ARGUMENT: "contract one of the two edges with C5, then the
# other becomes a self-connection ⇒ 0 via C1". This reasoning is INVALID: an
# edge may only be merged if one of its sides carries a single colour. For two
# mixed 1/3 nodes this contraction is not even a legal merge operation to
# begin with, so that intermediate picture never arises. The result 0 is NOT a
# consequence of C5+C1 here, but its own rule (the monochromatic bigon).
#
# WHICH ALSO ANSWERS THE FOLLOW-UP: why isn't the DIFFERENT-coloured double
# connection a needle case, given that contracting one edge would turn the
# other into a self-connection? Because that same contraction is illegal for
# the same reason. The rule set is thus NOT order-dependent at this spot,
# and the registry order (this rule/C12/C10 before C5) is an optimization,
# not a correctness crutch.
#
# In the driver, C5 called
# directly can pick an edge with `_circular_conn_count >= 2`, but doesn't
# fire there because an earlier rule matches first — so the
# order is de facto load-bearing in the driver even though it needn't be
# mathematically.
#
# WHY A SEPARATE RULE, NOT AN EXTENSION OF C1 (`_fr_needle`): C1 matches an
# EDGE with BOTH ends on the SAME node (self-loop, one node, one edge). Here
# there are TWO nodes and TWO edges — structurally the C12/C10 pattern ("two
# nodes, two connections"), not C1's. Placing it next to C12 keeps the
# comparison (only the colour guard differs) visible at a glance.
#
# ORDER: MUST COME BEFORE C12 (`_fr_two_adjacent_merge`) — C12 restricts to
# DIFFERENT colours (see there); the two rules are complementary at the same
# spot, so their relative order is strictly arbitrary; it stands first
# anyway as the "cheaper" of the two checks.
#
# SCOPE: adjacent only. The rule matches the same SPOT as C12, the two
# connections cyclically adjacent on both nodes. A SEPARATED same-coloured
# double connection — the analogue of C10's separated `{1,3}` case — matches
# no rule.
"""
    _fr_mono_double(g::CircularGraph) -> Union{Nothing, CircularComboR}

Two different non-`:braid` nodes `A`, `B`, connected by EXACTLY TWO
connecting edges that sit at cyclically adjacent slots `(k, k+1)` on `A` and
likewise at cyclically adjacent slots on `B`, wired in opposite orientation
(the same SPOT as [`_fr_two_adjacent_merge`](@ref), C12) — AND both
connecting edges carry the SAME colour. Result: the empty `CircularComboR()` (the
diagram is 0) — the monochromatic bigon. `nothing` if no such pair
exists, or if the two connections carry DIFFERENT colours (then C12
applies).

The precondition is `arm_count >= 3` on BOTH nodes (variant B:
"monochromatic-ness can surely be implemented for any arm count") — a
same-coloured double connection is 0 at ANY valid arm count, including the
3-armed node (class (ii) from the `_fr_merge` block comment,
C10/CircularRules.jl:1007). By the `CircularNode` constructor invariant a
node has either EXACTLY 1 arm (a dot) or AT LEAST 3 (`arm_count == 2`, a
"bead", is not a valid end state, `CircularGraph.jl:120`) — with only 1 arm
there aren't two distinct slots for two adjacent connections anyway, so the
condition excludes this case implicitly.

MUST come before C12, see the file header comment.
"""
function _fr_mono_double(g::CircularGraph)
    n = length(g.nodes)
    for i in 1:n, j in 1:n
        i < j || continue
        ndi, ndj = g.nodes[i], g.nodes[j]
        _circular_braidlike(ndi) && continue
        _circular_braidlike(ndj) && continue
        di, dj = arm_count(ndi), arm_count(ndj)
        (di >= 3 && dj >= 3) || continue    # any valid arm count except a dot

        at = _circular_conn_slots(g, i, j)
        length(at) == 2 || continue    # EXACTLY two connections between i and j

        (si1, si2) = sort(collect(keys(at)))
        (ei1, sj1) = at[si1]
        (ei2, sj2) = at[si2]

        # Adjacency on i, adjacency + opposite orientation on j — exactly the
        # same two guards as in C12 below.
        mod1(si1 + 1, di) == si2 || continue
        mod1(sj2 + 1, dj) == sj1 || continue

        # THE colour guard: SAME colour on both connections (the counterpart
        # to C12's "different colour" guard below).
        g.edges[ei1].colour == g.edges[ei2].colour || continue

        return CircularComboR()
    end
    return nothing
end

# ---- C12: _fr_two_adjacent_merge — EXACTLY TWO adjacent, DIFFERENT-coloured
#           connections
#
# THE RULE: for exactly two connections both connecting edges get deleted, so
# `1....(k,k+1)....n` connected to `1.....(l,l+1)....m`
# becomes `1....k-1, l+2.....n.....l-1, k+2....n`.
#
# COLOUR GUARD: "if the double connection has different colours, i.e.
# 1+3, then it gets merged. With 2 equal colours it's needle and 0." The
# same-coloured case is `_fr_mono_double`'s (above, MUST come before this
# rule) and returns 0 there instead of merging.
#
# RELATION TO C10: C10 (`_fr_commutation_merge`) does not require adjacency,
# only that the two connecting edges carry colours `{1,3}`
# (`_circular_adjacent_conn_pairs` sorts adjacent pairs first but also matches
# separated ones). C12 here is STRUCTURALLY defined (adjacency on BOTH
# nodes) plus the "different" colour guard — for `{1,3}`-coloured adjacent
# connections, C10 and C12 give the same result — pure overlap,
# not competition. C12 is placed BEFORE C10 as the more specific rule
# (requires adjacency on BOTH nodes); C10 remains responsible for the
# SEPARATED `{1,3}` case (`double_leaf(12321,…)`-style layouts).
#
# WORKED EXAMPLES (see the `merge_at_edges` docstring, circular/CircularMerge.jl, for
# full verification):
#
#   1. A=[1,1,3,3,3] (n=5), B=[1,1,3,3,3] (m=5), pair at A-slots (1,2), B-slots
#      (1,2), opposite orientation (A.1↔B.2, A.2↔B.1, both colour 1). This
#      example is SAME-coloured ⇒ C12 does not match (colour guard);
#      `_fr_mono_double` matches and returns 0 (empty `CircularComboR`).
#   2. A=[1,3,1,1,3] (n=5), B=[3,1,1,3,1] (m=5), pair at A-slots (1,2) with
#      colours (1,3), wired opposite to B-slots (2,1) (colours (1,3)) —
#      DIFFERENT colours at the two connections, as C10 would also see, but
#      here additionally checked STRUCTURALLY (adjacent on both sides). C12
#      matches (colour guard satisfied). Result `[1,1,3,1,3,1]` (A from slot
#      3: 1,1,3; B from slot 3: 1,3,1).
#   3. Four-armed nodes (n=m=4): A=[1,1,3,3], B=[1,1,3,3], pair at
#      (1,2)/(1,2), opposite orientation. SAME-coloured (both colour 1) ⇒
#      C12 does not match, `_fr_mono_double` returns 0.
#
# OPPOSITE ORIENTATION (required, machine-checked, see the `merge_at_edges`
# docstring in circular/CircularMerge.jl): the pair `(k,k+1)` on `A` must connect to
# `(l+1,l)` on `B` — NOT `(l,l+1)` (same orientation). With unequal colours
# on the two connections (the only case left after the colour guard),
# colour-faithfulness often enforces this on its own; the same-orientation
# pairing would leave a NON-adjacent self-connection after contracting one
# edge — the same trap as the `double_leaf(12321,…)` case.

"""
    _fr_two_adjacent_merge(g::CircularGraph) -> Union{Nothing, CircularComboR}

C12. Two different non-`:braid` nodes `A`, `B` are connected by
EXACTLY TWO connecting edges that sit at cyclically ADJACENT slots `(k,
k+1)` on `A` and likewise at cyclically adjacent slots on `B`, wired in
OPPOSITE orientation (`A`-slot `k` ↔ `B`-slot `l+1`, `A`-slot `k+1` ↔
`B`-slot `l`), AND the two connecting edges carry DIFFERENT colours. Both
connecting edges disappear, the rest merges into ONE node: `A` from `k+2`
forward (`n−2` arms), then `B` from `l+2` forward (`m−2` arms) — see
[`merge_at_edges`](@ref) for the precise rule and examples. `nothing` if no
such pair exists (including when the two connections are NOT adjacent on
both nodes — that's C10's territory — or wired same-orientation instead of
opposite), or if the two connections carry the SAME colour (then
[`_fr_mono_double`](@ref) applies, result 0).

MUST come before C10 (`_fr_commutation_merge`): C10 only sees the `{1,3}`
colour case. MUST come after `_fr_mono_double`, see there.
"""
function _fr_two_adjacent_merge(g::CircularGraph)
    n = length(g.nodes)
    for i in 1:n, j in 1:n
        i < j || continue
        ndi, ndj = g.nodes[i], g.nodes[j]
        _circular_braidlike(ndi) && continue
        _circular_braidlike(ndj) && continue
        di, dj = arm_count(ndi), arm_count(ndj)
        (di >= 4 && dj >= 4) || continue    # after removing 2 slots, >= 2 arms remain

        # all connecting edges between i and j: slot(i) -> (edge index, slot(j))
        at = _circular_conn_slots(g, i, j)
        length(at) == 2 || continue    # EXACTLY two connections between i and j

        (si1, si2) = sort(collect(keys(at)))
        (ei1, sj1) = at[si1]
        (ei2, sj2) = at[si2]

        # Adjacency on i.
        mod1(si1 + 1, di) == si2 || continue

        # Opposite orientation: A-slot si1 (the "first" of the pair) must
        # connect to B-slot sj2, A-slot si2 to B-slot sj1, AND sj1/sj2 must
        # themselves be adjacent on j with sj1 = sj2 + 1 (cyclic) — i.e. the
        # pair on j read as (l, l+1) with l = sj1, l+1 = sj2 would be SAME
        # orientation; the required orientation is (b a): the j-slot
        # belonging to si1 is the LARGER one (l+1), the one belonging to si2
        # the SMALLER (l).
        mod1(sj2 + 1, dj) == sj1 || continue
        l, lnext = sj2, sj1     # l = sj2, l+1 = sj1 = the side belonging to si1

        # COLOUR GUARD: C12 only merges when the two connecting edges
        # carry DIFFERENT colours. Same colour is `_fr_mono_double`'s
        # business (above, MUST come before this rule) and returns 0 there
        # instead of merging.
        g.edges[ei1].colour != g.edges[ei2].colour || continue

        try
            merged = merge_at_edges(g, ei1, ei2)
            return CircularComboR(merged)
        catch err
            err isa ArgumentError || rethrow()
            continue    # e.g. colour-faithfulness violated despite structural adjacency — no match
        end
    end
    return nothing
end

# ---- C10: _fr_commutation_merge — the adjacent 1-3 pair -----
#
# BACKGROUND. In the Circular picture there are only `12`- and `23`-braids;
# a 1-3 crossing is not a `:braid` node, but an ordinary `:mixed` node with
# arms `[1,3,1,3]` (circular/CircularConvert.jl). So R12 ("commutation there
# and back ↦ identity") has no object of its own in the Circular picture —
# the statement has to follow from the 1/3 nodes themselves. The naive
# `_fr_merge` does NOT do that: it contracts exactly ONE edge and leaves
# every further connection as a self-connection, which `needle` then pulls
# to 0 — even in the two-strand case, where R12 demands the identity.
#
# THE RULE: "if there is a 1 and a 3 connection one after
# another, then the rest gets merged and these two edges vanish. Not ALL
# connections." So exactly ONE adjacent connecting pair is consumed (one
# edge colour 1, one edge colour 3, cyclically adjacent arms on BOTH nodes);
# every other connection stays and becomes a self-connection — which is why
# the THREE-strand case is still 0 (`needle` on the leftover edge), but the
# TWO-strand case is the identity.
#
# NO CASE SPLIT BY ARM COUNT:
#
#   "we have 2 13-nodes, say they look like (i1i2......ab) and
#    (baj1....jm), where a and b are connected (going clockwise on one,
#    counter-clockwise on the other, i.e. opposite). When they merge, the
#    new node has the order (i1......inj1....jm) and the i's and j's are
#    still connected exactly as before"
#
# So always the same principle: `a`/`b` disappear, the rest is concatenated
# in its cyclic order, and every remaining edge keeps its far end. Where the
# pair sits in the node is irrelevant (start, end, or middle) — a CircularNode is
# only defined up to rotation anyway, and the matcher searches cyclically.
#
# A `1133`-node arising this way (the case with two arms left on each side)
# is a REGULAR outcome, not a special case: it simply stays as a 4-armed
# node.
#
# MUST FIRE BEFORE `_fr_merge`: otherwise the single-edge contraction there
# creates the self-connection that would wrongly pull the two-strand case to
# 0.

"""
Arm pairs `(s1, s2)` at node `ni` whose edges both go to the same node `nj`
and whose colours are `{1, 3}`.

ADJACENCY IS NOT REQUIRED. The specification is by colour, not by
adjacency:

> a 1133-node falls apart. With ONE connection it gets merged. With TWO
> connections of DIFFERENT colours it also gets merged. With a double
> connection of ONE colour (adjacent) or THREE connections in a row it's 0.

On a 4-armed node both readings coincide (two slots of different colour ARE
adjacent there). On a 5-armed node like `[3,1,3,1,1]` the two connections
can be separated by a third arm — then C10 wouldn't match, `_fr_merge` would
fire its single-edge contraction, leave a self-connection, and `needle`
would wrongly pull the diagram to 0 (found at the DL
`double_leaf(12321,[1,0,1,0,0]²)`).

ADJACENT PAIRS COME FIRST — that's C10's match order, and the adjacent case
is the already-secured one (R12 in the 4-armed picture).
"""
function _circular_adjacent_conn_pairs(g::CircularGraph, ni::Int, nj::Int)
    nd = g.nodes[ni]
    d = arm_count(nd)
    # slot → (edge index, slot on the other node), only for edges ni–nj
    at = _circular_conn_slots(g, ni, nj)
    adj  = Tuple{Int,Int,Int,Int,Int,Int}[]   # (s1,s2, ei1,ei2, t1,t2)
    far  = Tuple{Int,Int,Int,Int,Int,Int}[]
    slots = sort(collect(keys(at)))
    for s in slots, s2 in slots
        s == s2 && continue
        (ei1, t1) = at[s]; (ei2, t2) = at[s2]
        ei1 == ei2 && continue
        # The two colours must be 1 and 3 ("two connections of different colours").
        Set((arm_colour(nd, s), arm_colour(nd, s2))) == Set((1, 3)) || continue
        push!(mod1(s + 1, d) == s2 ? adj : far, (s, s2, ei1, ei2, t1, t2))
    end
    return vcat(adj, far)
end

"""
    _fr_commutation_merge(g::CircularGraph) -> Union{Nothing, CircularComboR}

C10. Two different non-`:braid` nodes of the 1/3 world are
connected by a CYCLICALLY ADJACENT pair of connecting edges, one of colour
1, one of colour 3, sitting on BOTH nodes at adjacent arms. With
`A = (i₁…iₙ a b)` and `B = (b a j₁…jₘ)` (the reversed order of the pair is
the opposite orientation), the result is the ONE node `(i₁…iₙ j₁…jₘ)`; all
`i` and `j` keep their existing wiring. No special case by arm count — if
two arms remain on each side, a regular 4-armed `1133` node results (see
the file header comment).

Any connections between the two nodes that do NOT belong to the pair stay
put and become self-connections (⇒ `needle` ⇒ 0 — the three-strand case).
`nothing` if no such pair exists.
"""
function _fr_commutation_merge(g::CircularGraph)
    for i in 1:length(g.nodes), j in 1:length(g.nodes)
        i < j || continue
        ndi, ndj = g.nodes[i], g.nodes[j]
        _circular_braidlike(ndi) && continue
        _circular_braidlike(ndj) && continue
        # pure 1/3 nodes, at least 4 arms (a crossing or larger)
        di, dj = arm_count(ndi), arm_count(ndj)
        (di >= 4 && dj >= 4) || continue
        all(c -> c == 1 || c == 3, arms(ndi)) || continue
        all(c -> c == 1 || c == 3, arms(ndj)) || continue

        for (si, si2, ei1, ei2, tj, tj2) in _circular_adjacent_conn_pairs(g, i, j)
            # OPPOSITE ORIENTATION, not just adjacency ("clockwise on
            # one, counter-clockwise on the other, i.e. opposite"): on `i`
            # the pair reads `(a, b)` at slots `(si, si2 = si+1)`; on `j` it
            # must be `(b, a)`, i.e. the far end of `si` is the SECOND of
            # the two j-slots. Formally: `tj = tj2 + 1` (cyclic). If it runs
            # same-orientation (`tj + 1 == tj2`), the two nodes aren't
            # oriented as in the spec and the rule does NOT apply — otherwise
            # a `1313` node would result instead of the identity.
            # OPPOSITE ORIENTATION, stated generally: walking on `i` from the
            # first to the second pair-slot forward, one must walk on `j`
            # from the corresponding second to the first — the pair order is
            # reversed. For ADJACENT pairs (`si2 == si+1`) that's exactly the
            # earlier condition `tj == tj2 + 1`; here it's phrased so it also
            # holds for separated slots.
            mod(si2 - si, di) == mod(tj - tj2, dj) || continue
            sj, sj2 = tj2, tj          # `(b a …)`: sj is the first j-slot

            # The surviving arms are the ARCS between the two pair slots —
            # not simply "from si2+1, di-2 of them". For an adjacent pair
            # that's the same arc (the other is empty); for separated slots
            # each node boundary splits into two arcs, and the planar splice
            # takes the arc of `i` AFTER its pair and the arc of `j` AFTER
            # its pair.
            arc(nd, from, to, dd) =
                [mod1(from + k, dd) for k in 1:(mod(to - from, dd) - 1)]
            slotsA = arc(ndi, si2, si, di)     # from si2 forward up to before si
            slotsB = arc(ndj, sj2, sj, dj)     # analogous on j
            length(slotsA) + length(slotsB) == (di - 2) + (dj - 2) || continue
            restA = [arm_colour(ndi, s) for s in slotsA]
            restB = [arm_colour(ndj, s) for s in slotsB]
            # the ports on the remaining arms, in the same order
            portA = [NodePort(i, s) for s in slotsA]
            portB = [NodePort(j, s) for s in slotsB]

            keep = Edge[]
            for (ei, e) in enumerate(g.edges)
                (ei == ei1 || ei == ei2) && continue      # the consumed pair
                push!(keep, e)
            end

            # UNIFORM, no case split by arm count:
            # `A = (i₁…iₙ a b)`, `B = (b a j₁…jₘ)` ⇒ `(i₁…iₙ j₁…jₘ)`. The `i`
            # and `j` keep their wiring; only `a`/`b` disappear. `restA` = the
            # iₖ (from the slot AFTER the pair, read cyclically), `restB` =
            # the jₖ analogously — concatenated, that's exactly the sequence
            # above.
            new_arms = vcat(restA, restB)
            # Fewer than 3 arms left is NOT a valid CircularNode. The case with 2
            # arms left on each side lands here: it is NOT resolved by C10,
            # it produces a 4-armed `1133` node (which stays, see the file
            # header comment). Fewer than 3 arms cannot occur anyway with two
            # nodes of arm count >= 4 (2+2 = 4).
            length(new_arms) >= 3 || continue
            newnd = circular_node(new_arms)
            # Slot assignment: new slots 1..|restA| to portA, then portB.
            remap = Dict{NodePort,Int}()
            for (k, p) in enumerate(portA); remap[p] = k; end
            for (k, p) in enumerate(portB); remap[p] = length(portA) + k; end

            # The new node takes i's place; j disappears.
            newnodes = CircularNode[]
            idxmap = Dict{Int,Int}()
            for (ni, nd) in enumerate(g.nodes)
                ni == j && continue
                push!(newnodes, ni == i ? newnd : nd)
                idxmap[ni] = length(newnodes)
            end
            fixp(p::Leaf) = p
            fixp(p::Circle) = p
            function fixp(p::NodePort)
                if p.node == i || p.node == j
                    haskey(remap, p) || return NodePort(idxmap[i], 1)   # shouldn't happen
                    return NodePort(idxmap[i], remap[p])
                end
                return NodePort(idxmap[p.node], p.slot)
            end
            newedges = [Edge(e.colour, fixp(e.a), fixp(e.b)) for e in keep]
            return CircularComboR(CircularGraph(g.word, newnodes, newedges))
        end
    end
    return nothing
end

"The far end of the edge at port `p` (excluding edges in `skip`): (port, colour, edge index)."
function _circular_far_end_of(g::CircularGraph, p::NodePort, skip::Set{Int})
    for (ei, e) in enumerate(g.edges)
        ei in skip && continue
        if e.a isa NodePort && e.a.node == p.node && e.a.slot == p.slot
            return (e.b, e.colour, ei)
        elseif e.b isa NodePort && e.b.node == p.node && e.b.slot == p.slot
            return (e.a, e.colour, ei)
        end
    end
    return nothing
end

# ---- MERGE LOCK for mixed 1/3 nodes -------------------------------------
#
# THE PROBLEM. A general Circular 1/3 node implicitly assumes every strand of a
# colour is fully connected. Single-edge contraction violates that: two
# `[1,3,1,3]` nodes joined by ONE edge give a 6-armed node — the exact same
# one a 4-valent 3-node with a passing-through 1-line would also give. Two
# different diagrams collapse to the same result — the merge is not
# information-preserving there.
#
# THE RULE: "If 2 circular 13 nodes are connected, they only get
# merged on ONE connection if one of the nodes is purely one colour. If both
# have both colours, no merge. With a double connection it's merged as
# before."
#
# THE SAME CONDITION, RESTATED:
#
#   "One connection: only merge if at least one node contains only one
#    colour.
#    Two connections: a) different colours ⇒ merge · b) same colours ⇒
#    needle = 0."
#
# STATEMENT 1 IS EXACTLY WHAT THE CODE BELOW DOES: the lock fires at
# `_circular_conn_count == 1` when BOTH nodes are `_circular_is_bicoloured` — so a
# merge only happens when at least one carries a single colour.
#
# STATEMENT 2 (two connections) belongs to the two rules responsible for it,
# BOTH of which come before C5: `_fr_mono_double` (same colour ⇒ 0, part b)
# and `_fr_two_adjacent_merge`/C12 resp. `_fr_commutation_merge`/C10
# (different colour ⇒ merge, part a). C5 itself checks nothing at ≥ 2
# connections.
#
# ⚠ NO KNOWN INSTANCE: the SAME-orientation `{1,3}` double connection — two
# mixed 1/3 nodes with exactly two DIFFERENT-coloured connections wired
# same-orientation instead of opposite — matches none of `_fr_mono_double`,
# C12 and C10 (C12/C10 want opposite orientation, `_fr_mono_double` equal
# colour) and would be cleaned up by C5 alone (a single-edge contraction
# leaving a self-connection, then `needle` ⇒ 0), where statement 2a calls for
# a merge. No such figure has been produced: the exhaustive two-node search
# finds no planar wiring, and over 700 double leaves every one of the 50
# `{1,3}` double connections is opposite-oriented. The planar case is the
# opposite one, and that is the case C12 merges: with `1 3` on one node and
# `3 1` on the other, the merged node carries the arms before the pair on the
# first node, then those after the pair on the second, then those before it,
# then those after the pair on the first. A blanket lock at
# `_circular_conn_count >= 2` is not the fix — it creates new fixed points
# where the registry would otherwise make progress (with the lock in place, no
# rule in the registry matches afterwards) — so the lock is OFF and the spot is
# order-dependent (see the ordering traps in
# the CircularDriver.jl header). Does this case need its own rule following the (a)/(b) scheme?
#
# THE CRITERION is "does the node carry both colours 1 AND 3". Locked are
# e.g. `[1,3,1,3]` against `[1,3,1,3]` or
# `[1,3,3,1]`; monochromatic and still mergeable are `[1,1,1]`, `[3,3,3]`,
# `[2,2,2]`, and every trivalent node. A node mixing 1/3 with colour 2 never
# occurs — `circular_node` (CircularGraph.jl:125-134) disallows it outright, so
# "carries both colours" and "pure 1/3 node with both colours" coincide.

"""
    _fr_merge(g::CircularGraph) -> Union{Nothing, CircularComboR}

Contracts the first edge between two DIFFERENT non-`:braid` nodes that are
NOT both degree 1 (that case is `_fr_barbell`, C3 — must be checked before
this rule). `nothing` if no such edge exists.

**THE MERGE CONDITION:**

> One connection: only merge if at least one node contains only one colour.
> Two connections: a) different colours ⇒ merge · b) same colours ⇒
> needle = 0.

**The single-connection branch of this rule is exactly the first sentence**
(merge lock): if the edge is the ONLY connection between its two nodes and
BOTH nodes carry both colour 1 and colour 3 ([`_circular_is_bicoloured`](@ref)),
it is **not** contracted; if at least one of them is monochromatic, it
merges.

**The two-connection case (statement 2) does not belong here**, it belongs
to [`_fr_mono_double`](@ref) (statement b, ⇒ 0) and
[`_fr_two_adjacent_merge`](@ref)/C12 resp. [`_fr_commutation_merge`](@ref)/C10
(statement a, ⇒ merge) — all three sit in the registry BEFORE C5. C5 checks
nothing at ≥ 2 connections, and must not. A same-orientation `{1,3}` double
connection would match none of the three and be cleaned up by C5 alone (⇒ 0
via `needle`) where a merge is called for — but no such figure is known to
exist; see the block comment above.
"""
function _fr_merge(g::CircularGraph)
    for (ei, e) in enumerate(g.edges)
        e.a isa NodePort && e.b isa NodePort || continue
        e.a.node == e.b.node && continue
        na, nb = g.nodes[e.a.node], g.nodes[e.b.node]
        # A braid-like node with `2k >= 8` arms ("gbraid") inherits the
        # `:braid` merge lock (`_circular_braidlike`) — otherwise
        # `merge_at_edge` THROWS here once C13 has produced one.
        _circular_braidlike(na) && continue
        _circular_braidlike(nb) && continue
        (arm_count(na) == 1 && arm_count(nb) == 1) && continue   # barbell, not here
        # Merge lock: ONE connection between two mixed 1/3 nodes. This is
        # statement 1 of the merge condition ("only merge if at least one
        # node contains only one colour"). Statement 2 (two
        # connections) is not here, see the block comment.
        if _circular_is_bicoloured(na) && _circular_is_bicoloured(nb) &&
           _circular_conn_count(g, e.a.node, e.b.node) == 1
            continue
        end
        return CircularComboR(merge_at_edge(g, ei))
    end
    return nothing
end
