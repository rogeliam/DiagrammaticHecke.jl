# circular/rules/CircularBraidChannels.jl — channel normalization at braid pairs
# (derived and measured).
#
# THE PROBLEM. C7 (`_fr_braid_back`, R5) and C8 (`_fr_braid_relation`, R9) mirror
# the analogous rules for plain (non-circular) diagrams: C7 requires THREE DIRECT
# edges between the two braids, C8 two direct plus one PURE trivalent. Since
# circular nodes exist, a general-13 node can sit in such a channel instead —
# then neither rule fires, even though the situation is the same. A minimal
# example is known.
#
# THE FIX: given two connected braids, look at the
# neighbouring arms on the braids. If a general-13 node sits there: with only
# 2 arms, force the 13-node apart, which reduces to the classic case; with more
# than 2 arms, pull a trivalent out. So `13133` becomes `1313` connected with
# `333`.
#
# The node is normalized FIRST — the existing rules apply unchanged
# afterward.
#
# ⚠ WHY THE NORMALIZATION LIVES INSIDE THE RULE and not as its own step in the
# registry: on the split-open diagram, `_fr_merge` (C5) fires FIRST — it sits
# before C7/C8 in `CIRCULAR_RULES` — and immediately merges the fresh trivalent back
# into the 13-node. As a standalone step the split would therefore not just be
# useless but a CYCLE (split → merge → split → …). That's why both happen here
# in ONE step; no intermediate state where C5 could interfere ever exists.
#
# TERMINATION, measured: `circular_weight` decreases despite the split, because its
# first component is the number of braids and the rule consumes one braid — at
# the test bench `(2,3,17) ↦ (1,3,13)`. The assert in the driver
# (`CircularDriver.jl`) still holds.
#
# FORWARD REFERENCE: `expand_circular_merge` lives in `circular/rules/CircularInverseRules.jl`,
# included well after this file (the registry `CIRCULAR_RULES` needs the functions
# here already at its own definition). The call sits inside a function body and
# is therefore resolved only at run time, by which point the module is fully
# loaded. No need to move anything, but keep this in mind when reordering includes.

"""
    _circular_braid_channel(g::CircularGraph, a::Int, b::Int) -> nothing | NamedTuple

A node `v` that sits between the `:braid` nodes `a` and `b` in ONE channel:
exactly one arm into `a`, exactly one arm into `b`, both arms the same colour
`c`, and the two slots cyclically ADJACENT.

Adjacency isn't a convenience assumption, it's measured: over all planar
arrangements with `n = 2, 3, 4` arms of the channel colour, the
two connecting slots are without exception cyclically adjacent. The channel's
bigon region lies between them, and no further arm can reach outside it.

Returns `(v, k1, k2, colour, nc, d)`: `k1 → k2` is the clockwise arc between
the two connecting arms, `nc` the total number of `v`'s arms in the channel
colour (the classification criterion), `d` the arm count of `v`.
"""
function _circular_braid_channels(g::CircularGraph, a::Int, b::Int)
    out = NamedTuple[]
    for (v, nd) in enumerate(g.nodes)
        (v == a || v == b) && continue
        _circular_braidlike(nd) && continue
        d = arm_count(nd)
        sa = Int[]; sb = Int[]
        for (_, pv, po) in _circular_edges_at_node(g, v)
            po isa NodePort || continue
            po.node == a && push!(sa, pv.slot)
            po.node == b && push!(sb, pv.slot)
        end
        (length(sa) == 1 && length(sb) == 1) || continue
        p, q = sa[1], sb[1]
        c = arm_colour(nd, p)
        c == arm_colour(nd, q) || continue
        k1, k2 = mod1(p + 1, d) == q ? (p, q) :
                 mod1(q + 1, d) == p ? (q, p) : (0, 0)
        k1 == 0 && continue                      # not adjacent ⇒ not applicable
        nc = count(s -> arm_colour(nd, s) == c, 1:d)
        push!(out, (v = v, k1 = k1, k2 = k2, colour = c, nc = nc, d = d))
    end
    return out
end

"The first channel node, or `nothing`. Convenience wrapper over `_circular_braid_channels`."
_circular_braid_channel(g::CircularGraph, a::Int, b::Int) =
    (chs = _circular_braid_channels(g, a, b); isempty(chs) ? nothing : chs[1])

"""
    _circular_pull_trivalent(g::CircularGraph, a::Int, b::Int) -> nothing | CircularGraph

"More than 2 arms ⇒ pull a trivalent out." If a node with `nc ≥ 3` arms of the
channel colour sits between the braids `a` and `b`, a PURE trivalent `[c,c,c]`
carrying the two braid connections is split off from it; the rest stays behind
as its own node. `[1,3,1,3,3] ↦ [3,3,3] + [1,3,1,3]`.

This is literally `C5⁻¹` applied at exactly this spot, i.e. `expand_circular_merge`
with the arc `k1..k2` — the function itself checks that merging back gives the
input again. Nothing here is guessed.
"""
function _circular_pull_trivalent(g::CircularGraph, a::Int, b::Int)
    for ch in _circular_braid_channels(g, a, b)
        ch.nc >= 3 || continue
        # `conn_slot_a`/`conn_slot_b` are NOT a mathematical choice, but
        # `circular_canonical_key` sees them — and `expand_circular_merge` checks with
        # `check = true` that merging back reproduces `g`. Only a few
        # combinations close the round trip (see its docstring), so they're
        # tried in turn; the first one that passes the check is a correct
        # split. The trivalent has 3 arms, the rest `d - 1`.
        #
        # ⚠ SECOND, the COMPLEMENTARY arc is also offered. `expand_circular_merge`
        # builds its second arc as `[s for s in 1:d if s ∉ firstarc]`, i.e.
        # ASCENDING rather than cyclic; if the first arc runs across the slot
        # boundary (e.g. `{3,4}` at `d = 5` ⇒ cyclic rest `5,1,2`, but
        # ascending `1,2,5`), the round trip doesn't close and the function
        # throws. Offering the same cut the other way round — first arc = the
        # long one — makes the rest `{k1,k2}`, which is guaranteed to be
        # correct. Which of the two new nodes ends up the trivalent doesn't
        # matter to the rule.
        cuts = ((ch.k1, ch.k2),
                (mod1(ch.k2 + 1, ch.d), mod1(ch.k1 - 1, ch.d)))
        for (j1, j2) in cuts, csa in 1:ch.d, csb in 1:ch.d
            h = try
                expand_circular_merge(g, ch.v, j1, j2, ch.colour;
                                 conn_slot_a = csa, conn_slot_b = csb)
            catch
                nothing        # round trip doesn't close ⇒ this slot choice doesn't fit
            end
            h === nothing && continue
            return h
        end
    end
    return nothing
end

"""
    _circular_dissolve_through(g::CircularGraph, a::Int, b::Int) -> nothing | CircularGraph

"Only 2 arms ⇒ force the 13-node apart." If the channel node carries exactly
TWO arms of the channel colour, it's a pass-through line of that colour with a
foreign-coloured line beside it. Both colours commute (`13 ↔ 31`), so the node
decomposes into a direct channel edge and a free arc — leaving the classic C7
case.

With EXACTLY TWO foreign-coloured arms both parts would be **beads** (2-armed
nodes), which `circular_node` forbids and `_fr_bead` (C4) normalizes to edges anyway
— so the edges are built directly here instead of going through
`expand_circular_merge`.

**More than two foreign arms:** then, besides the channel edge, a genuine,
SINGLE-coloured node with `d − 2` arms remains. Its slot geometry is
**measured** (fixture `[1,1,1,3,3]` between two 23-braids): out of all
`3! = 6` assignments of the remaining arms to the new node's slots, boundary
word, `circular_degree`, `euler == 2` and
`check_wiring` survive for **exactly three** — and those are exactly the three
cyclic rotations of the order-preserving assignment (the three reversals fail
at `euler = 0` + `:reversed_node`). Since a `CircularNode`'s identity is its
CYCLIC arm sequence, the three are the same object: all three give bit-identical
C7 results (same `circular_canonical_key`). The arc built is therefore the one
**behind `k2`**, i.e. the remaining arms in their cyclic order.

Still `nothing` if the remaining arms aren't single-coloured (can't happen for
a two-coloured `CircularNode`, but stays as a guard) or a self-loop hangs off the
channel node.
"""
function _circular_dissolve_through(g::CircularGraph, a::Int, b::Int)
    for ch in _circular_braid_channels(g, a, b)
        ch.nc == 2 || continue
        nd = g.nodes[ch.v]
        # The remaining arms in CYCLIC order, starting behind `k2`.
        # (Taking ascending `1:d` is the same only when the arc `k1→k2` doesn't
        # run across the slot boundary — the same pitfall as in
        # `_circular_pull_trivalent`.)
        rest = [mod1(ch.k2 + i, ch.d) for i in 1:(ch.d - 2)]
        length(rest) >= 2 || continue
        other_colour = arm_colour(nd, rest[1])
        all(s -> arm_colour(nd, s) == other_colour, rest) || continue

        opposite = Dict{Int,Port}()
        drop = Set{Int}()
        for (ei, pv, po) in _circular_edges_at_node(g, ch.v)
            opposite[pv.slot] = po
            push!(drop, ei)
        end
        length(opposite) == ch.d || continue        # self-loop or similar

        if length(rest) == 2
            keep = Edge[e for (ei, e) in enumerate(g.edges) if !(ei in drop)]
            push!(keep, Edge(ch.colour, opposite[ch.k1], opposite[ch.k2]))
            push!(keep, Edge(other_colour, opposite[rest[1]], opposite[rest[2]]))
            return _circular_delete_nodes(g, Set([ch.v]), keep)
        end

        gplus = CircularGraph(g.word, vcat(g.nodes, [circular_node(fill(other_colour, length(rest)))]),
                         g.edges)
        new_node = length(gplus.nodes)
        keep = Edge[e for (ei, e) in enumerate(gplus.edges) if !(ei in drop)]
        push!(keep, Edge(ch.colour, opposite[ch.k1], opposite[ch.k2]))
        for (i, s) in enumerate(rest)
            push!(keep, Edge(other_colour, opposite[s], NodePort(new_node, i)))
        end
        return _circular_delete_nodes(gplus, Set([ch.v]), keep)
    end
    return nothing
end

# The braid pairs C7/C8 can apply to at all.
function _circular_braid_pairs(g::CircularGraph)
    bs = [i for (i, nd) in enumerate(g.nodes) if nd.kind === :braid && arm_count(nd) == 6]
    [(a, b) for a in bs for b in bs if a < b &&
        Set(arms(g.nodes[a])[1:2]) == Set(arms(g.nodes[b])[1:2])]
end

# Resolve all pass-through channels of a braid pair (at most two; node numbers
# shift as we go, so re-search each time).
function _circular_dissolve_all_through(g::CircularGraph)
    out = g
    for _ in 1:2
        progress = false
        for (a, b) in _circular_braid_pairs(out)
            h = _circular_dissolve_through(out, a, b)
            h === nothing && continue
            out = h; progress = true; break
        end
        progress || break
    end
    return out === g ? nothing : out
end

"""
    _fr_braid_back_gen(g::CircularGraph) -> Union{Nothing, CircularComboR}

C7 with channel normalization: try the unmodified rule first, otherwise
resolve the pass-through nodes between two braids and try again. The classic
case runs bit-identically to `_fr_braid_back` as before.
"""
function _fr_braid_back_gen(g::CircularGraph)
    r = _fr_braid_back(g)
    r === nothing || return r
    gg = _circular_dissolve_all_through(g)
    gg === nothing && return nothing
    return _fr_braid_back(gg)
end

"""
    _fr_braid_relation_gen(g::CircularGraph) -> Union{Nothing, CircularComboR}

C8 with channel normalization: try the unmodified rule first; otherwise, at
each braid pair, pull a trivalent out of the general-13 node (and, if the
other channel carries a pass-through, resolve it too) and try again. The
classic case runs bit-identically to before.
"""
function _fr_braid_relation_gen(g::CircularGraph)
    r = _fr_braid_relation(g)
    r === nothing || return r
    for (a, b) in _circular_braid_pairs(g)
        gg = _circular_pull_trivalent(g, a, b)
        gg === nothing && continue
        r2 = _fr_braid_relation(gg)
        r2 === nothing || return r2
        hh = _circular_dissolve_all_through(gg)          # second channel, if it's a pass-through
        hh === nothing && continue
        r3 = _fr_braid_relation(hh)
        r3 === nothing || return r3
    end
    return nothing
end
