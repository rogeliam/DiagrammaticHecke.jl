# circular/rules/CircularGen12Null.jl — C22 `gen12_on_gen12_null`: two
# braid-like nodes sharing FOUR OR MORE edges are the ZERO element.
#
# THE EQUATION
#
#     :braid (na arms) + :braid (nb arms), k >= 4 shared edges   ↦   0
#
# Empty term list, no right-hand side.
#
# WHY. A fourfold connection contains a threefold one, so the three-edge rule
# applies, and the surplus edge closes on itself afterwards:
#
#   1. Three of the `k` shared edges sit cyclically adjacent with descending
#      partner slots — the pattern of C15/C16. Their surgery
#      (`_circular_glue3_rewire_at`, CircularBraidOnGen12.jl) applies: ONE term,
#      coefficient 1, both braid-like nodes replaced by ONE new node plus two
#      trivalents.
#   2. The SURPLUS shared edges join two ports that both land on that ONE new
#      node — so each becomes a SELF-LOOP, i.e. a needle.
#   3. C1 `_fr_needle` makes a needle 0. Hence the original figure is 0.
#
# NOTHING IS ASSUMED HERE. The matcher does not merely count edges: it RUNS the
# C15/C16 surgery on the three-edge window and checks that the result is wired
# (`check_wiring`) and really carries a self-loop. Only then is the term 0, so
# every hit of this branch is a proof by construction rather than an
# extrapolation. Every valid three-edge window of the figures below behaves the
# same way — surgery, then `reduce_circular` reaching `[:needle]`:
#
#   | figure    | k | windows | after the surgery          |
#   |-----------|--:|--------:|----------------------------|
#   | `8 + 8`   | 4 |       2 | `braid/8`  + 2 trivalents  |
#   | `8 + 6`   | 4 |       2 | `braid/6`  + 2 trivalents  |
#   | `10 + 8`  | 4 |       2 | `braid/10` + 2 trivalents  |
#   | `8 + 8`   | 5 |       3 | `braid/8`  + 2 trivalents  |
#   | `10 + 8`  | 5 |       3 | `braid/10` + 2 trivalents  |
#   | `12 + 8`  | 5 |       3 | `braid/12` + 2 trivalents  |
#   | `10 + 10` | 6 |       4 | `braid/12` + 2 trivalents  |
#
# THE `6 + 6` CASE. The surgery guards itself with `p = na + nb − 6 >= 8` (a
# 4-armed braid does not exist), so at `6 + 6` it declines, and the three-edge
# rule there is C7 `braid_back` (id + pinch) instead. The same argument goes
# through: apply C7 to three of the four edges, with the fourth cut open and
# glued back, and
#
#   term `id`    → the fourth edge closes into a FREE CIRCLE   → C1 needle → 0
#   term `pinch` → the fourth edge joins the two new trivalents: a DIGON
#                  → `mono_double` (or merge + needle)         → 0
#
# Hence the SECOND branch of the matcher: two braid-like nodes with a FOUR-EDGE
# WINDOW — four shared edges, cyclically consecutive on both nodes, `b`-slots
# ascending and `v`-slots descending (the C15/C16 shape extended by one edge) —
# are 0 for every arm count. The window is what makes the argument local: three
# consecutive edges carry the three-edge rule, the adjacent fourth becomes the
# needle resp. the digon.
#
# ⚠ SCOPE: a surplus edge that is NOT adjacent to a three-edge window (something
# enclosed between two shared edges) is covered only by the FIRST branch, i.e.
# where the surgery runs (`na + nb >= 14`) and the self-loop is actually seen.
# At `6 + 6` such a figure would need a pinch-term argument of its own — the
# fourth edge lands on the DOT legs and gives a barbell, not a digon — and is
# not claimed here.
#
# TERMINATION: trivial — the rule returns the EMPTY term list, so there is
# nothing left to rewrite. `circular_weight` is not consulted (`_lift` handles
# an empty `CircularComboR` without touching the weight assert).

"""
    _circular_shared_slot_pairs(g::CircularGraph, v::Int, b::Int) -> Vector{Tuple{Int,Int}}

The slot pairs `(slot at v, slot at b)` of the edges connecting the nodes `v`
and `b`, in edge order.
"""
function _circular_shared_slot_pairs(g::CircularGraph, v::Int, b::Int)
    pairs = Tuple{Int,Int}[]
    for e in g.edges
        (e.a isa NodePort && e.b isa NodePort) || continue
        if e.a.node == v && e.b.node == b
            push!(pairs, (e.a.slot, e.b.slot))
        elseif e.b.node == v && e.a.node == b
            push!(pairs, (e.b.slot, e.a.slot))
        end
    end
    return pairs
end

"""
    _circular_glue3_windows(g::CircularGraph, v::Int, b::Int)
        -> Vector{Vector{Tuple{Int,Int}}}

All C15/C16-shaped **three-edge windows** among the shared edges of the
braid-like nodes `v` and `b`: three edges whose `b`-slots run cyclically
ascending and whose `v`-slots run cyclically descending — the criterion
[`_circular_gen12_braid_glue`](@ref) applies, but without its "exactly three
shared edges" restriction.

Empty vector if there is no such window (or fewer than three shared edges).
"""
function _circular_glue3_windows(g::CircularGraph, v::Int, b::Int)
    out = Vector{Tuple{Int,Int}}[]
    nv, nb = arm_count(g.nodes[v]), arm_count(g.nodes[b])
    pairs = _circular_shared_slot_pairs(g, v, b)
    length(pairs) >= 3 || return out
    for start in pairs
        b1, j1 = start[2], start[1]
        window = Tuple{Int,Int}[]
        for i in 1:3
            want = (mod1(j1 - (i - 1), nv), mod1(b1 + i - 1, nb))
            want in pairs || break
            push!(window, want)
        end
        length(window) == 3 && allunique(window) && push!(out, window)
    end
    return out
end

"""
    _circular_has_glue4_window(g::CircularGraph, v::Int, b::Int) -> Bool

Is there a **four-edge window** among the shared edges of `v` and `b`: four
edges whose `b`-slots run cyclically ascending and whose `v`-slots run
cyclically descending — the shape of [`_circular_glue3_windows`](@ref) with one
more edge? Such a figure is 0 for every arm count: three of its edges carry the
three-edge rule, and the fourth closes into a needle resp. a digon.
"""
function _circular_has_glue4_window(g::CircularGraph, v::Int, b::Int)
    nv, nb = arm_count(g.nodes[v]), arm_count(g.nodes[b])
    pairs = _circular_shared_slot_pairs(g, v, b)
    length(pairs) >= 4 || return false
    for start in pairs
        b1, j1 = start[2], start[1]
        window = Tuple{Int,Int}[(mod1(j1 - (i - 1), nv), mod1(b1 + i - 1, nb)) for i in 1:4]
        (all(w in pairs for w in window) && allunique(window)) && return true
    end
    return false
end

"Does `h` carry a needle — an edge with both ends on the SAME node, or a free circle?"
function _circular_has_self_loop(h::CircularGraph)
    for e in h.edges
        (e.a isa Circle || e.b isa Circle) && return true
        (e.a isa NodePort && e.b isa NodePort && e.a.node == e.b.node) && return true
    end
    return false
end

"""
    _fr_gen12_on_gen12_null(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C22 `gen12_on_gen12_null`.** Two braid-like nodes sharing **four or more**
edges are the zero element: the C15/C16 surgery on three of them leaves the
surplus edge as a **needle**. Returns the EMPTY `CircularComboR` (= 0) on a
hit, `nothing` otherwise.

Two branches (see the file header): (1) the C15/C16 surgery on a three-edge
window is actually carried out and the result checked for a self-loop
(`na + nb >= 14`); (2) a FOUR-EDGE WINDOW — four consecutive shared edges in the
C15/C16 shape — is 0 for every arm count, in particular at `6 + 6`.
"""
function _fr_gen12_on_gen12_null(g::CircularGraph)
    braids = [i for (i, nd) in enumerate(g.nodes)
                  if nd.kind === :braid && arm_count(nd) >= 6]
    for ia in eachindex(braids), ib in eachindex(braids)
        ia == ib && continue
        v, w = braids[ia], braids[ib]
        length(_circular_shared_slot_pairs(g, v, w)) >= 4 || continue
        for window in _circular_glue3_windows(g, v, w)
            h = _circular_glue3_rewire_at(g, v, w, window)
            h === nothing && continue
            isempty(check_wiring(h)) || continue
            _circular_has_self_loop(h) || continue
            return CircularComboR()          # the empty sum = 0
        end
        _circular_has_glue4_window(g, v, w) && return CircularComboR()
    end
    return nothing
end
