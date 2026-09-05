# circular/rules/CircularGen12TwoEdges.jl — C18 `gen12_two_edges`: an 8-armed
# generalized braid node (a `:braid` node with eight arms, "gbraid") and an
# ordinary `:braid`, connected by edges or by a channel node.
#
# APPROACH. The `:gbraid` is expanded via `expand_gbraid` (C13⁻¹) onto one of
# its **four frozen preimages**, and the ORDINARY rules fire on the result. No
# search is needed for the preimages themselves — they are literals
# (`GBRAID_PREIMAGES`) — and an expansion is kept only when the subsequent
# reduction is genuinely lighter.
#
# THE TRIGGERING CASE. Of the
# 16 preimage combinations of the pattern "two `:gbraid`, three edges" (C16), 14
# land on the C16 right-hand side; two stay stuck at
# `:braid` + `:gbraid` + two trivalents — weight `(14,4,20)`, a fixed point.
# Counting their connections: exactly **two** edges sit between `:braid` and
# `:gbraid`. This is the pattern C18 catches.
#
# Over both anchor slots `k = 1, 2` and all four outliers:
#
#   * 32 planar expansions per outlier (4 variants × 8 rotations),
#   * of which **16 lower the weight**, and these 16 all give ONE single
#     `circular_canonical_key` — the right-hand side is therefore well defined and
#     does not depend on which preimage is chosen;
#   * the rules that fire are `[:merge, :braid_relation, :gen12_merge]` resp.
#     `[:merge, :braid_relation, :merge, :gen12_merge]` — C5, C8, C13, all
#     existing rules, no new surgery;
#   * the result is **identical to the C16 right-hand side** for the same
#     `k` (`:gbraid` + two trivalents, weight `(8,3,14)`), by
#     `circular_canonical_key`.
#
# Together with C16, this rule covers all 16 preimage combinations, not just 14.
#
# TERMINATION. The rule only returns results whose `circular_weight` is STRICTLY
# smaller than the input's — checked in the body (`w1 < w0` per term). Example:
# `(14,4,20) ↦ (8,3,14)`. The expansion step itself
# raises the weight (`expand_gbraid` is an expanding rule), but it is only an
# INTERMEDIATE state — only the reduced result is returned.
#
# NO RECURSION. The intermediate reduction runs over `CIRCULAR_RULES` WITHOUT C18
# itself (filtered by name), so it cannot call itself. The reference to
# `CIRCULAR_RULES` sits inside a function body and is resolved only at call time —
# the same forward reference as `expand_circular_merge` in `CircularBraidChannels.jl`.
#
# ⚠️ SCOPE: an 8-armed `:gbraid`, partner is either a 6-armed `:braid` or a
# second 8-armed `:gbraid`, with **one, two, or four** shared edges, same
# colour pair. EXACTLY THREE edges are excluded — that is the
# pattern of C15 (`:braid`, contiguous slots) resp. C16 (`:gbraid`); letting
# C18 match there too would duplicate their right-hand sides. Three SCATTERED
# edges (which C15's connectivity check rejects) remain open.
#
# The synthetic minimal
# patterns — a bare `:gbraid` + `:braid` with one or four edges, `:gbraid` +
# `:gbraid` with two edges, and also `:gbraid` + `:braid` with two edges but
# WITHOUT the real outlier's trivalents — are all `CIRCULAR_RULES` fixed points,
# and NONE of the 32 planar expansions per case lowers the weight. Widening
# the scope is therefore a no-op on these bare patterns (the weight guard
# rejects everything, they stay normal forms) and only helps in embedded
# situations where an expansion genuinely uncovers something lighter — as with
# the two-edge outlier, which needs its two trivalents to fall.
#
# ⚠️ TWO ORDINARY BRAIDS DOUBLY CONNECTED are a different case and stay as
# they are: a surgery with its own right-hand side would change the DL
# inventory. C18 does NOT do that — it invents no right-hand side, only lets
# the existing rules compute on a preimage, and the result is the C16 normal
# form that already exists.
#
# TERMINOLOGY: "gbraid" means a braid-like node with `2k >= 8` arms, i.e.
# `kind === :braid && arm_count >= 8`; there is no node kind of its own for
# it. The matchers here check `kind === :braid && arm_count == 8`.
#
# ---- CHANNEL SITES AND TARGETED PREIMAGE CHOICE ----------------------------
#
# At the four earlier outliers there is
# structurally an F2-3 channel
# (`nc = 3`, `t = 1`) between `:gbraid` and `:braid`, and
# `_fr_braid_relation_gen` fires at none of them — because `_circular_braid_pairs`
# (`CircularBraidChannels.jl`) requires `kind === :braid && arm_count == 6` at BOTH
# channel ends.
#
#   * The four outliers are already caught by the EXISTING C18 (each has two
#     direct `:gbraid`–`:braid` edges): under the full `CIRCULAR_RULES`, each
#     reduces in ONE step via `[:gen12_two_edges]`, `(14,4,20) ↦ (8,3,14)`.
#   * Per case, of the 32 planar expansions exactly **8 are F2-3-capable**
#     (`_fr_braid_relation` matches directly).
#   * All accepted computations run via `[:merge, :braid_relation, :merge,
#     :gen12_merge]` and give, per colour case, ONE `circular_canonical_key`
#     (`:gbraid` + two trivalents, `(8,3,14)`) — the right-hand side is well
#     defined and preimage-independent.
#
# TWO CHANGES follow from this:
#  (1) **Site search widened:** a `:gbraid` and a braid-like partner also
#      count as a site when connected via a CHANNEL NODE
#      (`_circular_gen12_channel_pairs`), not only via 1/2/4 direct edges. This
#      closes the matcher gap structurally: `_circular_braid_channels` already
#      recognizes channels generically, only the upstream site search was too
#      narrow. ⚠️ This is NOT channel normalization applied to `:gbraid` —
#      the channel is read only as a SITE; the computation still runs via
#      `expand_gbraid` + existing rules.
#  (2) **Targeted preimage choice:** expansions on which F2-3 fires directly
#      are tried FIRST (`_circular_gen12_expansions`, F2-3-capable ones up front).
#      The "first weight-lowering" path is the fallback behind them.
#
# ⚠️ SCOPE: `_circular_braid_channels` requires cyclically ADJACENT connection
# slots. Channels with non-adjacent slots are not seen by the widened site
# search either.
#
# ---- FACTORED OUT ----------------------------------------------------------
#
# The mechanism "try site × expansion, run `reduce_circular` WITHOUT itself, accept
# only if strictly lighter and label-free" now lives as the reusable combinator
# `circular_via_preimage` (circular/rules/CircularViaPreimage.jl) — future rules following the
# same scheme (C6/C7/C8 candidates) should use it too. `_fr_gen12_two_edges`
# below supplies only the sites (`pairs`) and the expansions
# (`_circular_gen12_expansions`) and calls the combinator.

"""
    _circular_gen12_partner_pairs(g::CircularGraph) -> Vector{Tuple{Int,Int,Int}}

The triples `(v, b, n)` of an 8-armed `:gbraid` `v` and a braid-like partner
`b` (6-armed `:braid` or a second 8-armed `:gbraid`) with the same colour
pair, connected by `n ∈ {1, 2, 4}` direct edges. `n == 3` is
excluded (the C15/C16 pattern, see file header). `:gbraid`–`:gbraid` pairs
appear in both orders — the body expands `v` in each case, so both sides get
tried.
"""
function _circular_gen12_partner_pairs(g::CircularGraph)
    out = Tuple{Int,Int,Int}[]
    for (v, nv) in enumerate(g.nodes)
        (nv.kind === :braid && arm_count(nv) == 8) || continue
        for (b, nb) in enumerate(g.nodes)
            b == v && continue
            ((nb.kind === :braid && arm_count(nb) == 6) ||
             (nb.kind === :braid && arm_count(nb) == 8)) || continue
            Set(arms(nv)[1:2]) == Set(arms(nb)[1:2]) || continue
            n = count(e -> e.a isa NodePort && e.b isa NodePort &&
                           Set((e.a.node, e.b.node)) == Set((v, b)), g.edges)
            n in (1, 2, 4) && push!(out, (v, b, n))
        end
    end
    return out
end

"""
    _circular_gen12_braid_two_edges(g::CircularGraph) -> Vector{Tuple{Int,Int}}

The original C18 matcher: `:braid` partners with **exactly two** edges only.
Today a filter over `_circular_gen12_partner_pairs`; still used by the tests as
the fixed point of the old semantics.
"""
_circular_gen12_braid_two_edges(g::CircularGraph) =
    Tuple{Int,Int}[(v, b) for (v, b, n) in _circular_gen12_partner_pairs(g)
                   if n == 2 && g.nodes[b].kind === :braid && arm_count(g.nodes[b]) == 6]

"""
    _circular_gen12_channel_pairs(g::CircularGraph) -> Vector{Tuple{Int,Int}}

The pairs `(v, b)` of an 8-armed `:gbraid` `v` and a braid-like partner `b`
(6-armed `:braid` or 8-armed `:gbraid`) with the same colour pair, connected
via at least one CHANNEL NODE — the sites where F2-3 structurally applies but
`_circular_braid_pairs` does not look.

Channel recognition is `_circular_braid_channels` (`CircularBraidChannels.jl`)
unchanged: it does not restrict `a`/`b` to `:braid`, only the channel node
itself must not be braid-like. The narrowness was entirely in
`_circular_braid_pairs`, which stays untouched (C7/C8 must not change).
"""
function _circular_gen12_channel_pairs(g::CircularGraph)
    out = Tuple{Int,Int}[]
    for (v, nv) in enumerate(g.nodes)
        (nv.kind === :braid && arm_count(nv) == 8) || continue
        for (b, nb) in enumerate(g.nodes)
            b == v && continue
            ((nb.kind === :braid && arm_count(nb) == 6) ||
             (nb.kind === :braid && arm_count(nb) == 8)) || continue
            Set(arms(nv)[1:2]) == Set(arms(nb)[1:2]) || continue
            isempty(_circular_braid_channels(g, v, b)) && continue
            push!(out, (v, b))
        end
    end
    return out
end

"""
    _circular_gen12_expansions(g::CircularGraph, v::Int) -> Vector{CircularGraph}

The planar expansions of the `:gbraid` `v` (4 preimages × 8 rotations),
**F2-3-capable ones first**: expansions on which `_fr_braid_relation` resp.
`_fr_braid_relation_gen` matches DIRECTLY come first.

This is the targeted preimage choice: *if the right preimage is chosen for a
gbraid/braid pair, F2-3 can be applied directly; for others it cannot — so the
right preimage has to be found locally.* 8 of the 32 planar
expansions per case are F2-3-capable (file header). The rest stay at the back,
so cases C18 already solved via another path keep running unchanged.
"""
function _circular_gen12_expansions(g::CircularGraph, v::Int)
    front = CircularGraph[]; back = CircularGraph[]
    for variant in 1:4, rot in 0:7
        e = expand_gbraid(g, v, variant; rot = rot)
        e === nothing && continue
        (is_wired(e) && euler(e) == 2 && isempty(check_wiring(e))) || continue
        f23 = _fr_braid_relation(e) !== nothing || _fr_braid_relation_gen(e) !== nothing
        push!(f23 ? front : back, e)
    end
    return vcat(front, back)
end

"""
    _fr_gen12_two_edges(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C18 `gen12_two_edges`.** If an 8-armed `:gbraid` is connected to a
braid-like partner (`:braid` with 6 or `:gbraid` with 8 arms; three edges
belong to C15/C16, see file header) by **one, two, or four** direct edges —
or via a channel node (`_circular_gen12_channel_pairs`) — the `:gbraid` is expanded
onto one of its four preimages (`expand_gbraid`, C13⁻¹) and the result is
reduced with the remaining rules. Expansions are tried in the order of
`_circular_gen12_expansions` (**F2-3-capable first**); the first one whose
reduction is **strictly lighter** than the input is taken.

The right-hand side does not depend on the choice of
preimage: across all four cases, the 16 weight-lowering expansions
give the same `circular_canonical_key`, namely the C16 normal form (see the file
header).

`nothing` if the pattern is absent or no expansion yields anything lighter.
"""
function _fr_gen12_two_edges(g::CircularGraph)
    # Sites: direct edges (the original C18 case) AND channel connections
    # (see file header). Order: direct sites first, then channel sites — so
    # the pre-existing cases stay bit-identical.
    pairs = unique(vcat([v for (v, _, _) in _circular_gen12_partner_pairs(g)],
                        [v for (v, _) in _circular_gen12_channel_pairs(g)]))
    # The mechanism (try site × expansion, run `reduce_circular` WITHOUT itself,
    # accept only if strictly lighter and label-free) lives in
    # `circular_via_preimage` (CircularViaPreimage.jl) — here we only supply the sites
    # and, per site, the expansions (F2-3-capable first).
    return circular_via_preimage(g, pairs, v -> _circular_gen12_expansions(g, v);
                             except_rule = :gen12_two_edges)
end
