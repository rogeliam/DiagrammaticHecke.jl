# circular/rules/CircularRegionRules.jl — REGION-based rules (the second step of the
# region matcher).
#
# IDEA: a rule matches not on node shape but on the BOUNDARY WORD of a region
# (circular/CircularBoundaryWord.jl) — invariant under merging, so the case split of
# the C-series drops out at match time. The SURGERY still uses the existing
# bodies (`_fr_braid_back_at`, `_fr_braid_relation_at`) — no second copy; the
# region only supplies candidates GEOMETRICALLY instead of via an O(n²) node
# scan.
#
# MATCHING is up to ROTATION of the boundary word only, NEVER reflection —
# the same choice as `circular_canonical_key` (circular/CircularGraph.jl): equality holds up
# to rotation, `flip` is a separate operation. A mirrored variant would be its
# own rule.
#
# NO PATTERN DSL: each rule below uses a one-off letter predicate; the
# invariance comes from the region scan, not from a pattern language.
#
# ZAMO DIRECTION — lives in its own file, `circular/rules/CircularZamoRegion.jl`.
# Zamo is characterized by its **6** inner regions: all four
# Zamo sides have 6 inner regions of boundary-word length 3, whose corners
# cover all 7 nodes.
#
# Zamo is NOT registered in `CIRCULAR_REGION_RULES` below and gets no own entry:
# it is a RELATION, not a reduction — `circular_weight` does not fall (7 ↦ 7
# nodes), and the assert in CircularDriver.jl forbids that in any `CircularRule`
# registry. Its place is a dedicated step in `reduce_to_circular_leave`, next to
# `circular_zamo_step`. The construction there is the same as here: locate the
# STELLE geometrically (there over MULTIPLE regions and their shared edges,
# `_fzr_try`), delegate the SURGERY to the existing implementation
# (`apply_zamo_rule_combo`).
#
# PLANARITY PRECONDITION: region rules read the face tracer; on non-planar
# wiring (`euler != 2`) faces are meaningless and the rules simply don't match
# (conservative — the C fallbacks in `CIRCULAR_RULES_REGION` still apply).
# Existing C7/C8 test fixtures from test/circular.jl are wired with
# `euler = −4`: the ALGEBRAIC rules don't care, but the tracer does — the
# region tests (test/circularregionrules.jl) therefore use planar fixtures.
#
# TERMINATION: the three rules below call exactly the C1/C7/C8 surgeries — the
# delta table in circular/CircularWeight.jl applies verbatim.

# ---- the scan ------------------------------------------------------------

"""
    _frr_scan(g::CircularGraph, len::Int, try_match) -> Union{Nothing, CircularComboR}

The shared core of the region rules: iterates over all INNER regions of `g`
with boundary-word length `len` and all rotations of the boundary word, and
calls `try_match(letters)` (the rotated letter sequence). The first
non-`nothing` return wins (first-match, as everywhere in the driver).
"""
function _frr_scan(g::CircularGraph, len::Int, try_match::Function)
    for w in circular_boundary_words(g)
        length(w) == len || continue
        is_interior(w) || continue
        for r in 0:(len - 1)
            letters = [w.letters[mod1(j + r, len)] for j in 1:len]
            res = try_match(letters)
            res === nothing || return res
        end
    end
    return nothing
end

_frr_corner_node(l::CircularBoundaryLetter) =
    l.corner.vertex[1] === :node ? l.corner.vertex[2] : 0

# ---- C1 in the region registry --------------------------------------------
#
# C1 needs no region version. Its pattern is an EDGE pattern — `e.a.node ==
# e.b.node`, regardless of what the edge encloses — and that is the right
# formulation: by planarity two arms of a node can only be self-connected if
# they are adjacent and everything between them bottoms out on dots, so a
# self-loop never encloses anything substantial. A boundary-word scan
# restricted to length 1 would be a narrower reading of the same rule, not a
# genuinely different one. The region registry therefore takes `_fr_needle`
# itself, under its own name.

# ---- Pilot 2: braid_back_region (region version of C7) --------------------

"""
    _frr_braid_back(g::CircularGraph) -> Union{Nothing, CircularComboR}

Region version of C7 (`_fr_braid_back`): candidates `(a, b)` come from a
BIGON region (boundary-word length 2) between two distinct `:braid` nodes;
guards and surgery are unchanged, `_fr_braid_back_at`.
"""
function _frr_braid_back(g::CircularGraph)
    return _frr_scan(g, 2, function (ls)
        a, b = _frr_corner_node(ls[1]), _frr_corner_node(ls[2])
        (a != 0 && b != 0 && a != b) || return nothing
        (g.nodes[a].kind === :braid && arm_count(g.nodes[a]) == 6 &&
         g.nodes[b].kind === :braid && arm_count(g.nodes[b]) == 6) || return nothing
        return _fr_braid_back_at(g, min(a, b), max(a, b))
    end)
end

# ---- Pilot 3: braid_relation_region (region version of C8) ----------------

"""
    _frr_braid_relation(g::CircularGraph) -> Union{Nothing, CircularComboR}

Region version of C8 (`_fr_braid_relation`): candidates `(b1, b2, tv)` come
from a TRIANGLE region (boundary-word length 3) with corners at one pure
trivalent and two `:braid` nodes (the path-A triangle); guards and surgery
are unchanged, `_fr_braid_relation_at`. Replaces the O(braids²·trivs) scan
with a geometric candidate search.
"""
function _frr_braid_relation(g::CircularGraph)
    return _frr_scan(g, 3, function (ls)
        tv = _frr_corner_node(ls[1])
        b1 = _frr_corner_node(ls[2])
        b2 = _frr_corner_node(ls[3])
        (tv != 0 && b1 != 0 && b2 != 0) || return nothing
        (tv != b1 && tv != b2 && b1 != b2) || return nothing
        _circular_is_pure_trivalent(g.nodes[tv]) || return nothing
        (g.nodes[b1].kind === :braid && arm_count(g.nodes[b1]) == 6 &&
         g.nodes[b2].kind === :braid && arm_count(g.nodes[b2]) == 6) || return nothing
        return _fr_braid_relation_at(g, min(b1, b2), max(b1, b2), tv)
    end)
end

# ---- Registry ---------------------------------------------------------------

"""
The three region rules as `CircularRule`s (via `_lift`, label transfer as
everywhere).
"""
const CIRCULAR_REGION_RULES = CircularRule[
    CircularRule(:needle,                _lift(_fr_needle)),
    CircularRule(:braid_back_region,     _lift(_frr_braid_back)),
    CircularRule(:braid_relation_region, _lift(_frr_braid_relation)),
]

# `CIRCULAR_RULES_REGION` — the cross-check registry: `CIRCULAR_RULES` with
#   * `:needle` REPLACED by the region version (no interaction with the
#     ordering traps in CircularDriver.jl — C2/C3 before C5, C12 before C10, C13
#     last stay untouched);
#   * `:braid_back_region`/`:braid_relation_region` inserted BEFORE their
#     C7/C8 counterparts, which stay as fallback (the channel versions
#     `_fr_*_gen` cover the general-13 channel, which the bigon/triangle
#     region does not see — there the region is NOT a bigon, because the
#     channel node sits in between).
# `CIRCULAR_RULES` itself stays unchanged (the cross-check oracle).
function _build_circular_rules_region()
    out = CircularRule[]
    for r in CIRCULAR_RULES
        if r.name === :braid_back
            push!(out, CircularRule(:braid_back_region, _lift(_frr_braid_back)))
            push!(out, r)
        elseif r.name === :braid_relation
            push!(out, CircularRule(:braid_relation_region, _lift(_frr_braid_relation)))
            push!(out, r)
        else
            push!(out, r)
        end
    end
    return out
end

"Cross-check registry: `CIRCULAR_RULES` with the region versions inserted (see above)."
const CIRCULAR_RULES_REGION = _build_circular_rules_region()
