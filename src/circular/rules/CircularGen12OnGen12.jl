# circular/rules/CircularGen12OnGen12.jl — C16 `gen12_on_gen12`: TWO gbraid nodes
# joined by THREE edges collapse into ONE gbraid plus two trivalent nodes.
#
# THE EQUATION
#
#     gbraid (nv arms) + gbraid (nb arms), three shared edges
#         ↦  one gbraid (nv + nb − 8 arms) + two trivalent nodes
#
# ONE term, coefficient 1. Both sides have the same boundary and
# `circular_degree = −4`.
#
# THE SAME COMPUTATION AS C15, just with different arm counts. The patched
# boundary has `(8−3) + (8−3) = 10` connections instead of eight; its colour
# sequence again contains exactly two adjacent same-coloured pairs, each
# getting a trivalent node of that colour, and the remaining `10 − 2 = 8`
# strands alternate — hence one gbraid and no `:braid`. C16 therefore
# shares its surgery with C15: `_circular_glue3_rewire` in `CircularBraidOnGen12.jl`
# handles both cases.
#
# WHERE THE TRIVALENT NODES GO is, as with C15, forced and readable off the
# boundary word: `1212212122` (splice at an odd leaf) has its two doublings
# in colour 2, `1121211212` (even leaf) has them in colour 1 — and that
# colour is exactly what the two trivalent nodes carry. Nothing to choose.
#
# At every candidate site the pattern is a `CIRCULAR_RULES` FIXED POINT (weight
# `(16,2,16)`), while the 16 combinations of preimages (`expand_gbraid` per
# node) fall into FOUR classes — the same situation as before C15, just one
# level up:
#
#   | Class | Combinations | Result                                     | Weight      |
#   |---|---|---|---|
#   | 1 | 8 | gbraid + gbraid (= the starting pattern)         | `(16,2,16)` |
#   | 2 | 6 | gbraid + two trivalent nodes                     | `(8,3,14)`  |
#   | 3 | 1 | `:braid` + gbraid + two trivalent nodes          | `(14,4,20)` |
#   | 4 | 1 | ditto, different key                              | `(14,4,20)` |
#
# All four have degree −4, the same boundary word, and all four are fixed
# points. Class 2 is the lightest — and it is exactly the right-hand side of
# this rule.
#
# TERMINATION. The structural identity: the single new node carries
# `nv + nb − 8` arms and the two trivalents are not braid-like. So the BRAID
# ARM SUM — component 1 of `circular_arm_weight` — goes `nv + nb ↦ nv + nb − 8`,
# a drop of exactly 8 at every pair of arm counts, not just at a measured one.
# Under the `sep` and `chain` readings that `CIRCULAR_WEIGHT_MODE = :auto`
# selects today, component 1 falls by 2 instead. Either way `circular_weight`
# falls and the assert in `CircularDriver.jl` holds. Measured at the arm-count
# pairs `(8,6)`, `(10,6)`, `(12,6)`, `(8,8)`, `(10,8)`, `(10,10)` and `(12,8)`:
# the rule fires, the result is planar and wiring-clean, `circular_degree` and
# the boundary word are unchanged, and the weight falls under all three
# readings.
#
# At `(8, 8)` the arm reading is `(16,2,16) ↦ (8,3,14)`.
#
# ARM COUNTS: any two gbraids (`arm_count >= 8` each). The surgery is not tied
# to a count — `_circular_glue3_rewire_at` derives everything from
# `p = nv + nb − 6` and rejects a pair that does not fit. What the rule DOES
# fix is the number of shared edges: EXACTLY THREE. One, two and four or more
# are the business of C18 and C22.
#
# TERMINOLOGY: "gbraid" means a braid-like node with `2k >= 8` arms, i.e.
# `kind === :braid && arm_count >= 8`; there is no node kind of its own for it.

"""
    _fr_gen12_on_gen12(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C16 `gen12_on_gen12`.** Two gbraid nodes (`arm_count >= 8` each) joined by
**three** shared edges are replaced by **one** gbraid and **two trivalent
nodes** —
at the two same-coloured neighbour pairs of the patched boundary. One term,
coefficient 1.

The surgery is the shared core [`_circular_glue3_rewire`](@ref) (see
`CircularBraidOnGen12.jl`); only the pattern lives here. Since both nodes have the
same kind, the pair is tried in **both orders** — which of the two carries
the descending slot sequence is decided by the embedding, not by the
numbering.

Derivation and measurements: file header.
"""
function _fr_gen12_on_gen12(g::CircularGraph)
    gb = [i for (i, nd) in enumerate(g.nodes)
              if nd.kind === :braid && arm_count(nd) >= 8]
    for a in 1:length(gb), b in 1:length(gb)
        a == b && continue
        h = _circular_glue3_rewire(g, gb[a], gb[b])
        h === nothing && continue
        return CircularComboR(h)
    end
    return nothing
end
