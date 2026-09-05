# circular/rules/CircularWeight.jl — the termination weight for the circular-reduction driver.
#
# WHY. A pure node-count argument is not sufficient: it is demonstrably false
# for R4/R5 — there the node count RISES (one braid becomes two trivalents
# plus two dots, R5's "pinch" term). `circular_weight` provides a weight that
# provably falls strictly under every one of the 8 rules (C1–C8; C1–C5 in
# CircularRules.jl, C6–C8 in CircularBraidRules.jl) — see the delta table
# below.
#
# THE WEIGHT (lexicographic; Julia's tuple `<` compares exactly this way):
#
#   circular_weight(g) = (Σ arm_count over BRAID-LIKE nodes, node count,
#                    Σ arm_count over all nodes)
#
# Braid-like means `:braid`, including the 2k-armed gbraid. Component 1 is
# the ARM SUM over these nodes rather than their count, because a generalized
# braid node with `2k ≥ 8` arms does not count as an ordinary `:braid` — C14
# `dot_on_gen12` turns `gbraid(8) + dot(1)` into `braid(6) + trivalent(3)`,
# where node count and arm sum stay put but a plain node-count component would
# have risen. On any diagram without such generalized nodes, component 1 is
# exactly `6 · #:braid` — the same ordering as counting `:braid` nodes alone.
#
# Component 1 (braid arm sum) catches C6–C8: these three rules break down a
# `:braid` node (m=3) but may raise node count AND arm sum while doing so
# (R4/R5-style splits into several trivalents/dots). Components 2 (node count)
# and 3 (arm sum) catch C1–C5: these rules leave the braid weight unchanged
# (they don't touch a `:braid` node) and either lower the node count (nodes
# vanish: needle, barbell, bead, merge) or — in the one case where node count
# stays put (`_fr_dot_walks_on`'s main branch) — strictly lower the arm sum
# (one pair of arms is replaced by a single new dot arm).
#
# THIS IS THE TERMINATION PROOF for `reduce_circular` (CircularDriver.jl), not a debug
# aid: the driver `@assert`s after every rule step that each result term has
# strictly smaller `circular_weight` than the input, and that assert stays active
# in production.
#
# ---------------------------------------------------------------------------
# DELTA TABLE. `circular_weight` = (#braid, #nodes, Σarm_count).
#
# | Rule                         | before       | after (per term)          | strictly falls in |
# |-------------------------------|--------------|----------------------------|--------------------|
# | C1 `_fr_needle`               | (0, 1, 3)    | empty `CircularComboR`          | — (trivial: no terms) |
# | `_fr_mono_double`              | e.g. (0,2,10)| empty `CircularComboR`          | — (trivial, like C1) |
# | C2 `_fr_dot_walks_on` (main, `[1,3,3,1,3,3]`+dot) | (0, 2, 7) | (0, 2, 5) | arm sum (node count unchanged) |
# | C2 `_fr_dot_walks_on` (two-arm collapse, `[1,3,1,3]`+dot) | (0, 2, 5) | (0, 1, 1) | node count + arm sum |
# | C3 `_fr_barbell`               | (0, 2, 2)    | (0, 0, 0)                  | node count |
# | C4 `_fr_bead` (isolated)       | (0, 1, 2)    | (0, 0, 0)                  | node count |
# | C4 `_fr_bead` (in context)     | (0, 2, 5)    | (0, 1, 3)                  | node count + arm sum |
# | C5 `_fr_merge` (splice collapse, `[1,1,1]`+dot) | (0, 2, 4) | (0, 0, 0) | node count + arm sum |
# | C5 `_fr_merge` (general branch, `[1,1,1,3,3]`+dot) | (0, 2, 6) | (0, 1, 4) | node count + arm sum |
# | C6 `_fr_dot_into_braid`       | (1, 2, 7)    | term A (0,3,5), term B (0,1,1) | braid weight (node count/arm sum are allowed to rise or vary — term A does rise) |
# | C7 `_fr_braid_back`           | (2, 2, 12)   | term id (0,0,0), term pinch (0,4,8) | braid weight (pinch term: node count 2→4, arm sum 12→8, both allowed to move either way) |
# | C8 `_fr_braid_relation`       | (2, 3, 15)   | (1, 2, 9)                  | braid weight (one braid node remains; node count 3→2 also falls here) |
#
# C10/C12, the double-connection rules (circular/rules/CircularMergeRules.jl):
#
# | Rule                          | before (varies)   | after (varies)       | delta (CONSTANT)  |
# |--------------------------------|--------------------|----------------------|--------------------|
# | C12 `_fr_two_adjacent_merge`  | e.g. (0,2,10),(0,2,8) | e.g. (0,1,6),(0,1,4) | **(0, −1, −4)**, no exception |
# | C10 `_fr_commutation_merge`   | e.g. (0,10,26)     | e.g. (0,9,22)         | **(0, −1, −4)**, no exception |
#
# For both rules the delta is EXACTLY `(0, −1, −4)`, independent of the arm
# counts involved: component 1 stays 0 (both rules exclude braid-like nodes),
# component 2 (node count) falls by exactly 1 (two nodes become one),
# component 3 (arm sum) falls by exactly 4 (the four slots of the two consumed
# connecting edges).
#
# RESULT: the weight falls strictly and lexicographically for all 9 rules, on
# every result term.
#
# C12 matches only differently-coloured double connections; the same-coloured
# case is `_fr_mono_double` (⇒ empty combo, trivially terminating like C1).

# ---------------------------------------------------------------------------
# THE SECOND READING: (arm sum, separating edges, node count).
#
# The reading above is `circular_arm_weight`. It carries C14 `dot_on_gen12` and
# CANNOT carry C19 `trivalent_into_gen12`, the inverse direction. Since
# `C19 ∘ C19 ∘ C14 ∘ C14 = id` bit-exactly, no weight falls under both, so at
# most one of the two may be registered — and each needs its own reading.
#
# WHAT C19 DOES TO THE ARM SUM. C19 swallows a trivalent (3 arms) into a
# `2k`-node and leaves a DOT (1 arm) on the new `2(k+1)`-node:
#
#     (2k) + 3   ↦   (2k + 2) + 1        total arm sum UNCHANGED
#
# — while the BRAID arm sum rises `2k ↦ 2k+2`, which is why the first reading
# cannot carry it. Every OTHER rule strictly lowers the total arm sum. So the
# total arm sum orders all rules but C19, and C19 is decided by the second
# component:
#
#     separating edges = edges whose two darts lie in DIFFERENT regions
#                        (the edges of the dual graph).
#
# A DOT edge never separates — the face walk runs into the excursion and back
# out, so the same region lies on both sides. C19 turns a separating edge (the
# one to the trivalent) into a dot edge, hence the count falls by one. C14,
# running the other way, raises it by one at unchanged arm sum, which is
# precisely why it may not be registered once C19 is.
#
# Component 3 (node count) is a tie-breaker only.
#
# THE THIRD READING: nodes as their trivalent chain. The `m`-valent merge
# (C19 on a node with `m` arms, `2k + m ↦ 2(k+m−2)` arms plus `m−2` dots) is
# bit-identical to `m−2` successive trivalent merges, but under the plain arm
# sum the one-step version RISES by `2(m−3)` while the step-by-step version is
# constant: the plain arm sum is not invariant under writing one node as its
# trivalent chain (`m` against `3(m−2)`). `circular_chain_arms` removes exactly
# that difference, and then the `m`-valent merge falls in component 2, like the
# `m = 3` case always does.
#
# `CIRCULAR_WEIGHT_MODE` pins one reading without a source edit, so a
# measurement can put them side by side. ⚠ With `:arms` pinned while the merge
# direction is on, C19 trips the driver's assert — that mode is for measuring,
# not for running.
# ---------------------------------------------------------------------------

"""
    CIRCULAR_TRIVALENT_MERGE

**The direction switch.** `true` (the default) means trivalents are MERGED INTO
braid-like nodes: C19 `trivalent_into_gen12` is live, C14 `dot_on_gen12` is not,
and [`circular_weight`](@ref) is [`circular_sep_weight`](@ref). `false` reverses
all three — C14 live, C19 dead, [`circular_arm_weight`](@ref).

ONE switch for all three, because they are one decision: C14 and C19 are
inverse, so at most one may run, and each needs its own weight. Setting the
pieces inconsistently is what this switch exists to prevent — read it, do not
work around it.

The rules read it at RUNTIME (their matchers return `nothing` when their
direction is off), so a measurement can flip it in a live session.
"""
const CIRCULAR_TRIVALENT_MERGE = Ref(true)

"""
    CIRCULAR_MVALENT_MERGE

Whether C19 also fires on a one-coloured node with `m > 3` arms — the
`m`-valent merge, `2k + m ↦ 2(k+m−2)` arms plus `m−2` dots. Default **`true`**.

It is one switch for the rule and its weight, because the two belong together:
the `m`-valent step is bit-identical to `m−2` successive trivalent merges, but
under the plain arm sum the one-step version RISES while the step-by-step
version is constant. So with this on, `circular_weight` reads as
[`circular_chain_weight`](@ref), which weighs a node with `m ≥ 3` arms as the
chain of `m−2` trivalents it stands for and removes exactly that difference;
with it off, the plain [`circular_sep_weight`](@ref). With
[`CIRCULAR_TRIVALENT_MERGE`](@ref) off it has no effect at all.
"""
const CIRCULAR_MVALENT_MERGE = Ref(true)

"""
    CIRCULAR_WEIGHT_MODE

Which reading [`circular_weight`](@ref) uses. `:auto` (the default) follows the
direction switches — [`circular_chain_weight`](@ref) when
[`CIRCULAR_MVALENT_MERGE`](@ref) is on, else [`circular_sep_weight`](@ref) when
[`CIRCULAR_TRIVALENT_MERGE`](@ref) is on, else [`circular_arm_weight`](@ref).
`:arms`, `:sep` and `:chain` pin one reading regardless of the direction; that
combination is for MEASURING. All three are 3-tuples compared
lexicographically.
"""
const CIRCULAR_WEIGHT_MODE = Ref(:auto)

"""
    circular_weight(g::CircularGraph) -> Tuple{Int,Int,Int}
    circular_weight(d::CircularDecorated) -> Tuple{Int,Int,Int}

The termination weight for the circular-reduction driver (`reduce_circular`,
CircularDriver.jl), selected by [`CIRCULAR_WEIGHT_MODE`](@ref) over the three
readings [`circular_arm_weight`](@ref), [`circular_sep_weight`](@ref) and
[`circular_chain_weight`](@ref). For `CircularDecorated`, only the underlying
`CircularGraph` is weighed — the labels carry no information relevant to
termination.

**This is the termination proof, not a debug aid**: the driver checks via an
active `@assert` after every rule step that each result term has strictly
smaller `circular_weight` than the input.
"""
function circular_weight(g::CircularGraph)
    m = CIRCULAR_WEIGHT_MODE[]
    m === :arms  && return circular_arm_weight(g)
    m === :sep   && return circular_sep_weight(g)
    m === :chain && return circular_chain_weight(g)
    CIRCULAR_TRIVALENT_MERGE[] || return circular_arm_weight(g)
    return CIRCULAR_MVALENT_MERGE[] ? circular_chain_weight(g) : circular_sep_weight(g)
end
circular_weight(d::CircularDecorated) = circular_weight(d.graph)

"""
    circular_arm_weight(g) -> Tuple{Int,Int,Int}

`(Σ arm_count over braid-like nodes, node count, Σ arm_count)` — the reading
described in the first header block. It carries C14 `dot_on_gen12` and cannot
carry C19 `trivalent_into_gen12`.
"""
circular_arm_weight(g::CircularGraph) = (
    sum(nd -> _circular_braidlike(nd) ? arm_count(nd) : 0, g.nodes; init = 0),
    length(g.nodes),
    sum(arm_count, g.nodes; init = 0),
)
circular_arm_weight(d::CircularDecorated) = circular_arm_weight(d.graph)

"""
    circular_separating_edges(g::CircularGraph) -> Int

Edges whose two darts lie in DIFFERENT regions — the edges of the dual graph.
A DOT edge never separates (the face walk runs into the excursion and back out,
so the same region lies on both sides), which is why swallowing a trivalent into
a node and leaving a dot in its place lowers this count by one. Component 2 of
[`circular_sep_weight`](@ref) and the component that decides C19.
"""
function circular_separating_edges(g::CircularGraph)
    dr, t = _circular_dart_region(g)
    sep = 0
    for e in g.edges
        haskey(t.port_dart, e.a) || continue
        d  = t.port_dart[e.a]
        d2 = t.darts[d].rev
        R1 = get(dr, d, 0); R2 = get(dr, d2, 0)
        (R1 == 0 || R2 == 0 || R1 == R2) && continue
        sep += 1
    end
    return sep
end
circular_separating_edges(d::CircularDecorated) = circular_separating_edges(d.graph)

"""
    circular_sep_weight(g) -> Tuple{Int,Int,Int}

`(Σ arm_count over all nodes, separating edges, node count)`, lexicographic.
Every rule but C19 falls in component 1; C19 keeps it constant and falls in
component 2. See the second header block.
"""
circular_sep_weight(g::CircularGraph) =
    (sum(arm_count, g.nodes; init = 0), circular_separating_edges(g), length(g.nodes))
circular_sep_weight(d::CircularDecorated) = circular_sep_weight(d.graph)

"""
    circular_chain_arms(nd::CircularNode) -> Int

What a single node contributes to component 1 of [`circular_chain_weight`](@ref):
a NON-braid node with `m ≥ 3` arms counts `3(m−2)` — as much as the chain of
`m−2` trivalents it stands for — while braid-like nodes and dots or beads count
their arm count. For `m = 3` both readings agree (`3 = 3·1`).

⚠ It must be every non-braid node, not just the single-coloured ones: C2
`dot_walks_on` takes a MIXED `133133` node to a MONO `3333` node, and counting
only the mono side as a chain would make that step rise.
"""
circular_chain_arms(nd::CircularNode) =
    (arm_count(nd) >= 3 && !_circular_braidlike(nd)) ? 3 * (arm_count(nd) - 2) :
                                                       arm_count(nd)

"""
    circular_chain_weight(g) -> Tuple{Int,Int,Int}

`(Σ circular_chain_arms, separating edges, node count)` — [`circular_sep_weight`](@ref)
with non-braid nodes weighed as their trivalent chain, see
[`circular_chain_arms`](@ref). The reading `circular_weight` uses when
[`CIRCULAR_MVALENT_MERGE`](@ref) is on.
"""
circular_chain_weight(g::CircularGraph) =
    (sum(circular_chain_arms, g.nodes; init = 0), circular_separating_edges(g),
     length(g.nodes))
circular_chain_weight(d::CircularDecorated) = circular_chain_weight(d.graph)
