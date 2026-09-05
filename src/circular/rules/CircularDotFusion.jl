# circular/rules/CircularDotFusion.jl — the GENERAL D4 RULE:
# the dot analogue of the general parallel rule (CircularRexFusion.jl).
#
# Analogous to the general parallel rule: if a dot sits in a region where the
# word can end on the dot's letter, a new relation must be derived here too.
#
# THE TRIGGER is the SAME criterion as for the parallel rule: the shortest
# unreduced transition (`circular_unreduced_transition`) — here in its `:dot`
# kind: an `s`-dot inside region `R` whose reduced region word `w` has `s` as a
# right descent (the word can end on the dot's letter).
#
# THE ONE NEW IDEA (everything else is CircularRexFusion.jl verbatim): the cut
# curve is EXTENDED THROUGH THE DOT'S CAP EDGE. The tree path to `R` crosses
# edges of colours `w`; appending the dot's edge (colour `s`) gives a curve of
# colours `w·s` — exactly the word the derived relation
# `circular_rex_relation(w, s)` lives on. Splicing its `rhs` replaces the local
# identity on `w·s` by the lower terms; the dot stays on the far side of the
# last cut and caps the `s`-structure of each term. No fusion-partner gate is
# needed: the dot IS the partner.
#
# HOW THE CASES LINE UP with the existing special rules:
#   :direct    zero braid moves — the rex cycle is trivial and the relation is
#              NOT derivable (measured, CircularRexFusion.jl file head). This is
#              the classic D4 / D4-node case; the step declines.
#   :dihedral  one braid move — the derived version of
#              what `circular_regionword_dihedral_step` (gated dot slide) does
#              by hand surgery.
#   :long      several braid moves — the derivation runs here directly.
#
# THIS IS LABEL-FREE (same policy as CircularRexFusion/CircularDotSlide):
# non-trivial labels REFUSE, they are not guessed through the splice.
#
# TERMINATION is not claimed here: no driver hook, switch default OFF, the
# step is callable independently. Candidate measure: the length of the
# shortest unreduced region word.
#
# THE RESULT IS THE RAW SPLICE: in every term one dot ends
# up sitting directly on a node — that is not a drawing bug but the unreduced
# splice. The rule that cleans it up exists and fires: `:merge` (C5,
# CircularMergeRules.jl), exactly ONE step per term, and every term loses
# exactly one node — nothing grows.
# The `reduce = true` keyword below runs exactly this clean-up
# (`reduce_circular_combo`); default `false` keeps the raw splice inspectable.

"""
    CIRCULAR_DOT_FUSION_ENABLED

Switch for the driver hook of [`circular_dot_fusion_step`](@ref) in
`reduce_to_circular_leave`, default **`true`**. The step is callable
independently of it.
"""
const CIRCULAR_DOT_FUSION_ENABLED = Ref(true)

"""
    circular_dot_fusion_step(fdm::CircularDecoratedMorphism; reduce = false,
                             kwargs...)
        -> Union{Nothing, NamedTuple}

The general D4 rule (file head): if the shortest unreduced transition is a
`:dot` transition — an `s`-dot in a region whose reduced word `w` has `s` as a
right descent — splice the derived rex-cycle relation on `w·s` into the tree
path to that region EXTENDED by the dot's cap edge. Returns

* `combo`      — the sum that replaces the diagram (label-free terms),
* `transition` — the transition of [`circular_unreduced_transition`](@ref),
* `curve`      — the spliced edges (colours `w·s`; last one = the dot's edge),
* `dot`        — the dot node (`transition.to`),
* `sides`      — the assignment anchored on the identity,
* `pair`       — the i-pair that carried the fusion inside the relation,
* `simple_curve` — the step-3a classification (reported only, `_crf_simple`).

`nothing` whenever a gate declines: non-trivial labels, no transition, an
`:edge` transition (the rex fusion's case), `case == :direct` (classic D4's
case, and the rex cycle is trivial there), a dot edge already on the tree
path, an identity splice that does not reproduce the graph, or a term that
does not splice. Raises only if the relation itself is not derivable.

`reduce = true` cleans every term of `combo` with
[`reduce_circular_combo`](@ref) before returning — one `:merge` (C5) step per
term takes the dot off its node (file head). Default `false`: the RAW splice
stays inspectable, and nothing existing changes.

Remaining `kwargs` are passed to [`circular_rex_relation`](@ref) (`dir`,
`recompute`, `maxpairs`).
"""
function circular_dot_fusion_step(fdm::CircularDecoratedMorphism;
                                  reduce::Bool = false, kwargs...)
    all(isone, fdm.region_labels) && isone(fdm.outer_label) || return nothing
    m = fdm.m
    g = m.graph

    T = circular_unreduced_transition(m)
    T === nothing && return nothing
    T.kind === :dot || return nothing           # `:edge` is CircularRexFusion's case
    T.case === :direct && return nothing        # classic D4; relation underivable

    p = circular_path_edges(m, T.region)
    p === nothing && return nothing
    T.edge in p.edges && return nothing         # cap edge must extend, not repeat
    curve = vcat(p.edges, T.edge)
    cw = vcat(p.word, T.colour)
    [g.edges[e].colour for e in curve] == cw || return nothing

    rel = circular_rex_relation(p.word, T.colour; kwargs...)

    # anchor the assignment on the identity, which must give `g` back
    anchor = circular_splice(g, curve, circular_splice_identity(cw))
    anchor === nothing && return nothing
    circular_canonical_key(anchor.graph) == circular_canonical_key(g) || return nothing

    n = length(cw)
    out = CircularComboR()
    for (d, coeff) in pairs_of(rel.rhs)
        (all(isone, d.region_labels) && isone(d.outer_label)) || return nothing  # label-free
        tm = CircularMorphismGraph(d.graph, 0, n)             # cuts of the derivation
        (bottom(tm) == cw && top(tm) == cw) || return nothing
        r = circular_splice(g, curve, tm; sides = anchor.sides)
        r === nothing && return nothing
        out = out + coeff * CircularComboR(circular_decorated(r.graph))
    end
    reduce && (out = reduce_circular_combo(out))
    return (combo = out, transition = T, curve = curve, dot = T.to,
            sides = anchor.sides, pair = rel.pair,
            simple_curve = _crf_simple(g, curve))
end
