# circular/rules/CircularDriver.jl — the Circular-reduction driver.
#
# Runs FULLY decorated: rules fire on `CircularDecorated`, the driver works on
# `CircularComboR` (= `CircularCombo{SoergelPoly}` over `CircularDecorated` representatives,
# see CircularDecorated.jl). There is no undecorated fast path and no fusion in the
# registry (`_fr_fusion`/`_circular_fuse_edge` are directed tools taking an edge
# argument, not an `apply(fd) -> Union{Nothing,CircularComboR}` contract).
#
# REGISTRY ORDER IS SEMANTIC — do not "clean up" or alphabetize. Documented
# ordering traps:
#   - `_fr_dot_walks_on` (C2) MUST come before `_fr_merge` (C5).
#   - `_fr_barbell` (C3) MUST come before `_fr_merge` (C5) (two degree-1 nodes
#     next to each other would otherwise make `_fr_merge` throw instead of
#     producing α_i).
#   - `_fr_mono_double` MUST come before `_fr_two_adjacent_merge` (C12): both
#     match the same site (two non-braid-like nodes with exactly two
#     adjacent, opposite-direction connections), complementary by color —
#     `_fr_mono_double` takes the SAME-color case (⇒ 0), C12 the
#     DIFFERENT-color case (⇒ merge). See the CircularRules.jl header comment
#     before `_fr_mono_double`.
#   - `_fr_two_adjacent_merge` (C12) MUST come before `_fr_commutation_merge`
#     (C10): C10 only covers the `{1,3}` color case. See the CircularRules.jl
#     header comment before C12 for details.
#
# The remaining order (C1 first = cheapest test; C6–C8 last = the three braid
# rules) is likewise taken from CircularRules.jl — see there.

# ---- A1: Registry -------------------------------------------------------

"""
    CircularRule

A local circular-graph rule: `apply(fd::CircularDecorated) -> Union{Nothing, CircularComboR}`
returns the linear combination `fd` rewrites to (0, 1, or several
`CircularDecorated` terms), or `nothing` if the rule does not match. 1:1 analogue
of `GraphRule`.
"""
struct CircularRule
    name::Symbol
    apply::Function
end

# ---- A2: lifting the eight CircularGraph rules to CircularDecorated ---------------
#
# Eight of the nine rules (CircularRules.jl) are written on `CircularGraph` (they don't
# know about labels); `_fr_barbell` (C3) already takes `CircularDecorated` and
# produces the α_i coefficient itself. `_lift` is the single wrapper for the
# other eight — not eight separate copies.

"""
    _lift(rule_fn) -> Function

Lifts a `CircularGraph` rule `rule_fn(g::CircularGraph) -> Union{Nothing,CircularComboR}`
(a `CircularComboR` over `CircularGraph` terms) to a `CircularDecorated` rule
`fd::CircularDecorated -> Union{Nothing,CircularComboR}` (over `CircularDecorated` terms):
calls `rule_fn` on `fd.graph`, then carries the existing labels onto each result
term via `_circular_transfer_labels`.

Errors from `_circular_transfer_labels` propagate unchanged — no `try`/`catch`, no
guessing, no silent fallback. The transfer only
throws when a merge would combine two NON-trivial labels; "non-trivial wins"
is the rule otherwise. There is therefore exactly ONE wrapper, with no special
case for a "bare" (all-labels-1) input — that path runs through unchanged too
(it is cheap and never throws there).
"""
function _lift(rule_fn::Function)
    return function (fd::CircularDecorated)
        res = rule_fn(fd.graph)
        res === nothing && return nothing
        out = CircularComboR()
        for (g2, coeff) in pairs_of(res)
            labels2 = _circular_transfer_labels(fd.graph, g2, fd.region_labels)
            # The outer label is carried along unchanged (the lifted rules
            # are purely structural, and the outer region stays the outer
            # region through any rewrite — it can neither merge nor split).
            _add!(out, CircularDecorated(g2, labels2, fd.outer_label), coeff)
        end
        return out
    end
end

const CIRCULAR_RULES = CircularRule[
    CircularRule(:needle,              _lift(_fr_needle)),               # C1
    CircularRule(:pitchfork,           _lift(_fr_pitchfork)),            # C21 — right behind C1: the same statement (term = 0) for the general figure
    CircularRule(:dot_walks_on,        _lift(_fr_dot_walks_on)),          # C2 — before merge
    CircularRule(:barbell,             _fr_barbell),                      # C3 — before merge, already takes CircularDecorated
    CircularRule(:bead,                _lift(_fr_bead)),                  # C4
    # ⚠ The house convention: the normal form is the NODE, not the identity —
    # two parallel edges of colour 1 and 3 that don't touch, connected with
    # matching colours, merge (P1, `circular_parallel_merge_step`), and {1,3}
    # nodes are allowed to stay together. Dissolving such a node back into two
    # strands is a legitimate relation, and it maps the normal forms one to one,
    # so the leaf count is the same either way; the direction here is the
    # convention, not a mathematical necessity.
    #
    # `_fr_mono_double`: same-colored double connection ↦ 0 — MUST come
    # before C12 (see the CircularRules.jl header comment before it);
    # complementary to C12 at the same site.
    CircularRule(:mono_double,         _lift(_fr_mono_double)),          # before C12
    CircularRule(:two_adjacent_merge,  _lift(_fr_two_adjacent_merge)), # C12 — before C10
    CircularRule(:commutation_merge,   _lift(_fr_commutation_merge)),  # C10 — MUST come before merge
    CircularRule(:merge,               _lift(_fr_merge)),                 # C5
    CircularRule(:dot_into_braid,      _lift(_fr_dot_into_braid)),        # C6
    # C7/C8 run through the CHANNEL formulation
    # (circular/rules/CircularBraidChannels.jl): first the unchanged rule, otherwise
    # normalize the general-13 node in the channel and try once more. The
    # classical case is bit-identical — `_fr_braid_back`/`_fr_braid_relation`
    # are called first and remain exported and separately tested.
    CircularRule(:braid_back,          _lift(_fr_braid_back_gen)),        # C7
    CircularRule(:braid_relation,      _lift(_fr_braid_relation_gen)),    # C8
    # C20 `:two_adjacent_dots` (circular/rules/CircularDotMerge.jl):
    # two dots on ADJACENT arms of a `2k`-node give the naked `2(k-1)`-node —
    # C14 + C5 in one step, verified bit-identical for 2k = 8/10 and across the
    # seam. Registered because it lowers arm sum AND node count (so the
    # termination measure of CircularWeight.jl is untouched), and because the
    # dot policy (`circular_dot_may_pass`) is NOT confluent without it: it
    # would leave „8-armed + two adjacent dots" as a second normal form of the
    # naked 6-armed node. Sits before C15.
    CircularRule(:two_adjacent_dots,   _lift(_fr_two_adjacent_dots)),      # C20
    # C14 `:dot_on_gen12` is gated by the dot policy
    # (`circular_dots_reducible`): it only fires where one of two cases holds
    # — two adjacent dots, or a colour left with fewer than three arms. It is
    # the ONLY rule that performs case (b) at `2k >= 10`, where the generic
    # preimage stage does not reach. C20 sits before it and does case (a) in one
    # step. Switch the gate off with `CIRCULAR_DOT_POLICY[] = false` to get the
    # ungated behaviour.
    #
    # ⚠ The entry is the COMPOSITE `_fr_dot_on_gen12_collapse`, not the single
    # step: one C14 step RISES by a separating edge, and only the fully
    # collapsed node is lighter. The composite does step and collapse in one go
    # and accepts only a strictly lighter result, so the driver's termination
    # assert holds. It is INDEPENDENT of `CIRCULAR_TRIVALENT_MERGE`, which
    # governs C19 alone: the C14↔C19 circle cannot start because C19 refuses to
    # build a node the policy would immediately collapse (the loop guard in
    # `_circular_mono_into_gen12_at`). See CircularDotOnGen12.jl.
    CircularRule(:dot_on_gen12,        _lift(_fr_dot_on_gen12_collapse)),  # C14
    CircularRule(:braid_on_gen12,      _lift(_fr_braid_on_gen12)),        # C15
    # C16 directly after C15: the same equation one level up (two gbraid
    # nodes instead of gbraid + `:braid`), the same surgery. It also
    # lowers the arm sum 16 ↦ 8, so it carries the same weight argument.
    # circular/rules/CircularGen12OnGen12.jl.
    CircularRule(:gen12_on_gen12,      _lift(_fr_gen12_on_gen12)),        # C16
    # C22 directly after C16: the SAME site with one shared edge too many. C16
    # needs exactly three, so four or more would be a fixed point for it; the
    # surgery on three of them leaves the surplus edge as a needle, hence 0.
    # Standing AFTER C16 keeps the three-edge case under its own name — C22 only
    # ever sees what C16 declined. circular/rules/CircularGen12Null.jl.
    CircularRule(:gen12_on_gen12_null, _lift(_fr_gen12_on_gen12_null)),   # C22
    # C18 after C16 and before C13: the TWO-EDGE gbraid–`:braid` case
    # (circular/rules/CircularGen12TwoEdges.jl). It has no surgery of its own — it
    # opens the gbraid via `expand_gbraid` (C13⁻¹) into a preimage and
    # lets the remaining rules compute, so it MUST come after any rule that
    # could handle the case itself (C15/C16) and before C13 (which would
    # otherwise just rebuild the gbraid). It only accepts its result when
    # `circular_weight` strictly drops (checked in the body), and calls the
    # registry WITHOUT itself — no recursion.
    CircularRule(:gen12_two_edges,     _lift(_fr_gen12_two_edges)),       # C18
    # ---- C13 COMES LAST -----------------------------------------------------
    # ⚠ The 2-braid + 2-trivalent → gbraid merge must not
    # have high priority; other rules go first: at two of four
    # preimages, giving C13 high priority makes it fire again immediately and
    # rebuild the gbraid, where the reduction then stalls, while the
    # other two preimages reach `:braid` + 2 trivalent nodes. C13 last
    # means exhausting everything else first, building the gbraid only
    # when nothing else applies. circular/rules/CircularGen12Merge.jl.
    CircularRule(:gen12_merge,         _lift(_fr_gen12_merge)),           # C13
    # ---- C19 DIRECTLY AFTER C13 ---------------------------------------------
    # C19 `:trivalent_into_gen12` (circular/rules/CircularDotMerge.jl) is the
    # INVERSE of C14 `dot_on_gen12`: a trivalent hanging on a braid-like
    # `2k`-node is swallowed, the arm becomes three arms, and the middle one
    # gets a dot. It RAISES the braid arm sum, so it cannot run under
    # `circular_arm_weight`; under `circular_sep_weight` it falls, because
    # swallowing the trivalent turns a separating edge into a dot edge, and a
    # dot edge has the same region on both sides. Which reading is in force is
    # decided by `CIRCULAR_TRIVALENT_MERGE`, the one switch that also decides
    # whether C14 or C19 runs — see CircularWeight.jl.
    #
    # It sits next to C13 for the same reason C13 comes last: other rules go
    # first when a large braid node is being BUILT.
    CircularRule(:trivalent_into_gen12,
                 _lift(g -> CIRCULAR_TRIVALENT_MERGE[] ?
                            _fr_trivalent_into_gen12(g) : nothing)),      # C19
    # C23 directly after C19: the same merge, on a trivalent the rule has to
    # free from a mixed {1,3} node first. It sits after C19 so that a trivalent
    # already lying there is taken first, and a node is only broken open when
    # nothing cheaper applies. It is the one rule in
    # `CIRCULAR_WEIGHT_ASSERT_EXEMPT` — see there and the box before
    # `CIRCULAR_MIXED_MERGE` (circular/rules/CircularDotMerge.jl).
    CircularRule(:mixed_into_gen12,    _lift(_fr_mixed_into_gen12)),      # C23
]

# ---- THE GENERIC PREIMAGE STAGE -------------------------------------------
#
# ONE driver step covers the widenings of C6/C7/C8: if nothing else matches,
# open up a circular node and try again. The fixed
# rule set does not recognize every case where two diagrams are equal — that
# is why the new-kind nodes exist. Working new rules for those nodes means
# looking at various preimages, so sometimes `reduce` has to be locally
# undone, to then apply other rules, with the goal that things get simpler
# afterward.
#
# THE STAGE IS NOT IN `CIRCULAR_RULES`, for three reasons:
#   1. **No existing behavior changes.** Every existing rule keeps priority;
#      the stage is only consulted once the loop over `rules` returns
#      `nothing` — i.e. only at diagrams that are fixed points for
#      `CIRCULAR_RULES`.
#   2. **No ordering trap.** It would otherwise have to be placed among
#      C14/C15/C16/C18 (all four touch the same gbraid); as a fallback the
#      question doesn't arise.
#   3. **No recursion.** Since it isn't in `rules`, the inner `reduce_circular`
#      can't call it again. Belt and braces: `circular_via_preimage` always passes
#      a fresh rule array (so the `rules === CIRCULAR_RULES` identity check below
#      fails there), and `_CIRCULAR_PREIMAGE_DEPTH` counts on top of that. Depth
#      stays at 1.
#
# SITES: only the 8-armed gbraid — the only node kind with a frozen
# preimage (`GBRAID_PREIMAGES` is fixed to 8, CircularGen12Expand.jl). The site
# list is capped by `maxsites`, because each site costs up to 32 expansions,
# each with a full `reduce_circular`.
#
# EXPANSIONS: reuses `_circular_gen12_expansions` (CircularGen12TwoEdges.jl) unchanged —
# four preimages × eight rotations, filtered by `is_wired`/`euler ==
# 2`/`check_wiring`, F2-3-capable ones first. The same pre-sorting as C18: it
# costs nothing and typically finds a usable expansion immediately.
#
# WHY THIS PRODUCES THE SAME RIGHT-HAND SIDES AS THE SPECIAL-CASE RULES:
#   * **C6 is contained in the preimage path.** An 8-armed gbraid with a
#     dot on each of its eight boundary leaves, preimage path without
#     `:dot_on_gen12`: 256 planar expansions, ALL weight-decreasing, one
#     `circular_canonical_key` per leaf, and 8/8 match the C14 result. Rule
#     sequences `[:dot_into_braid, :merge, …, :needle]` resp.
#     `[:merge, :braid_relation]`.
#   * **C15 is reproduced exactly.** 16 cases, 512 expansions, 384
#     weight-decreasing, one key per case, 16/16 match the C15 right-hand
#     side.
# The three acceptance conditions (weight strictly drops · no foreign labels
# · result independent of the chosen preimage) are unchanged from
# `circular_via_preimage`.

"""
    CIRCULAR_GENERIC_PREIMAGE_MAXSITES

How many expansion sites the generic preimage stage tries at most (default
4). Caps the runtime: each site costs up to 32 expansions, each with a full
`reduce_circular`. Kept as a `Ref` so measurements can tune it without editing the
source.
"""
const CIRCULAR_GENERIC_PREIMAGE_MAXSITES = Ref(4)

# Recursion guard, second lock (see box above): while > 0, an intermediate
# reduction of the stage is already running and the fallback is skipped.
# Not `task_local_storage` — the circular driver runs single-threaded.
const _CIRCULAR_PREIMAGE_DEPTH = Ref(0)

"""
    _fr_generic_preimage(g::CircularGraph; maxsites = CIRCULAR_GENERIC_PREIMAGE_MAXSITES[])
        -> Union{Nothing, CircularComboR}

The **generic preimage stage** at `CircularGraph` level: opens an 8-armed
gbraid via `expand_gbraid` (C13⁻¹) into one of its four frozen preimages
and lets the full `CIRCULAR_RULES` compute on it. The result is accepted only if
it is non-empty, EVERY term is strictly lighter (`circular_weight`), and no term
carries foreign labels — the three acceptance conditions of
`circular_via_preimage`.

`nothing` if there is no 8-armed gbraid, or no expansion reveals anything
lighter. That is the normal case: the stage is only consulted at diagrams
that are already fixed points for `CIRCULAR_RULES`.

**Not in `CIRCULAR_RULES`** — it is the fallback in `_circular_first_rule_match`, see
the box above.
"""
function _fr_generic_preimage(g::CircularGraph;
                               maxsites::Int = CIRCULAR_GENERIC_PREIMAGE_MAXSITES[])
    sites = Int[v for (v, nd) in enumerate(g.nodes)
                  if nd.kind === :braid && arm_count(nd) == 8]
    # THE DOT POLICY: this stage, like C14, would push a dot through an
    # 8-armed node. A site carrying a
    # dot it may not pass is therefore skipped — a DOTLESS 8-armed node stays
    # expandable (that is C13 work, not dot pushing).
    if CIRCULAR_DOT_POLICY[]
        with_dot = Set(p.node for e in g.edges for p in (e.a, e.b)
                     if p isa NodePort && ((e.a isa NodePort && arm_count(g.nodes[e.a.node]) == 1) ||
                                           (e.b isa NodePort && arm_count(g.nodes[e.b.node]) == 1)))
        sites = Int[v for v in sites if !(v in with_dot) || circular_dots_reducible(g, v)]
    end
    isempty(sites) && return nothing
    length(sites) > maxsites && (sites = sites[1:maxsites])
    return circular_via_preimage(g, sites, v -> _circular_gen12_expansions(g, v);
                             except_rule = nothing)
end

"""
    circular_generic_preimage(fd::CircularDecorated; maxsites = CIRCULAR_GENERIC_PREIMAGE_MAXSITES[])
        -> Union{Nothing, CircularComboR}
    circular_generic_preimage(g::CircularGraph; maxsites = …)

The generic preimage stage as a `CircularDecorated` step (labels are transferred
via `_lift`, like any rule). The driver calls exactly this in the fallback of
`_circular_first_rule_match`; it is public so measurements and tests can query the
stage on its own, without going through the driver.
"""
circular_generic_preimage(fd::CircularDecorated;
                      maxsites::Int = CIRCULAR_GENERIC_PREIMAGE_MAXSITES[]) =
    _lift(g -> _fr_generic_preimage(g; maxsites = maxsites))(fd)
circular_generic_preimage(g::CircularGraph;
                      maxsites::Int = CIRCULAR_GENERIC_PREIMAGE_MAXSITES[]) =
    circular_generic_preimage(circular_decorated(g); maxsites = maxsites)

# ---- A3 hookup: the termination assert -----------------------------------
#
# Checks that every term of a rule application's result has strictly smaller
# `circular_weight` than the input graph (circular/CircularWeight.jl; termination proof for
# all 9 rules, see the delta table there). An empty `CircularComboR` (C1
# `_fr_needle` ⇒ 0) has no terms — the loop is then empty and the assert is
# trivially satisfied, which is correct.
#
# THIS ASSERT STAYS ACTIVE: it is the termination proof, not a debug
# artifact — do not remove it as "noise".
"""
    CIRCULAR_WEIGHT_ASSERT_EXEMPT

The rules whose termination does NOT come from `circular_weight`, and which
`_assert_circular_weight_drops` therefore skips. A rule belongs here only with
a termination argument of its own written down at its definition.

C23 `mixed_into_gen12` is the one entry: it raises the arm sum by six
(`2k + m ↦ 2k + m + 6`), so no reading of the weight catches it, and it
terminates instead because every piece of the broken-open node has strictly
fewer arms than the node it came from.
"""
const CIRCULAR_WEIGHT_ASSERT_EXEMPT = Set{Symbol}([:mixed_into_gen12])

function _assert_circular_weight_drops(rule_name::Symbol, fd::CircularDecorated, result::CircularComboR)
    rule_name in CIRCULAR_WEIGHT_ASSERT_EXEMPT && return nothing
    w0 = circular_weight(fd)
    for (t, _) in pairs_of(result)
        w1 = circular_weight(t)
        @assert w1 < w0 (
            "reduce_circular: rule :$rule_name does not strictly decrease circular_weight " *
            "($w0 -> $w1) — termination proof violated")
    end
    return nothing
end

# find the first matching rule at ALL (scans CIRCULAR_RULES in order), returns
# (name, result-CircularComboR) or `nothing`.
#
# FALLBACK: if no rule matches, the generic preimage stage
# (`_fr_generic_preimage`, box above) is consulted, and its result reported
# under the name `:generic_preimage` — it then shows up in `reduce_circular`'s
# `history` like any other rule name. Two conditions gate it:
#
#   * `rules === CIRCULAR_RULES` — the stage always computes internally with the
#     full registry; pulling it into a restricted run (`without-C18` & co. in
#     measurements and tests) would be a silent change of meaning. Callers
#     that restrict the registry keep getting exactly the restricted
#     computation.
#   * `_CIRCULAR_PREIMAGE_DEPTH[] == 0` — recursion guard. The first condition
#     already handles this (the intermediate reduction gets a fresh array);
#     the counter is the backstop that still holds if someone passes
#     `CIRCULAR_RULES` through directly.
function _circular_first_rule_match(fd::CircularDecorated, rules)
    for r in rules
        res = r.apply(fd)
        if res !== nothing
            _assert_circular_weight_drops(r.name, fd, res)
            return (r.name, res)
        end
    end
    if rules === CIRCULAR_RULES && _CIRCULAR_PREIMAGE_DEPTH[] == 0
        _CIRCULAR_PREIMAGE_DEPTH[] += 1
        res = try
            circular_generic_preimage(fd)
        finally
            _CIRCULAR_PREIMAGE_DEPTH[] -= 1
        end
        if res !== nothing
            _assert_circular_weight_drops(:generic_preimage, fd, res)
            return (:generic_preimage, res)
        end
    end
    return nothing
end

# ---- A4: driver functions (1:1 analogues) -------------

"""
    reduce_circular(fd::CircularDecorated; rules = CIRCULAR_RULES, maxsteps = 1000)
        -> (CircularComboR, history::Vector{Symbol})
    reduce_circular(g::CircularGraph; rules = CIRCULAR_RULES, maxsteps = 1000)
        -> (CircularComboR, history::Vector{Symbol})

Simplifies `fd` stepwise: finds the first matching rule for each term,
applies it once (→ one combination), recurses over the results until nothing
matches. Returns the final combination and the list of fired rule names.
Isomorphic (structurally AND label-equal) results are fused via
`circular_canonical_key(::CircularDecorated)` (CircularComboR). 1:1 analogue of
`reduce_graph`. After every rule step, the active
`circular_weight` termination assert runs (see above).
"""
function reduce_circular(fd::CircularDecorated; rules = CIRCULAR_RULES, maxsteps::Int = 1000)
    combo = CircularComboR(fd)
    history = Symbol[]
    for _ in 1:maxsteps
        fired = false
        next = CircularComboR()
        for (dd, coeff) in pairs_of(combo)
            hit = _circular_first_rule_match(dd, rules)
            if hit === nothing
                next = next + (coeff * CircularComboR(dd))
            else
                (name, res) = hit
                push!(history, name)
                fired = true
                next = next + (coeff * res)
            end
        end
        combo = next
        fired || break
    end
    return combo, history
end
reduce_circular(g::CircularGraph; rules = CIRCULAR_RULES, maxsteps::Int = 1000) =
    reduce_circular(circular_decorated(g); rules = rules, maxsteps = maxsteps)

"""
    reduce_circular_full(fd::CircularDecorated; rules = CIRCULAR_RULES,
                    memo = Dict{Any,CircularComboR}(), all_matches = false) -> CircularComboR
    reduce_circular_full(g::CircularGraph; rules, memo, all_matches = false) -> CircularComboR

Fully reduces `fd` to a combination of irreducible `CircularDecorated` terms,
memoized over `circular_canonical_key(fd::CircularDecorated)` (which includes the
labels — two structurally equal but differently decorated graphs are
different memo entries). Termination follows from every rule strictly
decreasing `circular_weight` (active assert, see above). With `all_matches = true`
(debug/tests), also reduces over EVERY applicable rule and checks that all
results agree — a first-match confluence check, analogous to
`_assert_confluent`. 1:1 analogue of `reduce_full`.
"""
function reduce_circular_full(fd::CircularDecorated; rules = CIRCULAR_RULES,
                          memo::Dict{Any,CircularComboR} = Dict{Any,CircularComboR}(),
                          all_matches::Bool = false)
    key = circular_canonical_key(fd)
    haskey(memo, key) && return memo[key]
    hit = _circular_first_rule_match(fd, rules)
    if hit === nothing
        result = CircularComboR(fd)                          # irreducible
    else
        (_, res) = hit
        result = CircularComboR()
        for (t, c) in pairs_of(res)
            result = result + c * reduce_circular_full(t; rules = rules, memo = memo)
        end
        if all_matches
            _assert_circular_confluent(fd, result, rules, memo)
        end
    end
    memo[key] = result
    return result
end
reduce_circular_full(g::CircularGraph; rules = CIRCULAR_RULES,
                memo::Dict{Any,CircularComboR} = Dict{Any,CircularComboR}(),
                all_matches::Bool = false) =
    reduce_circular_full(circular_decorated(g); rules = rules, memo = memo, all_matches = all_matches)

# Reduces over EVERY applicable rule (not just the first) and checks that all
# results agree. Analogue of `_assert_confluent`.
function _assert_circular_confluent(fd::CircularDecorated, first_result::CircularComboR, rules, memo)
    for r in rules
        res = r.apply(fd)
        res === nothing && continue
        _assert_circular_weight_drops(r.name, fd, res)
        alt = CircularComboR()
        for (t, c) in pairs_of(res)
            alt = alt + c * reduce_circular_full(t; rules = rules, memo = memo)
        end
        alt == first_result ||
            error("reduce_circular_full: not confluent at $(circular_canonical_key(fd)): " *
                  "rule $(r.name) gives a different result")
    end
end

"""
    circular_combo_by_weight(c::CircularComboR) -> Vector{Tuple{CircularDecorated,SoergelPoly}}

The (term, coefficient) pairs of `c`, sorted by DESCENDING `circular_weight` of the
representative (Julia tuple comparison, lexicographic) — most complex terms
first. Circular analogue of `combo_by_complexity`, here directly via the
termination-weight candidate `circular_weight` instead of separate polynomial
degree/node count.
"""
function circular_combo_by_weight(c::CircularComboR)
    ps = pairs_of(c)
    return sort(ps; by = p -> circular_weight(p[1]), rev = true)
end

"""
    reduce_circular_combo(c::CircularComboR; rules = CIRCULAR_RULES) -> CircularComboR

Fully reduces every term of `c` (one shared memo across all terms), in
`circular_combo_by_weight` order (most complex term first). 1:1 analogue of
`reduce_combo`.
"""
function reduce_circular_combo(c::CircularComboR; rules = CIRCULAR_RULES)
    memo = Dict{Any,CircularComboR}()
    out = CircularComboR()
    for (fd, coeff) in circular_combo_by_weight(c)
        out = out + coeff * reduce_circular_full(fd; rules = rules, memo = memo)
    end
    return out
end
