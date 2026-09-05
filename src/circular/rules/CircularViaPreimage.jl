# circular/rules/CircularViaPreimage.jl — the combinator `circular_via_preimage`,
# the body of the scheme also used by C18 `_fr_gen12_two_edges`
# (CircularGen12TwoEdges.jl).
#
# THE SCHEME, in three steps: recognize an EDGE pattern (the
# site), locally EXPAND it (reduce⁻¹), let the existing rules run on it via
# `reduce_circular` (WITHOUT itself), and accept the result only if it (a) is not
# empty, (b) is STRICTLY lighter than the input (`circular_weight`, the
# termination argument), and (c) carries no foreign labels (the intermediate
# reduction runs undecorated). C18 is the only rule that does exactly this;
# this combinator is its body.
#
# WHY HERE AND NOT FURTHER UP. `circular_via_preimage` needs `reduce_circular` and
# `CIRCULAR_RULES` (CircularDriver.jl) as well as `circular_weight` (CircularWeight.jl), both of
# which appear textually LATER in the include list. That's harmless: the
# reference sits inside a function body and is only resolved at RUNTIME
# (Julia resolves global names in a function body only when it's called, not
# when the method is defined) — the same forward reference as
# `expand_circular_merge` in CircularBraidChannels.jl. This file sits right
# AFTER CircularGen12Expand.jl (which supplies the expansion tools, though
# this file doesn't need them itself — it's node-kind-agnostic) and BEFORE
# CircularGen12TwoEdges.jl, which uses it. Future rules following the same
# scheme (C6/C7/C8) line up after this file too.
#
# NO RECURSION on its own — that's the caller's job: `except_rule` must be the
# caller's own rule name, and the combinator then filters it out of the
# intermediate reduction, exactly as C18 does.
#
# ⚡ `except_rule` IS OPTIONAL (`nothing` = run the full
# registry), for the generic preimage stage in the driver
# (`_fr_generic_preimage`, CircularDriver.jl). The reason is a side
# effect of the filter: `except_rule` goes by RULE NAME, and a rule name in the
# registry always covers the whole entry — for the braid rules that includes
# the non-circular path of the same rule too (`:braid_relation` =
# `_fr_braid_relation_gen`, which internally tries the classical
# `_fr_braid_relation` first). A caller protecting a single rule via
# `except_rule = <its own name>` thus also disables, in the intermediate
# reduction, exactly the ordinary path it needs after expanding. The
# combinator is therefore unsuitable for single-rule extensions of C6/C7/C8;
# the generic stage doesn't have this problem because it is not a registry
# rule at all and gets its recursion guard elsewhere (a depth counter plus
# registry identity in `_circular_first_rule_match`, see CircularDriver.jl).

"""
    circular_via_preimage(g::CircularGraph, sites, expansions::Function;
                      except_rule::Union{Nothing,Symbol} = nothing,
                      w0 = circular_weight(g)) -> Union{Nothing, CircularComboR}

Scheme combinator, also used by C18: for each site `v ∈ sites`, try each
expansion `e ∈ expansions(v)` — `expansions` returns, for a site, the
(already checked to be planar) expansions as `Vector{CircularGraph}`, e.g.
`_circular_gen12_expansions`. For each expansion,
`reduce_circular(circular_decorated(e); rules = CIRCULAR_RULES minus `except_rule`)`
is computed, and the result is accepted ONLY if

  (a) it is not empty (`CircularComboR()` means "is zero";
      `isempty(pairs_of(c))` is discarded like a non-match, the same as in
      C18),
  (b) EVERY result term has `circular_weight` STRICTLY below `w0` (the
      termination argument, also required by the `CircularDriver.jl` assert —
      checked here up front so the assert never fires and the rule cleanly
      declines instead), and
  (c) NO result term carries non-trivial labels (the intermediate reduction
      runs on an UNDECORATED graph; if its result carries labels, they
      aren't the caller's).

The FIRST combination (site, then expansion, in the order supplied by
`sites`/`expansions`) satisfying all three conditions is taken; its terms are
copied UNDECORATED (just the `CircularGraph`s, with coefficient) into a new
`CircularComboR` — the caller (`_lift` in CircularDriver.jl) assigns labels afterward.

`nothing` if `sites` is empty or no site/expansion combination satisfies all
three conditions.

`except_rule = nothing` (the default) runs the **full** `CIRCULAR_RULES` — this is the
case for the generic preimage stage (`_fr_generic_preimage`, CircularDriver.jl),
which is not itself a registry rule and gets its recursion guard elsewhere.
A caller that IS ITSELF in `CIRCULAR_RULES` MUST pass its own rule name; doing so
also excludes the non-circular path of that same rule from the intermediate
reduction (see the file header).

These three conditions are taken verbatim from `_fr_gen12_two_edges` (not
invented, not generalized) — see CircularGen12TwoEdges.jl for the derivation.
"""
function circular_via_preimage(g::CircularGraph, sites, expansions::Function;
                           except_rule::Union{Nothing,Symbol} = nothing,
                           w0 = circular_weight(g))
    isempty(sites) && return nothing
    # ALWAYS A COPY, even when `except_rule === nothing`:
    # `_circular_first_rule_match` (CircularDriver.jl) only appends the generic
    # preimage stage when `rules === CIRCULAR_RULES` IDENTICALLY. A fresh array
    # here means the intermediate reduction never sees that stage — the
    # second guard against recursion, alongside the depth counter.
    other_rules = [r for r in CIRCULAR_RULES if except_rule === nothing || r.name !== except_rule]

    for v in sites
        for e in expansions(v)
            c, _ = reduce_circular(circular_decorated(e); rules = other_rules)
            isempty(pairs_of(c)) && continue
            # All terms must be strictly lighter — otherwise the expansion
            # was just a detour (typically: the rule's own merge rebuilds
            # the site immediately).
            all(circular_weight(t.graph) < w0 for (t, _) in pairs_of(c)) || continue
            # The intermediate reduction ran on an UNDECORATED graph; if its
            # result carries labels, they wouldn't be the caller's. Better
            # not to apply it then.
            all(all(isone, t.region_labels) for (t, _) in pairs_of(c)) || continue

            out = CircularComboR()
            for (t, coeff) in pairs_of(c)
                _add!(out, t.graph, coeff)
            end
            return out
        end
    end
    return nothing
end
