# rules/ZamoTermRules.jl — the two Zamolodchikov rules as a WHOLE-TERM
# replacement.
#
# WHAT THIS IS. `rules/ZamoRules.jl` (Z1) replaces SITED, at one place in the
# diagram. Here stands the coarser but sum-capable version: a rule matches a
# WHOLE term of a `CircularComboR` by its MORPHISM KEY (`_zamo_mkeyf` =
# `circular_canonical_key` after rotating the leaves onto cut1 — the exact match,
# not modulo rotation/swap13; why that is needed for termination: modulo
# swap13, Z2 would hit its own output and oscillate) and replaces
# it by the sum of its right-hand sides:
#
#     Z1:  rev(Zamo(8,1))  ->  Zamo(1,8)                 1 term
#     Z2:  rev(Zamo(9,2))  ->  Zamo(2,9) + X15 - Y15     3 terms
#
# DERIVED, NOT COPIED (the principle of ZamoRules.jl): the
# right-hand sides of Z2 are RECOMPUTED on first access —
# `X = rev(Zamo(1,2)) ∘ Zamo(1,8) ∘ Zamo(8,9)` and
# `Y = rev(Zamo(1,2)) ∘ rev(Zamo(8,1)) ∘ Zamo(8,9)` reduced under
# `CIRCULAR_RULES`, the terms identified by their morphism keys. If the rule set
# changes, X15/Y15 change with it; frozen fixtures could drift away from it. The
# price is compute time on the FIRST `zamo_term_rules()` (the two reductions;
# order of minutes) — hence lazy with a cache, and
# `zamo_term_rules(z2 = false)` returns Z1 alone without that cost.
#
# NOT in CIRCULAR_RULES (the `circular_weight` assert forbids 7 -> 7 nodes). The
# intended place is a step of its own in `reduce_to_circular_leave`, after
# P1/2parallel and before D4. Termination
# there: by the exact key matching no rule hits its own output (checked by
# `_zamo_term_rules_terminate` at build time), every rule fires at most
# finitely often per term, after which the weight-decreasing steps take over
# again.

"""
    ZamoTermRule(name, lhs, rhs)

A whole-term rule: if a term of a `CircularComboR` carries the morphism key
`lhs` (with respect to the `cut1` handed in), it is replaced by the sum of the
`rhs` graphs with their coefficients. The `rhs` graphs are stored NORMALISED
(leaves rotated onto `cut1 = 0`) and are rotated back when substituted.
"""
struct ZamoTermRule
    name::Symbol
    lhs::Any                            # morphism key of the left-hand side
    rhs::Vector{Tuple{CircularGraph,Int}}    # (normalised graph, coefficient)
end

# Morphism key: rotate the leaves so that cut1 lies on 0, then take the
# canonical key (the `mkeyf` convention from DLBasisWriter).
_zamo_mkeyf(g::CircularGraph, c1::Int) =
    circular_canonical_key(_dlb_rot_leaves(g, mod(c1, length(letters(g.word)))))
_zamo_norm(g::CircularGraph, c1::Int) =
    _dlb_rot_leaves(g, mod(c1, length(letters(g.word))))
_zamo_denorm(g::CircularGraph, c1::Int) =
    _dlb_rot_leaves(g, -mod(c1, length(letters(g.word))))

# ---- the rules, derived lazily ------------------------------------------------

const _ZAMO_TERM_RULES_Z1 = Base.RefValue{Any}(nothing)
const _ZAMO_TERM_RULES_Z2 = Base.RefValue{Any}(nothing)

function _zamo_build_z1()
    m = Zamo(1, 8)
    ZamoTermRule(:zamo_1_8, _zamo_mkeyf(circular(flip(Zamo(8, 1)).graph), flip(Zamo(8, 1)).cut1),
                 [(_zamo_norm(circular(m.graph), m.cut1), 1)])
end

function _zamo_build_z2()
    A = Zamo(1, 8); B = flip(Zamo(8, 1))
    X = compose(flip(Zamo(1, 2)), compose(A, Zamo(8, 9)))
    Y = compose(flip(Zamo(1, 2)), compose(B, Zamo(8, 9)))
    cX, _ = reduce_circular(circular_decorated(circular(X.graph)))
    cY, _ = reduce_circular(circular_decorated(circular(Y.graph)))
    kz = _zamo_mkeyf(circular(Zamo(2, 9).graph), Zamo(2, 9).cut1)
    kr = _zamo_mkeyf(circular(flip(Zamo(9, 2)).graph), flip(Zamo(9, 2)).cut1)
    tX = [t for (t, _) in pairs_of(cX)]; tY = [t for (t, _) in pairs_of(cY)]
    # X reduces to Zamo(2,9) + X15,
    # Y reduces to rev(Zamo(9,2)) + Y15 — two terms each, coefficients 1.
    X24 = [t for t in tX if _zamo_mkeyf(t.graph, X.cut1) == kz]
    X15 = [t for t in tX if _zamo_mkeyf(t.graph, X.cut1) != kz]
    Y15 = [t for t in tY if _zamo_mkeyf(t.graph, Y.cut1) != kr]
    (length(tX) == 2 && length(tY) == 2 &&
     length(X24) == 1 && length(X15) == 1 && length(Y15) == 1) || error(
        "zamo_term_rules: the X/Y reduction does not give two terms each, one " *
        "of them Zamo(2,9) resp. rev(Zamo(9,2)); Z2 cannot be derived. " *
        "Found: X -> $(length(tX)) terms, Y -> $(length(tY)) terms.")
    all(isone(c) for (_, c) in pairs_of(cX)) && all(isone(c) for (_, c) in pairs_of(cY)) ||
        error("zamo_term_rules: X/Y terms carry coefficients != 1")
    ZamoTermRule(:zamo_2_9, kr,
                 [(_zamo_norm(only(X24).graph, X.cut1),  1),
                  (_zamo_norm(only(X15).graph, X.cut1),  1),
                  (_zamo_norm(only(Y15).graph, Y.cut1), -1)])
end

# Termination check: no right-hand side may hit a left-hand one again.
function _zamo_term_rules_terminate(rules::Vector{ZamoTermRule})
    for R in rules, (g, _) in R.rhs
        any(S.lhs == _zamo_mkeyf(g, 0) for S in rules) && error(
            "zamo_term_rules: a right-hand side of $(R.name) carries a left " *
            "key again — the rule set does not terminate")
    end
    return rules
end

"""
    zamo_term_rules(; z2 = true) -> Vector{ZamoTermRule}

The whole-term rules, derived on first call and cached. With `z2 = true`
(default) Z2 as well — its right-hand sides cost the two X/Y reductions the
first time (minutes, see the file header); `z2 = false` returns Z1 alone,
without that cost.
"""
function zamo_term_rules(; z2::Bool = true)
    _ZAMO_TERM_RULES_Z1[] === nothing && (_ZAMO_TERM_RULES_Z1[] = _zamo_build_z1())
    z2 || return _zamo_term_rules_terminate([_ZAMO_TERM_RULES_Z1[]])
    _ZAMO_TERM_RULES_Z2[] === nothing && (_ZAMO_TERM_RULES_Z2[] = _zamo_build_z2())
    return _zamo_term_rules_terminate([_ZAMO_TERM_RULES_Z1[], _ZAMO_TERM_RULES_Z2[]])
end

# ---- the pass ----------------------------------------------------------------

"""
    zamo_pass(c::CircularComboR, cut1::Int; rules = zamo_term_rules())
        -> (CircularComboR, Vector{Symbol})

ONE pass over all terms of `c`: if a term carries the left morphism key of a
rule (with respect to `cut1`), it is replaced by that rule's right-hand sides
(coefficients multiply through); every rule fires at most once per term. Back
come the new combination and the names of the rules that fired (empty = nothing
matched).

Labels: the rule replaces WHOLE terms; decorations of the replaced term do not
come along (the right-hand sides are unlabelled). Terms with non-trivial region
labels are therefore NOT touched.
"""
# ---- the pipeline step ---------------------------------------------------------
#
# The Zamo rules go INTO the reduction pipeline
# — not into CIRCULAR_RULES (there the circular_weight assert forbids them,
# 7 ↦ 7 nodes), but as their own step in reduce_to_circular_leave: after the
# parallel-edge steps (P1/2parallel), before D4. Their termination argument: by the exact morphism-key matching a Zamo
# rule never hits its own output (checked at build time,
# `_zamo_term_rules_terminate`) — it can fire only finitely often per branch,
# after which the weight-decreasing steps take over again.
#
# COST PREFILTER. The Z2 derivation costs the two X/Y reductions the first time
# (~40 s warm). So that not every `reduce_to_circular_leave`
# call triggers them, the boundary word is checked FIRST: the two left-hand sides
# carry fixed 12-letter boundary words, obtainable WITHOUT any reduction from the
# fixtures. Only on a match are the rules (and hence possibly the derivation)
# touched at all.

const _ZAMO_LHS_WORDS = Base.RefValue{Any}(nothing)
function _zamo_lhs_words()
    _ZAMO_LHS_WORDS[] === nothing && (_ZAMO_LHS_WORDS[] =
        Set([join(letters(circular(flip(Zamo(8, 1)).graph).word)),
             join(letters(circular(flip(Zamo(9, 2)).graph).word))]))
    return _ZAMO_LHS_WORDS[]
end

"""
    circular_zamo_step(fd::CircularDecorated, cut1::Int) -> Union{Nothing, CircularComboR}

The Zamo step for `reduce_to_circular_leave`: if the WHOLE term `fd` carries the
left morphism key of Z1 or Z2 (with respect to `cut1`), the sum of the
right-hand sides comes back, otherwise `nothing`. Terms with non-trivial region
labels are not touched (the right-hand sides are unlabelled). The boundary-word
prefilter keeps the step free of cost for ordinary terms.
"""
function circular_zamo_step(fd::CircularDecorated, cut1::Int)
    all(isone, fd.region_labels) || return nothing
    join(letters(fd.graph.word)) in _zamo_lhs_words() || return nothing
    k = _zamo_mkeyf(fd.graph, cut1)
    rules = zamo_term_rules()
    r = findfirst(R -> R.lhs == k, rules)
    r === nothing && return nothing
    out = CircularComboR()
    for (g, s) in rules[r].rhs
        out = out + (s * CircularComboR(circular_decorated(_zamo_denorm(g, cut1))))
    end
    return out
end

# ---- the SITED version ---------------------------------------------------------
#
# `apply_zamo_rule` (rules/ZamoRules.jl) plants exactly ONE right-hand side and
# stays in the plain stack. Z2 needs the counterpart that plants THREE right-hand
# sides at the same outer connections and carries coefficients along — and since
# X15/Y15 are CIRCULAR graphs (CIRCULAR_RULES fixed points with merged nodes),
# the result necessarily lives in the circular stack: a `CircularComboR` over the
# converted host. `circular(::WordGraph)` makes that possible: the conversion
# keeps `word`, `edges`, all ports and node indices VERBATIM
# (circular/CircularConvert.jl) — so the plain match carries over verbatim.
#
# LEAF ASSIGNMENT. The right-hand sides are stored normalised to `cut1 = 0`;
# normalising the LHS fixture to the same convention makes the assignment the
# IDENTITY on morphism positions.
# So RHS leaf `k'` hangs on the outer connection at which the LHS leaf
# `mod1(k' + cut1_lhs, n)` hung (the same logic as `_zamo_leaf_alignment`, only
# via the normalisation instead of via two leaf lists).

"""
    apply_zamo_rule_combo(g::WordGraph, anchor::Int; pair = (2, 9), check = true)
        -> Union{Nothing, CircularComboR}

The SITED whole-term rule: if the LHS fixture `rev(Zamo(j,i))` of the pair
`pair` sits at `anchor` as a subdiagram (plain matcher `find_zamo_matches`), its
seven nodes are deleted and EVERY right-hand side of the corresponding
`ZamoTermRule` is planted at the same outer connections — the result is a
`CircularComboR` with one term per right-hand side (for Z2 three, coefficients
`+1, +1, −1`; for `pair = (1, 8)` one term — the sited Z1 in the circular
picture).

`nothing` if no match sits at `anchor`. The host `g` stays plain; the conversion
goes through `circular(g)` (ports/indices stay valid verbatim). **Do not build
this into a driver** — the intended place is the Zamo step in
`reduce_to_circular_leave` (file header).
"""
function apply_zamo_rule_combo(g::WordGraph, anchor::Int;
                               pair::Tuple{Int,Int} = (2, 9), check::Bool = true)
    ms = [m for m in find_zamo_matches(g; pair = pair) if m.anchor == anchor]
    isempty(ms) && return nothing
    m = ms[1]
    lhs, _ = zamo_rule_fixtures(pair[1], pair[2])
    nl = length(letters(lhs.graph.word))
    rulename = Symbol("zamo_$(pair[1])_$(pair[2])")
    R = only(r for r in zamo_term_rules(z2 = pair != (1, 8)) if r.name === rulename)

    host = circular(g)
    dead = Set{Int}(values(m.phi))
    drop = Set{Int}()
    for ni in dead, (ei, _, _) in _circular_edges_at_node(host, ni)
        push!(drop, ei)
    end

    out = CircularComboR()
    for (Rg, s) in R.rhs
        hostplus = CircularGraph(host.word, vcat(host.nodes, Rg.nodes), host.edges)
        off = length(host.nodes)
        keep = Edge[e for (j, e) in enumerate(hostplus.edges) if !(j in drop)]
        np(p::NodePort) = NodePort(off + p.node, p.slot)
        function np(p::Leaf)
            k = mod1(p.k + lhs.cut1, nl)        # normalised RHS leaf ↦ LHS leaf
            haskey(m.ext, k) || error(
                "apply_zamo_rule_combo: no outer connection for LHS leaf $k")
            return m.ext[k]
        end
        np(p::Circle) = p
        for e in Rg.edges
            push!(keep, Edge(e.colour, np(e.a), np(e.b)))
        end
        term = _circular_delete_nodes(hostplus, dead, keep)
        if check
            v = check_wiring(term)
            isempty(v) || error("apply_zamo_rule_combo: term (coeff $s) violates " *
                                "the wiring: $v")
        end
        out = out + (s * CircularComboR(circular_decorated(term)))
    end
    return out
end

"""
    apply_zamo_rule_combo(m::MorphismGraph, anchor::Int; kwargs...)

The same rule on a morphism; the cuts stay where they are. Back comes
`(combo::CircularComboR, cut1, cut2)` — the combination carries no cuts, so they
are handed back alongside.
"""
function apply_zamo_rule_combo(m::MorphismGraph, anchor::Int; kwargs...)
    out = apply_zamo_rule_combo(m.graph, anchor; kwargs...)
    out === nothing && return nothing
    return (out, m.cut1, m.cut2)
end

# ---- the pass -----------------------------------------------------------------

function zamo_pass(c::CircularComboR, cut1::Int; rules::Vector{ZamoTermRule} = zamo_term_rules())
    out = CircularComboR(); fired = Symbol[]
    for (dd, coeff) in pairs_of(c)
        r = all(isone, dd.region_labels) ?
            findfirst(R -> R.lhs == _zamo_mkeyf(dd.graph, cut1), rules) : nothing
        if r === nothing
            out = out + (coeff * CircularComboR(dd))
        else
            push!(fired, rules[r].name)
            for (g, s) in rules[r].rhs
                out = out + (coeff * (s * CircularComboR(circular_decorated(_zamo_denorm(g, cut1)))))
            end
        end
    end
    return out, fired
end
