# circular/rules/CircularZamoRegionStep.jl — the plain WordGraph-delegated step and the
# CircularGraph-native stage-2 surgery for the Zamo region rule (deletes the cluster
# nodes directly and splices in the right-hand side via outer_ports).
#
# Included directly after CircularZamoRegion.jl.

# ---- the step (surgery by delegation, plain) --------------------------------

"""
    circular_zamo_region_step(g::WordGraph; pair = (1, 8), check = true)
        -> Union{Nothing, CircularComboR}

The rule: find the site via the region pattern, insert the right-hand side
via the PLAIN implementation `apply_zamo_rule_combo` (one term for `(1,8)`,
three with coefficients `+1, +1, −1` for `(2,9)`). `nothing` if no site is
found OR the plain surgery does not fire at any anchor candidate.

**The second case is not an error**: there the region matcher saw a site whose
seven nodes do not sit in the fixture shape (merged) — for
`WordGraph` hosts this means the cluster cannot be expressed as a plain
sub-diagram at all. For `CircularGraph` hosts (terms from the reduction
pipeline) there is the CircularGraph-native version
[`circular_zamo_region_step(::CircularGraph)`](@ref) — see there.

⚠️ **Not in `CIRCULAR_RULES`** — see the file header. Its place is a step in
`reduce_to_circular_leave`, next to `circular_zamo_step`
(`circular/rules/CircularDecoratedRules.jl`, in the `zam === nothing` branch).
"""
function circular_zamo_region_step(g::WordGraph; pair::Tuple{Int,Int} = (1, 8),
                              check::Bool = true)
    for m in circular_zamo_region_matches(g; pair = pair), a in m.anchors
        out = apply_zamo_rule_combo(g, a; pair = pair, check = check)
        out === nothing || return out
    end
    return nothing
end

# ---- the surgery (CircularGraph-native, stage 2, restricted case) ---------------
#
# Unlike `apply_zamo_rule_combo` (plain, ZamoTermRules.jl), this version never
# goes through a `WordGraph`: it deletes the cluster nodes directly on the
# `CircularGraph` and splices in the right-hand side via the connections in
# `match.outer_ports` (instead of via `m.ext`). This requires that EVERY one
# of the 7 cluster nodes in the host has exactly the degree of its fixture
# counterpart — only then does `circular_zamo_region_matches` return a full
# `outer_ports` (12 entries). Otherwise (a genuine foreign-connection merge,
# e.g. three fixture nodes collapsed into one host node), this function
# returns `nothing` — a cleanly reported non-case, not a splitting surgery in
# the `circular_via_preimage` style.

"""
    zamo_inverse_terms(pair) -> (lhs::MorphismGraph, rhs::Vector{Tuple{CircularGraph,Int}})

The **backward direction** of the Zamo relation for `(i,j)`, COMPUTED from
the forward rule (same principle as the file header: derived, not hard-coded).
The forward rule `zamo_term_rules()` has the shape

    flip(Zamo(j,i))  ↦  1·Zamo(i,j)  +  Σ s_k · T_k

— the term carrying the Zamo key is the **main term** and has coefficient 1.
Solved for it:

    Zamo(i,j)  ↦  1·flip(Zamo(j,i))  −  Σ s_k · T_k

Concretely, for Z1 this is `Zamo(1,8) ↦ rev(Zamo(8,1))` (one term) and for Z2
`Zamo(2,9) ↦ rev(Zamo(9,2)) − X15 + Y15` (three terms, signs `+1,−1,+1`). Both
come out of the same computation here — the signs are nowhere hard-coded.

Throws if the forward rule does not have exactly ONE main term with
coefficient 1; then `zamo_term_rules` has changed and the inverse needs
re-checking.
"""
function zamo_inverse_terms(pair::Tuple{Int,Int})
    fwd_lhs, fwd_rhs = zamo_rule_fixtures(pair[1], pair[2])
    rulename = Symbol("zamo_$(pair[1])_$(pair[2])")
    R = only(r for r in zamo_term_rules(z2 = pair != (1, 8)) if r.name === rulename)

    # The main term: the one whose normalized key is that of `Zamo(i,j)`.
    kz = _zamo_mkeyf(_zamo_norm(circular(fwd_rhs.graph), fwd_rhs.cut1), 0)
    main = findall(t -> _zamo_mkeyf(t[1], 0) == kz, R.rhs)
    length(main) == 1 || error(
        "zamo_inverse_terms$(pair): $(length(main)) right-hand sides carry the " *
        "Zamo$(pair) key instead of exactly one — the forward rule has " *
        "changed, the inverse needs re-checking")
    h = only(main)
    R.rhs[h][2] == 1 || error(
        "zamo_inverse_terms$(pair): the main term has coefficient " *
        "$(R.rhs[h][2]) instead of 1 — division would be needed, not supported")

    inv_rhs = Tuple{CircularGraph,Int}[(_zamo_norm(circular(fwd_lhs.graph), fwd_lhs.cut1), 1)]
    for (k, (Rg, s)) in enumerate(R.rhs)
        k == h && continue
        push!(inv_rhs, (Rg, -s))
    end
    return (fwd_rhs, inv_rhs)
end

"""
    circular_zamo_region_surgery(g::CircularGraph, match; pair = (1, 8), check = true)
        -> Union{Nothing, CircularComboR}

The CircularGraph-native surgery for ONE hit `match` from
`circular_zamo_region_matches(g; pair)`. Deletes the 7 cluster nodes
(`values(match.nodemap)`) and splices in each right-hand-side term of the
associated `ZamoTermRule`, with the outer connections from
`match.outer_ports`.

`nothing` if `match.outer_ports` does not cover all 12 leaves (at least one
cluster node has a degree different from its fixture counterpart — the
restricted case handled by this stage does not apply, see the file header)
or the 7 host nodes in `match.nodemap` are not pairwise distinct.
"""
function circular_zamo_region_surgery(g::CircularGraph, match;
                                 pair::Tuple{Int,Int} = (1, 8), check::Bool = true,
                                 inverse::Bool = false)
    length(match.outer_ports) == 12 || return nothing
    dead = Set{Int}(values(match.nodemap))
    length(dead) == length(match.nodemap) || return nothing   # 7 pairwise distinct

    if inverse
        lhs, rhs_terms = zamo_inverse_terms(pair)
    else
        lhs, _ = zamo_rule_fixtures(pair[1], pair[2])
        rulename = Symbol("zamo_$(pair[1])_$(pair[2])")
        rhs_terms = only(r for r in zamo_term_rules(z2 = pair != (1, 8))
                         if r.name === rulename).rhs
    end
    nl = length(letters(lhs.graph.word))

    drop = Set{Int}()
    for ni in dead, (ei, _, _) in _circular_edges_at_node(g, ni)
        push!(drop, ei)
    end

    out = CircularComboR()
    for (Rg, s) in rhs_terms
        hostplus = CircularGraph(g.word, vcat(g.nodes, Rg.nodes), g.edges)
        off = length(g.nodes)
        keep = Edge[e for (j, e) in enumerate(hostplus.edges) if !(j in drop)]
        np(p::NodePort) = NodePort(off + p.node, p.slot)
        function np(p::Leaf)
            k = mod1(p.k + lhs.cut1, nl)        # normalized RHS leaf ↦ LHS leaf
            haskey(match.outer_ports, k) || error(
                "circular_zamo_region_surgery: no outer connection for LHS leaf $k")
            return match.outer_ports[k]
        end
        np(p::Circle) = p
        for e in Rg.edges
            push!(keep, Edge(e.colour, np(e.a), np(e.b)))
        end
        term = _circular_delete_nodes(hostplus, dead, keep)
        if check
            v = check_wiring(term)
            isempty(v) || error("circular_zamo_region_surgery: term (coeff $s) violates " *
                                "the wiring: $v")
        end
        out = out + (s * CircularComboR(circular_decorated(term)))
    end
    return out
end

"""
    circular_zamo_region_step(g::CircularGraph; pair = (1, 8), check = true)
        -> Union{Nothing, CircularComboR}

The CircularGraph-native version of [`circular_zamo_region_step(::WordGraph)`](@ref):
searches via [`circular_zamo_region_matches`](@ref), splices via
[`circular_zamo_region_surgery`](@ref) — both without a `WordGraph` detour.
`nothing` if no site is found OR the surgery does not fire on any hit (a
degree mismatch at at least one cluster node, see the header of the surgery
section — a restricted non-case, not an error).
"""
function circular_zamo_region_step(g::CircularGraph; pair::Tuple{Int,Int} = (1, 8),
                              check::Bool = true, inverse::Bool = false)
    for m in circular_zamo_region_matches(g; pair = pair, inverse = inverse)
        out = circular_zamo_region_surgery(g, m; pair = pair, check = check,
                                      inverse = inverse)
        out === nothing || return out
    end
    return nothing
end

"""
    circular_zamo_region_step(m::MorphismGraph; kwargs...)

The same rule on a morphism; the cuts stay where they are. Returns
`(combo::CircularComboR, cut1, cut2)` — the combination itself carries no cuts, so
they are passed alongside (as with `apply_zamo_rule_combo`).
"""
function circular_zamo_region_step(m::MorphismGraph; kwargs...)
    out = circular_zamo_region_step(m.graph; kwargs...)
    out === nothing && return nothing
    return (out, m.cut1, m.cut2)
end
