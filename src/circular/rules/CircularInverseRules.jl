# circular/rules/CircularInverseRules.jl — the CIRCULAR rules RUN BACKWARDS.
#
# Same house conventions:
#   * every function is named `expand_circular_*` and takes the SITE as an argument,
#   * none of them belongs in any registry — NEVER add to `CIRCULAR_RULES`, it would
#     grow the diagram and the driver would stop terminating,
#   * the backward direction is DERIVED: build a preimage, run the FORWARD rule
#     on it, its result supplies coefficients and correction terms.
# The full rationale (including which rules are NOT invertible and why) is not
# repeated here — only the circular-specific points are.
#
# CIRCULAR SPECIFICS: a node has no fixed arm count. That gives an inverse here that
# cannot exist in the plain stack — SPLITTING a circular node (`expand_circular_merge`,
# the reverse of C5 `_fr_merge`). This is exactly why the circular type exists:
# merging forgets how a cluster was bracketed, so splitting has to be TOLD the
# bracketing.

# ---- generic: solve for one term --------------------------------------------
#
# These three helpers are type-agnostic: they work on the `terms` dict that
# `DiagramCombo`, `DiagramComboR` and `CircularCombo` all have.

# The inverse of a coefficient, or `nothing` if it is not a unit. Int: only ±1
# (there are no other units in ℤ). SoergelPoly: a nonzero CONSTANT (degree 0);
# anything containing an α is not a unit in R.
_coeff_inverse(c::Integer) = (c == 1 || c == -1) ? c : nothing
function _coeff_inverse(c::SoergelPoly)
    iszero(c) && return nothing
    degree(c) == 0 || return nothing
    q = constant_term(c)
    return soergel_monomial((0, 0, 0), inv(q))
end

"""
    _solve_for_term(fwd, preimage_combo, tgtkey) -> Union{Nothing, typeof(fwd)}

Solve the forward identity `G = Σ_j c_j T_j` for the term with key `tgtkey`:

    T_k = c_k⁻¹ · ( G − Σ_{j≠k} c_j T_j )

`fwd` is the result of the forward rule (the right-hand side), `preimage_combo`
the preimage `G` as a one-element combination of the same type. `nothing` if
`tgtkey` does not occur at all or its coefficient is not a unit (then the
equation is not solvable over ℤ resp. over R).
"""
function _solve_for_term(fwd, preimage_combo, tgtkey)
    haskey(fwd.terms, tgtkey) || return nothing
    ck = fwd.terms[tgtkey][2]
    inv_ck = _coeff_inverse(ck)
    inv_ck === nothing && return nothing
    rest = fwd - _combo_of_key(fwd, tgtkey)        # Σ_{j≠k} c_j T_j
    return inv_ck * (preimage_combo - rest)
end

# The one-element combination `c_k T_k` out of `c`, given its key.
function _combo_of_key(c::C, key) where {C}
    out = C(typeof(c.terms)())
    haskey(c.terms, key) && (out.terms[key] = c.terms[key])
    return out
end

"""
    _check_expansion(fwd_result, target_key, where_) -> nothing

Self-check: the input diagram has to come back out under the forward rule.
Throws with an explanation if it does not — see the first-match trap below.
"""
function _check_expansion(fwd, tgtkey, where_::AbstractString)
    fwd === nothing && error(
        "$where_: the forward rule does not apply to the constructed preimage AT " *
        "ALL — the construction does not match the rule pattern (bug here, not in " *
        "the caller)")
    haskey(fwd.terms, tgtkey) || error(
        "$where_: the forward rule applies but does not return the input diagram. " *
        "Most common cause: the diagram contains a SECOND site where the same rule " *
        "matches earlier (first match). Force with `check = false` and look for " *
        "yourself.")
    return nothing
end

# ---- C5⁻¹: split a circular node back apart --------------------------------------
#
# Forward (`merge_at_edge`, circular/CircularMerge.jl): from A (degree `d_a`, connection
# at slot `s_a`) and B (degree `d_b`, connection at slot `s_b`) a single node
# with arm sequence
#     [A from s_a+1 forward (d_a-1 entries), B from s_b+1 forward (d_b-1 entries)] .
# Backward is a decomposition of the cyclic arm sequence into two arcs: A gets
# `[c, arms[k1..k2]]`, B gets the rest, with `c` the colour of the new
# connecting edge.
#
# NOT WLOG SLOT 1. `CircularNode` equality only holds up to rotation, but
# `circular_canonical_key` does NOT: rotating a single node's slot numbering (same
# cyclic arm sequence, same planarity, same region count) changes the key
# (e.g. on `circular(join_trivalents(1))`: `r = 0` gives the same key, `r = 1` and
# `r = 2` do not). That is why `conn_slot_a`/`conn_slot_b` are arguments — with them the
# round trip closes (2 of 144 combinations hit the original), without them
# never.

"""
    expand_circular_merge(g::CircularGraph, ni::Int, k1::Int, k2::Int, colour::Int;
                     conn_slot_a = 1, conn_slot_b = 1, check = true)
        -> Union{Nothing, CircularGraph}

Reverse of C5 (`_fr_merge`): split circular node `ni` into TWO nodes, joined by a
fresh edge of colour `colour`. The arms in the cyclic interval `k1..k2` go to
one node, all remaining ones to the other.

This is the genuine circular operation backwards: `merge` forgets HOW a cluster was
bracketed (exactly what the type is built for) — splitting has to be told the
bracketing. `k1`, `k2` and `colour` are that instruction; nothing is guessed.

`conn_slot_a`/`conn_slot_b` say on WHICH slot of the two new nodes the
connection sits. That is not a mathematical choice — the cyclic arm sequence
is the same — but `circular_canonical_key` sees the difference: rotating one node's
slot numbering changes the key even though arm sequence, planarity and region
count stay the same (e.g. on `circular(join_trivalents(1))`).
Anyone who wants a particular preimage back — e.g. in a round-trip test — must
specify the slots.

`nothing` if `ni` is not a valid node, the interval is empty or everything, or
either arm sequence does not yield a valid `CircularNode` (`circular_node` forbids
2-armed nodes, the "beads").

With `check = true` (default), checks that `merge_at_edge` on the new edge
gives back exactly `g` — the round trip closes.
"""
function expand_circular_merge(g::CircularGraph, ni::Int, k1::Int, k2::Int, colour::Int;
                          conn_slot_a::Int = 1, conn_slot_b::Int = 1,
                          check::Bool = true)
    (1 <= ni <= length(g.nodes)) || return nothing
    nd = g.nodes[ni]
    _circular_braidlike(nd) && return nothing          # braid-like nodes (incl. gbraid) are not mergeable
    d = arm_count(nd)
    (1 <= k1 <= d && 1 <= k2 <= d) || return nothing

    # the arm slots of the first arc, cyclically from k1 to k2
    firstarc = Int[]
    k = k1
    while true
        push!(firstarc, k)
        k == k2 && break
        k = mod1(k + 1, d)
        length(firstarc) > d && return nothing
    end
    secondarc = [s for s in 1:d if !(s in firstarc)]
    (isempty(firstarc) || isempty(secondarc)) && return nothing

    # arm sequence of a new node: the connection at `cs`, then the arc
    # clockwise. `slot_of[j]` = slot of the j-th arc arm.
    function build(arc::Vector{Int}, cs::Int)
        dnew = length(arc) + 1
        (1 <= cs <= dnew) || return (nothing, Int[])
        a = zeros(Int, dnew)
        a[cs] = colour
        slots = Int[]
        for (j, sold) in enumerate(arc)
            sl = mod1(cs + j, dnew)
            a[sl] = arm_colour(nd, sold)
            push!(slots, sl)
        end
        return (a, slots)
    end
    armsA, slotsA = build(firstarc, conn_slot_a)
    armsB, slotsB = build(secondarc, conn_slot_b)
    (armsA === nothing || armsB === nothing) && return nothing
    A = try; circular_node(armsA); catch; return nothing; end
    B = try; circular_node(armsB); catch; return nothing; end

    # Slot map: old slot -> (new node, new slot).
    newport = Dict{Int,Tuple{Int,Int}}()
    gplus = CircularGraph(g.word, vcat(g.nodes, [A, B]), g.edges)
    Ai, Bi = length(gplus.nodes) - 1, length(gplus.nodes)
    for (pos, s) in enumerate(firstarc);  newport[s] = (Ai, slotsA[pos]); end
    for (pos, s) in enumerate(secondarc); newport[s] = (Bi, slotsB[pos]); end

    drop = Set{Int}()
    keep = Edge[]
    for (ei, e) in enumerate(g.edges)
        touches = (e.a isa NodePort && e.a.node == ni) ||
                  (e.b isa NodePort && e.b.node == ni)
        touches || continue
        push!(drop, ei)
    end
    for (ei, e) in enumerate(gplus.edges)
        ei in drop && continue
        push!(keep, e)
    end
    remap(p::Port) = p
    remap(p::NodePort) = p.node == ni ?
        (haskey(newport, p.slot) ? NodePort(newport[p.slot]...) : p) : p
    for ei in sort(collect(drop))
        e = g.edges[ei]
        push!(keep, Edge(e.colour, remap(e.a), remap(e.b)))
    end
    conn = Edge(colour, NodePort(Ai, conn_slot_a), NodePort(Bi, conn_slot_b))
    push!(keep, conn)
    G = _circular_delete_nodes(gplus, Set([ni]), keep)

    if check
        # the connecting edge is the one pushed last
        back = try
            merge_at_edge(G, length(G.edges))
        catch err
            error("expand_circular_merge: cannot merge the new connecting edge back " *
                  "($err) — the split does not match `merge_at_edge`.")
        end
        circular_canonical_key(back) == circular_canonical_key(g) || error(
            "expand_circular_merge: merging back does NOT reproduce the input. " *
            "Construction bug here, not in the caller.")
    end
    return G
end

# ---- planarity helpers, shared by the expansions below ---------------------

# Euler test for CircularGraph: an ALIAS for the public `is_planar_embedding`
# (diagram/WiringCheck.jl, generic for CircularGraph too). The private name
# stays for existing call sites.
_circular_is_planar_embedding(g::CircularGraph) = is_planar_embedding(g)

_circular_planar_first(cands::Vector{CircularGraph}) =
    vcat([G for G in cands if _circular_is_planar_embedding(G)],
         [G for G in cands if !_circular_is_planar_embedding(G)])

# ---- C8⁻¹: ONE braid + trivalent -> TWO braids + trivalent -----------------
#
# Circular version of R9⁻¹ (`expand_braid_relation`) — the full derivation is
# not repeated here: which roles the seven outer connections carry,
# why the reading direction (2·2) and the rotation of the three preimage nodes
# (3·3·3) are enumerated, why planarity only sorts. In the
# gallery this is F2-3 `braid_relation`, forward `_fr_braid_relation`
# (circular/rules/CircularBraidRules.jl).
#
# Circular specifics:
#   * an m=3 braid is a 6-armed `:braid` node with alternating arm sequence
#     (how `circular_node` builds it), a pure trivalent a `[q,q,q]`
#     (`_circular_is_pure_trivalent`);
#   * `_fr_braid_relation` requires — like the plain forward rule — a PURE
#     trivalent as the third node. Generalising to arbitrary circular nodes is a
#     separate open point, not implemented here.

"""
    expand_circular_braid_relation(g::CircularGraph, tv::Int, bi::Int;
                              variant = 1, check = true) -> Union{Nothing, CircularGraph}

Reverse of C8 (`_fr_braid_relation`, gallery F2-3) on `CircularGraph`: from the
AFTER cluster (pure trivalent `tv` of colour `p` + 6-armed `:braid` `bi` with
colour pair `{p,q}`), rebuild the BEFORE cluster with two braids and a
trivalent of colour `q`. The seven outer connections may be boundary leaves OR
ports of other nodes.

AMBIGUOUS, just like the plain stack (`expand_braid_relation`): the forward
rule matches by role, not slot. `variant` picks one (planar first),
`expand_circular_braid_relation_variants` returns the whole list.
"""
function expand_circular_braid_relation(g::CircularGraph, tv::Int, bi::Int;
                                   variant::Int = 1, check::Bool = true)
    cands = expand_circular_braid_relation_variants(g, tv, bi)
    if isempty(cands)
        check && error(
            "expand_circular_braid_relation: no reading/rotation gives a preimage that " *
            "maps back to the input under C8 — the pattern does not match (pure " *
            "trivalent + 6-armed :braid as in C8's AFTER). A second C8 site in " *
            "the diagram is NOT a reason: the round trip is " *
            "checked SITED (`_fr_braid_relation_at`).")
        return nothing
    end
    (1 <= variant <= length(cands)) || throw(ArgumentError(
        "variant = $variant, but there are $(length(cands)) preimages"))
    return cands[variant]
end

"""
    expand_circular_braid_relation_variants(g::CircularGraph, tv::Int, bi::Int) -> Vector{CircularGraph}

All preimages for `expand_circular_braid_relation` — planar ones first, empty if
the pattern does not match.
"""
function expand_circular_braid_relation_variants(g::CircularGraph, tv::Int, bi::Int)
    none = CircularGraph[]
    (1 <= tv <= length(g.nodes) && 1 <= bi <= length(g.nodes)) || return none
    tv == bi && return none
    nt, nb = g.nodes[tv], g.nodes[bi]
    # NOTE: `_circular_is_pure_trivalent` requires `kind === :mixed` and so misses a
    # single-colour 2-trivalent — `circular_node([2,2,2])` is `:mono`. But C8's
    # AFTER trivalent carries colour p, which may well be 2. Hence the arm
    # condition directly here.
    (arm_count(nt) == 3 && length(unique(arms(nt))) == 1) || return none
    (nb.kind === :braid && arm_count(nb) == 6) || return none
    p = arms(nt)[1]
    p in arms(nb) || return none
    q = p == arms(nb)[1] ? arms(nb)[2] : arms(nb)[1]
    p == q && return none

    tslot = Dict{Int,Port}(); bslot = Dict{Int,Port}()
    for (_, pa, po) in _circular_edges_at_node(g, tv)
        pa isa NodePort || return none
        tslot[pa.slot] = po
    end
    for (_, pa, po) in _circular_edges_at_node(g, bi)
        pa isa NodePort || return none
        bslot[pa.slot] = po
    end
    (length(tslot) == 3 && length(bslot) == 6) || return none
    m = [k for k in 1:3 if tslot[k] isa NodePort && tslot[k].node == bi]
    n = [k for k in 1:6 if bslot[k] isa NodePort && bslot[k].node == tv]
    (length(m) == 1 && length(n) == 1) || return none
    m, n = m[1], n[1]
    tslot[m].slot == n && bslot[n].slot == m || return none
    arm_colour(nb, n) == p || return none

    tgt = circular_canonical_key(g)
    cands = CircularGraph[]
    # The round trip is checked SITED, not by first-match: the forward rule
    # only returns the FIRST hit in the whole diagram, and once the preimage
    # has a SECOND C8 site, first-match could hit the wrong one and drop a
    # valid candidate (e.g. `general_braid_cluster(3)`: depth-2
    # enumeration gives 36 instead of 57 preimages). `_fr_braid_relation_at`
    # (circular/rules/CircularBraidRules.jl) checks exactly the triple
    # `_circular_c8_preimage` just built.
    for dn in (1, -1), dt in (1, -1)
        B_mid  = bslot[mod1(n + dn, 6)]
        B_near = bslot[mod1(n + 2dn, 6)]
        e_T    = bslot[mod1(n + 3dn, 6)]
        S_near = bslot[mod1(n + 4dn, 6)]
        S_mid  = bslot[mod1(n + 5dn, 6)]
        B_far  = tslot[mod1(m + dt, 3)]
        S_far  = tslot[mod1(m - dt, 3)]
        for rS in (0, 2, 4), rB in (0, 2, 4), rt in 0:2
            G = _circular_c8_preimage(g, tv, bi, p, q,
                                 S_far, S_mid, S_near, e_T, B_near, B_mid, B_far,
                                 rS, rB, rt)
            G === nothing && continue
            S, T, B = _circular_c8_new_site(G)
            fwd = _fr_braid_relation_at(G, S, B, T)
            fwd !== nothing && haskey(fwd.terms, tgt) && push!(cands, G)
        end
    end
    return _circular_planar_first(cands)
end

# The three fresh nodes of `_circular_c8_preimage` — (S, T, B) — always end up at
# the back of the preimage: they are appended to `g.nodes` (slots m+1, m+2,
# m+3), and `_circular_delete_nodes` only strips the two old nodes `tv`, `bi`
# before them, keeping the order of the rest. That is how the round-trip
# check knows which site to run the forward rule on.
_circular_c8_new_site(G::CircularGraph) =
    (length(G.nodes) - 2, length(G.nodes) - 1, length(G.nodes))

# One rotation variant of the C8 preimage (path A); cf. `_r9_preimage`.
function _circular_c8_preimage(g::CircularGraph, tv::Int, bi::Int, p::Int, q::Int,
                          S_far::Port, S_mid::Port, S_near::Port, e_T::Port,
                          B_near::Port, B_mid::Port, B_far::Port,
                          rS::Int, rB::Int, rt::Int)
    gplus = CircularGraph(g.word,
                     vcat(g.nodes, [circular_node([p, q, p, q, p, q]),
                                    circular_node([q, q, q]),
                                    circular_node([q, p, q, p, q, p])]),
                     g.edges)
    S, T, B = length(gplus.nodes) - 2, length(gplus.nodes) - 1, length(gplus.nodes)
    dead = Set([tv, bi])
    drop = Set{Int}()
    for ni in (tv, bi), (ei, _, _) in _circular_edges_at_node(g, ni)
        push!(drop, ei)
    end
    keep = Edge[e for (j, e) in enumerate(gplus.edges) if !(j in drop)]
    sp(i) = NodePort(S, mod1(i + rS, 6))
    bp(i) = NodePort(B, mod1(i + rB, 6))
    tp(i) = NodePort(T, mod1(i + rt, 3))
    push!(keep, Edge(p, S_far,  sp(1)))
    push!(keep, Edge(q, S_mid,  sp(2)))
    push!(keep, Edge(p, S_near, sp(3)))
    push!(keep, Edge(q, sp(4), tp(1)))
    push!(keep, Edge(q, e_T,   tp(2)))
    push!(keep, Edge(q, sp(6), bp(1)))
    push!(keep, Edge(p, sp(5), bp(2)))
    push!(keep, Edge(q, tp(3), bp(3)))
    push!(keep, Edge(p, B_near, bp(4)))
    push!(keep, Edge(q, B_mid,  bp(5)))
    push!(keep, Edge(p, B_far,  bp(6)))
    return _circular_delete_nodes(gplus, dead, keep)
end

# ---- D4⁻¹ in the circular stack: the 1 -> 3 example -----------------------------

"""
    expand_circular_d4(fd::CircularDecorated, tv::Int, dot_arm::Int;
                  preimage_labels = nothing, check = true) -> Union{Nothing, CircularComboR}

Reverse of `apply_circular_d4`, solved for the TRIVALENT term — the circular version of
`expand_d4`, where the derivation is:

    T_triv = T_dotB + T_dotC - D

`tv` is a 3-armed circular node `[i,i,i]` in `fd`, `dot_arm ∈ 1:3` says which of its
arms was the DOT strand in the preimage; the other two are connected there by
edge `ei`.

LABELS as in the plain stack: `apply_circular_d4` OVERWRITES the label of the BC
region with `alpha(i)`, so the old one cannot be reconstructed — the preimage
is built unlabelled unless `preimage_labels` says otherwise.
"""
function expand_circular_d4(fd::CircularDecorated, tv::Int, dot_arm::Int;
                       preimage_labels::Union{Nothing,Vector{SoergelPoly}} = nothing,
                       check::Bool = true)
    g = fd.graph
    (1 <= tv <= length(g.nodes)) || return nothing
    nt = g.nodes[tv]
    _circular_is_pure_trivalent(nt) || return nothing
    (1 <= dot_arm <= 3) || return nothing
    at = _circular_edges_at_node(g, tv)
    length(at) == 3 || return nothing
    i = arms(nt)[1]
    byslot = Dict(pt.slot => po for (_, pt, po) in at)
    all(haskey(byslot, k) for k in 1:3) || return nothing
    others = [k for k in 1:3 if k != dot_arm]

    gplus = CircularGraph(g.word, vcat(g.nodes, [circular_node([i])]), g.edges)
    dt = length(gplus.nodes)
    drop = Set(ei for (ei, _, _) in at)
    keep = Edge[e for (j, e) in enumerate(gplus.edges) if !(j in drop)]
    push!(keep, Edge(i, byslot[dot_arm], NodePort(dt, 1)))            # dot strand
    push!(keep, Edge(i, byslot[others[1]], byslot[others[2]]))        # connection ei
    G = _circular_delete_nodes(gplus, Set([tv]), keep)
    new_dot = length(G.nodes)
    new_ei = length(G.edges)

    labels = preimage_labels === nothing ?
             fill(one(SoergelPoly), region_count(G)) : preimage_labels
    # Carry the TARGET's outer label along: `tgt` below
    # is `circular_canonical_key(fd)` and includes it — a preimage with `outer = 1`
    # could otherwise never hit the target term.
    D = CircularDecorated(G, labels, fd.outer_label)
    fwd = apply_circular_d4(D, new_dot, new_ei)
    tgt = circular_canonical_key(fd)
    if check
        _check_expansion(fwd, tgt, "expand_circular_d4")
    elseif fwd === nothing
        return nothing
    end
    return _solve_for_term(fwd, CircularComboR(D), tgt)
end
