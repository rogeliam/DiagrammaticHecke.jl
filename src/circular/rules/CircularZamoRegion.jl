# circular/rules/CircularZamoRegion.jl — the Zamolodchikov rule as a MULTI-REGION
# pattern.
#
# WHAT THIS IS. All region rules in `circular/rules/CircularRegionRules.jl` match on ONE
# region. Zamo is the first rule whose site is a pattern spanning MULTIPLE
# regions: it is characterized by its INNER REGIONS rather than by its nodes.
#
# THE COUNT: **6** inner regions. All four Zamo sides
# (`flip(Zamo(8,1))`, `Zamo(1,8)`, `flip(Zamo(9,2))`, `Zamo(2,9)`) have 7
# nodes, 24 edges, 18 regions, of which **6 are inner, all with boundary-word
# length 3**; their corners cover all 7 nodes. The six regions therefore
# determine the cluster completely.
#
# THE COUNT 6 IS ONLY A PRE-FILTER, NOT A GUARD. Both pairs have 6, all four
# sides have 6, and `Zamo(1,1)` (the full 14-step cycle) has 17 inner regions
# of length 3. The actual guard is the COLOURED ISOMORPHISM of the region
# pattern:
#
#   Zamo(1,8):  4-cycle with two pendants, cycle colours 1/3, pendant colour 2
#   Zamo(2,9):  HEXAGON, alternating colours (LHS 2/3, RHS 1/2)
#
# — i.e. TWO DIFFERENT PATTERNS, but one shared machinery parameterized by the
# pair.
#
# DERIVED, NOT HARD-CODED (same principle as ZamoRules.jl):
# the pattern is COMPUTED from the LHS fixture, not written down. If
# `Zamolodchikov.jl` changes, the pattern changes with it.
#
# SURGERY = DELEGATION, the same way `_frr_braid_back` delegates to
# `_fr_braid_back_at`: the right-hand side comes from `apply_zamo_rule_combo`
# (rules/ZamoTermRules.jl, the plain version), so it is identical to the
# oracle by construction; the region matcher's own job is to find the SAME
# sites as `find_zamo_matches`.
#
# WHERE IT RUNS: as the `:zamo` step of `reduce_to_circular_leave`, after
# P1/2parallel and before D4, switched by `CIRCULAR_ZAMO_ENABLED` (on by
# default) — like `circular_zamo_step`, and NOT in `CIRCULAR_RULES`. Zamo is a
# RELATION, not a reduction: `circular_weight` does not fall (7 ↦ 7 nodes, Z2
# even goes 1 term ↦ 3 terms), and the assert in CircularDriver.jl forbids
# that. The registry is the structural rules alone; a driver step is where a
# weight-neutral relation belongs.
# The termination measure there is not `circular_weight`: the rule runs in BOTH
# directions (Z1⁻¹/Z2⁻¹), so the measure is the direction rule below — backward
# only at ≥ 7 `zamo_regions` and a blocked forward step, accepted only if a
# non-Zamo rule then fires (which then carries the measure), plus the round-1
# guard `after_inverse` in `reduce_to_circular_leave`. Forward alone stays free
# of self-hits.
#
# PLANARITY: the rule reads the face tracer; on non-planar wiring faces are
# meaningless and the rule simply doesn't match (conservative, like all
# region rules).
#
# MATCHING UP TO ROTATION, NEVER REFLECTION — as in `circular_canonical_key`
# and CircularRegionRules.jl.
#
# ---- PLACE IN THE PIPELINE -------------------------------------------------
#
# The rule is hooked into `reduce_to_circular_leave` (CircularDecoratedRules.jl,
# after circular_zamo_step/P1, before D4; pre-filter `all(isone, labels) && nodes >= 7`).
# The rotation guard `colour_ok` (section "COLOUR SEQUENCE" below) breaks a
# three-step cycle between P1 and ZamoRegion; see that section for the
# mechanism.
#
# An end term can carry more than 6 inner regions.
#
# The runtime cost is the price of the region tracer, which
# `circular_zamo_region_matches` calls per term; the pipeline's node-count
# pre-filter keeps it away from small terms.

# ---- the pattern -----------------------------------------------------------

"""
    ZamoRegionPattern

The region pattern of a Zamo LHS fixture, derived from it (not hard-coded).
The pattern regions are `1:length(colours)` (the fixture's inner regions in
tracer order); per pattern region `p`:

* `colours[p][j]` — colour of the `j`-th boundary letter,
* `edges[p][j]`   — edge index in the fixture,
* `corners[p][j]` — corner node (fixture node index) after the `j`-th letter,
* `slots[p][j]`   — `(slot_in, slot_out)` at this corner node (stage-2
  surgery, see [`circular_zamo_region_matches`](@ref)).

Also `share[fixture_edge] = the pattern regions bordering it` (the
propagation structure), `nodes` = all cluster nodes of the fixture, `anchors`
= the anchor candidates where the PLAIN rule (`find_zamo_matches`) attaches,
`outer` = for each of the fixture's 12 outer leaves the triple
`(leaf_number, fixture_node, fixture_slot)`, and `degree[fixture_node] =
degree` (stage-2 surgery, `Zamo-Region` extension).
"""
struct ZamoRegionPattern
    pair::Tuple{Int,Int}
    len::Int
    colours::Vector{Vector{Int}}
    edges::Vector{Vector{Int}}
    corners::Vector{Vector{Int}}
    slots::Vector{Vector{Tuple{Int,Int}}}
    share::Dict{Int,Vector{Int}}
    nodes::Vector{Int}
    anchors::Vector{Int}
    outer::Vector{Tuple{Int,Int,Int}}
    degree::Dict{Int,Int}
    nodearms::Dict{Int,Vector{Int}}
end

Base.show(io::IO, P::ZamoRegionPattern) = print(io,
    "ZamoRegionPattern(", P.pair, ": ", length(P.colours), " inner regions of ",
    P.len, ", ", length(P.nodes), " nodes, anchors ", P.anchors, ")")

# The inner regions of a CircularGraph, as (region, letter list) — empty boundary
# words (free circles, Faces.jl) drop out.
function _fzr_inner(g::CircularGraph)
    out = Tuple{Int,Vector{CircularBoundaryLetter}}[]
    for w in circular_boundary_words(g)
        isempty(w.letters) && continue
        is_interior(w) || continue
        push!(out, (w.region, w.letters))
    end
    return out
end

const _FZR_PATTERNS = Dict{Tuple{Int,Int},ZamoRegionPattern}()
# Patterns for the BACKWARD direction (Z1⁻¹/Z2⁻¹): the same derivation, only
# from the OTHER fixture. Own cache, so `pair` can stay the key.
const _FZR_PATTERNS_INV = Dict{Tuple{Int,Int},ZamoRegionPattern}()

"""
    zamo_region_pattern(pair = (1, 8); inverse = false) -> ZamoRegionPattern

The region pattern of the LHS fixture `flip(Zamo(j,i))` for `pair = (i,j)`,
derived from the fixture and cached on first call. **Cheap** — it does NOT
touch `zamo_term_rules()`; the expensive Z2 derivation only happens in the
surgery.

`inverse = true` gives the pattern of the **backward direction** (Z1⁻¹/Z2⁻¹):
matches `Zamo(i,j)` instead of `flip(Zamo(j,i))`. `zamo_rule_fixtures` returns
both sides and the derivation is the same; only the cache is separate
(`_FZR_PATTERNS_INV`), so `pair` can stay the key.

Throws if the fixture does not have the expected shape (inner regions all the
same length, corners covering all nodes).
"""
function zamo_region_pattern(pair::Tuple{Int,Int} = (1, 8); inverse::Bool = false)
    cache = inverse ? _FZR_PATTERNS_INV : _FZR_PATTERNS
    haskey(cache, pair) && return cache[pair]
    fwd_lhs, fwd_rhs = zamo_rule_fixtures(pair[1], pair[2])
    # Backward direction: the two sides swap roles. `zamo_rule_fixtures`
    # returns `(flip(Zamo(j,i)), Zamo(i,j))`; for Z1⁻¹/Z2⁻¹ the left-hand
    # side to match on is therefore `Zamo(i,j)`.
    lhs = inverse ? fwd_rhs : fwd_lhs
    F = circular(lhs.graph)                       # node/edge indices VERBATIM
    inner = _fzr_inner(F)
    isempty(inner) && error(
        "zamo_region_pattern$(pair)$(inverse ? " (inverse)" : ""): the LHS fixture has no inner regions")
    len = length(inner[1][2])
    all(length(ls) == len for (_, ls) in inner) || error(
        "zamo_region_pattern$(pair)$(inverse ? " (inverse)" : ""): inner regions of differing boundary length " *
        "$(sort(unique(length(ls) for (_, ls) in inner))) — fixture has changed, " *
        "the pattern needs re-checking")

    colours = [[l.colour for l in ls] for (_, ls) in inner]
    edges   = [[l.edge   for l in ls] for (_, ls) in inner]
    corners = [[l.corner.vertex[2] for l in ls] for (_, ls) in inner]
    slots   = [[(l.corner.slot_in, l.corner.slot_out) for l in ls] for (_, ls) in inner]

    # The fixture's outer leaves: edges with exactly one leaf end. For the
    # stage-2 surgery (CircularGraph-native splicing), each leaf's adjacent
    # fixture node+slot must be known.
    outer = Tuple{Int,Int,Int}[]
    for e in F.edges
        aleaf, bleaf = e.a isa Leaf, e.b isa Leaf
        aleaf == bleaf && continue          # both leaf or both node port: internal/irrelevant
        k, p = aleaf ? (e.a.k, e.b) : (e.b.k, e.a)
        push!(outer, (k, p.node, p.slot))
    end
    length(outer) == 12 || error(
        "zamo_region_pattern$(pair)$(inverse ? " (inverse)" : ""): $(length(outer)) outer leaves instead of the " *
        "expected 12 — fixture has changed, the pattern needs re-checking")

    share = Dict{Int,Vector{Int}}()
    for (p, es) in enumerate(edges), e in es
        v = get!(share, e, Int[])
        p in v || push!(v, p)
    end
    nodes = sort(unique(vcat(corners...)))
    length(nodes) == length(F.nodes) || error(
        "zamo_region_pattern$(pair)$(inverse ? " (inverse)" : ""): the inner regions cover only " *
        "$(length(nodes)) of $(length(F.nodes)) nodes — the cluster is not " *
        "determined by them, the pattern does not carry")

    # Anchor candidates: the same choice as `find_zamo_matches` (strict
    # roles, otherwise "any 23-braid").
    strict = zamo_anchors(lhs.graph)
    anchors = isempty(strict) ? _zamo_braid23_anchors(lhs.graph) : strict
    isempty(anchors) && error(
        "zamo_region_pattern$(pair)$(inverse ? " (inverse)" : ""): the LHS fixture contains no 23-braid — " *
        "the rule is not formulated for this pair")

    degree = Dict{Int,Int}(v => arm_count(F.nodes[v]) for v in nodes)
    nodearms = Dict{Int,Vector{Int}}(v => arms(F.nodes[v]) for v in nodes)

    P = ZamoRegionPattern(pair, len, colours, edges, corners, slots, share, nodes,
                          anchors, outer, degree, nodearms)
    cache[pair] = P
    return P
end

# ---- the matcher -------------------------------------------------------------
#
# Propagation like `_zamo_cluster_match` (rules/ZamoRules.jl), just one
# dimension up: over REGIONS and their shared edges instead of nodes and their
# slots. From a seed (pattern region 1 ↦ host region, rotation) everything
# else is forced.

# One attempt from a seed. `hw[R]` = letters of inner host region `R` (only
# matching length), `hshare[e]` = list `(R, j)` of occurrences of edge `e` in
# those boundary words, `hostnodes` = `g.nodes` of the host (for the degree
# check of the stage-2 surgery: only cluster nodes whose degree exactly
# matches their fixture counterpart (`P.degree`) are accepted — genuine
# foreign-connection merges return `nothing`).
function _fzr_try(P::ZamoRegionPattern, hw::Dict{Int,Vector{CircularBoundaryLetter}},
                  hshare::Dict{Int,Vector{Tuple{Int,Int}}}, hr0::Int, r0::Int,
                  hostnodes::Vector{CircularNode})
    n = P.len
    regmap  = Dict{Int,Tuple{Int,Int}}()
    usedreg = Set{Int}()
    nodemap = Dict{Int,Int}(); nodeinv = Dict{Int,Int}()
    edgemap = Dict{Int,Int}(); edgeinv = Dict{Int,Int}()
    slotoffset = Dict{Int,Int}(); badslot = Set{Int}()
    queue = Tuple{Int,Int,Int}[(1, hr0, r0)]

    while !isempty(queue)
        (p, hr, r) = popfirst!(queue)
        if haskey(regmap, p)
            regmap[p] == (hr, r) || return nothing
            continue
        end
        hr in usedreg && return nothing
        haskey(hw, hr) || return nothing
        ls = hw[hr]
        all(j -> P.colours[p][j] == ls[mod1(j + r, n)].colour, 1:n) || return nothing
        regmap[p] = (hr, r); push!(usedreg, hr)

        for j in 1:n
            hl = ls[mod1(j + r, n)]
            fe, fv = P.edges[p][j], P.corners[p][j]
            he, hv = hl.edge, hl.corner.vertex[2]
            fsi, fso = P.slots[p][j]
            hsi, hso = hl.corner.slot_in, hl.corner.slot_out

            if haskey(edgemap, fe)
                edgemap[fe] == he || return nothing
            else
                get(edgeinv, he, fe) == fe || return nothing
                edgemap[fe] = he; edgeinv[he] = fe
            end
            if haskey(nodemap, fv)
                nodemap[fv] == hv || return nothing
            else
                get(nodeinv, hv, fv) == fv || return nothing
                nodemap[fv] = hv; nodeinv[hv] = fv
            end

            # Stage 2: degree, COLOUR SEQUENCE and slot rotation per cluster
            # node — NOT decisive for the match itself (the region match
            # stays degree-independent, see the docstring of
            # `circular_zamo_region_matches`). On a mismatch (wrong degree,
            # different colour sequence despite equal degree, inconsistent
            # offset across several regions of the same node), `fv` lands in
            # `badslot` and never re-enters `slotoffset` — the surgery
            # detects the gap and cleanly declines for this match
            # (`circular_zamo_region_surgery`: `length(match.outer_ports) == 12 ||
            # return nothing`).
            #
            # COLOUR SEQUENCE. The region match alone only sees the INNER
            # boundary words, not the node shape itself — a `[1,1,3,3]` node
            # (formed by P1 merging two parallel edges) carries exactly the
            # same colour MULTISET as the alternating `[1,3,1,3]` pattern of
            # the real fixture (`zamo_rule_fixtures(1,8)` nodes 3/7), but a
            # DIFFERENT cyclic order. Without this check the region matcher
            # would still match, the surgery would insert Zamo, and the
            # result would again be a P1 candidate — P1 and Zamo-Region would
            # then fire in a non-terminating three-step cycle
            # (`reduce_to_circular_leave` recurses on every rule step, no
            # explicit loop, so this would overflow the stack). Hence: the
            # arms of the host node, starting from the entry slot `hsi` (in
            # host rotation), must be ROTATION-EQUAL to the fixture's colour
            # sequence starting from `fsi` — not merely carry the same
            # multiset.
            if fv ∉ badslot
                deg = P.degree[fv]
                fixarms = P.nodearms[fv]
                hostarms = arms(hostnodes[hv])
                colour_ok = deg == arm_count(hostnodes[hv]) &&
                    all(t -> hostarms[mod1(hsi + t - fsi, deg)] == fixarms[t], 1:deg)
                off_in, off_out = mod(hsi - fsi, deg), mod(hso - fso, deg)
                ok = colour_ok && off_in == off_out &&
                     (!haskey(slotoffset, fv) || slotoffset[fv] == off_in)
                if ok
                    slotoffset[fv] = off_in
                else
                    delete!(slotoffset, fv); push!(badslot, fv)
                end
            end

            for q in P.share[fe]
                (q == p || haskey(regmap, q)) && continue
                jq = findfirst(==(fe), P.edges[q])
                jq === nothing && return nothing
                # the OTHER host region at this edge; ambiguous ⇒ reject
                cands = [c for c in get(hshare, he, Tuple{Int,Int}[]) if c[1] != hr]
                length(cands) == 1 || return nothing
                Rq, jjq = cands[1]
                Rq in usedreg && return nothing
                push!(queue, (q, Rq, mod(jjq - jq, n)))
            end
        end
    end

    length(regmap) == length(P.colours) || return nothing
    length(nodemap) == length(P.nodes)  || return nothing
    # `slotoffset` is allowed to be INCOMPLETE (degree/rotation mismatch at
    # individual nodes) — that is not a match failure, see above.
    return (regmap = regmap, nodemap = nodemap, edgemap = edgemap,
           slotoffset = slotoffset)
end

# For each outer fixture leaf `(k, fv, fs)` from `P.outer`: the host port
# `NodePort(nodemap[fv], slot)` is the CONNECTION OF THE CLUSTER NODE ITSELF
# (which gets deleted during surgery) — not the external end where the new
# RHS edge must attach. This function finds, in the host `g`, the edge
# hanging off that cluster port and returns its OTHER end. `nothing` in the
# result dict for `k` if no such edge exists (an unwired slot), or if both
# ends lie inside the cluster (two outer fixture leaves directly
# cross-connected — then the external connection is not uniquely located
# outside the cluster, a stage-2 non-case).
function _fzr_outer_ports(g::CircularGraph, P::ZamoRegionPattern,
                          nodemap::Dict{Int,Int}, slotoffset::Dict{Int,Int})
    dead = Set(values(nodemap))
    out = Dict{Int,Port}()
    for (k, fv, fs) in P.outer
        inner = NodePort(nodemap[fv], mod1(fs + slotoffset[fv], P.degree[fv]))
        for e in g.edges
            other = e.a == inner ? e.b : (e.b == inner ? e.a : nothing)
            other === nothing && continue
            (other isa NodePort && other.node in dead) && continue
            out[k] = other
            break
        end
    end
    return out
end

"""
    circular_zamo_region_matches(g::CircularGraph; pair = (1, 8)) -> Vector{NamedTuple}
    circular_zamo_region_matches(g::WordGraph; pair = (1, 8))
    circular_zamo_region_matches(m::MorphismGraph; pair = (1, 8))

All sites in `g` where the region pattern of the LHS fixture for `pair` sits.
Each hit is

    (regions, nodemap, edgemap, anchors, slotoffset, outer_ports)

with `regions` = the sorted host regions, `nodemap[fixture_node] =
host_node`, `edgemap` analogously, `anchors` = the host nodes where the plain
surgery can attach, `slotoffset[fixture_node]` = the slot rotation between
fixture and host node (defined only if both have the same degree — stage-2
surgery), and `outer_ports[leaf_number] = NodePort(host_node,
host_slot)` = the fixture's 12 outer connections in the host (empty if any
cluster node has a degree different from the fixture — then `slotoffset`
carries no valid value there and the stage-2 surgery cannot attach, see
`circular_zamo_region_surgery`).

**The region match itself checks NO node shape** — neither `kind` nor slot
colours. That is exactly what makes it invariant under merging (and
potentially WIDER than the plain matcher; a hit with no plain counterpart is
a finding). The DEGREE is checked by the slot
rotation guard, but only to populate `slotoffset`/`outer_ports` — a
hit's `regions`/`nodemap`/`edgemap`/`anchors` still comes back even on a
degree mismatch, just without `outer_ports`.

Deduplicated over the set of host regions. The pre-filter is the sheer count
(`< 6 inner regions of matching length ⇒ out`) — it is not a guard on its
own; the coloured isomorphism provides that.
"""
function circular_zamo_region_matches(g::CircularGraph; pair::Tuple{Int,Int} = (1, 8),
                                 inverse::Bool = false)
    P = zamo_region_pattern(pair; inverse = inverse)
    inner = [(R, ls) for (R, ls) in _fzr_inner(g) if length(ls) == P.len]
    length(inner) >= length(P.colours) || return NamedTuple[]

    hw = Dict{Int,Vector{CircularBoundaryLetter}}(R => ls for (R, ls) in inner)
    hshare = Dict{Int,Vector{Tuple{Int,Int}}}()
    for (R, ls) in inner, (j, l) in enumerate(ls)
        push!(get!(hshare, l.edge, Tuple{Int,Int}[]), (R, j))
    end

    out = NamedTuple[]; seen = Set{Vector{Int}}()
    for (R, _) in inner, r in 0:(P.len - 1)
        m = _fzr_try(P, hw, hshare, R, r, g.nodes)
        m === nothing && continue
        regs = sort([hr for (hr, _) in values(m.regmap)])
        regs in seen && continue
        push!(seen, regs)
        outer_ports = if length(m.slotoffset) == length(P.nodes)
            _fzr_outer_ports(g, P, m.nodemap, m.slotoffset)
        else
            Dict{Int,NodePort}()   # at least one cluster node has a different degree
        end
        push!(out, (regions = regs, nodemap = m.nodemap, edgemap = m.edgemap,
                    anchors = [m.nodemap[a] for a in P.anchors],
                    slotoffset = m.slotoffset, outer_ports = outer_ports))
    end
    return out
end

circular_zamo_region_matches(g::WordGraph; kwargs...) =
    circular_zamo_region_matches(circular(g); kwargs...)
circular_zamo_region_matches(m::MorphismGraph; kwargs...) =
    circular_zamo_region_matches(circular(m.graph); kwargs...)

