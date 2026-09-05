# Zamo fixtures and anchors on WordGraph — used by circular/rules/ZamoTermRules.jl and
# CircularZamoRegion.jl. The Zamolodchikov relation as a CALLABLE rule.
#
# WHAT THIS IS. The Zamolodchikov relation `rev(Zamo(8,1)) = Zamo(1,8)`
# exists nowhere else as a rule: `src/morphism/Zamolodchikov.jl` only builds
# the two sides, and `zamo_relation_endpoints_agree()` merely checks the
# endpoints. Zamo never occurs in any rule registry.
#
# ⚠️ NEVER put this into a rule registry. This is **hand tooling**: every function
# takes the SITE (the anchor) as an argument, so it does not search by itself, and
# it stands in NO table. `reduce_circular` does not see this file. Reasons:
#   1. The rule does NOT lower the node count (7 ↦ 7). A driver applying it does
#      not terminate — it can be applied back and forth arbitrarily often.
#   2. It is an ASSERTED equality: the existing local rules do not
#      prove it, it is input, not consequence. In a driver it would therefore be a
#      silent extra assumption.
#
# CONSTRUCTION: DERIVED FROM THE FIXTURE, NOT REBUILT. The rule does not
# state the pattern a second time. It takes
#
#     LHS fixture = flip(Zamo(j,i))      ("rev")
#     RHS fixture = Zamo(i,j)
#
# — both built by `Zamolodchikov.jl` — looks for a structural isomorphism of the
# LHS fixture in the input diagram and PLANTS the RHS fixture onto the same outer
# connections. If `Zamolodchikov.jl` changes, the rule changes with it; a hand-wired
# pattern could drift apart from it, this one cannot.
#
# THE LEAF ASSIGNMENT OF THE TWO FIXTURES. `flip` SWAPS THE CUTS: `Zamo(1,8)` has
# `cuts = (0,6)` (bottom = leaves 1..6), `flip(Zamo(8,1))` has `cuts = (6,0)`
# (bottom = leaves 7..12). So the leaf numbers of the two fixtures do NOT agree.
# The assignment therefore goes by MORPHISM POSITION: k-th bottom leaf ↦ k-th
# bottom leaf, k-th top leaf ↦ k-th top leaf (`_zamo_leaf_alignment`). That is the
# only place where the rule has to know something about the fixtures that is not in
# their edge list.

# ---- the fixtures ------------------------------------------------------------

"""
    zamo_rule_fixtures(i::Int, j::Int) -> (lhs::MorphismGraph, rhs::MorphismGraph)

The two sides of the Zamolodchikov rule for the antipodal index pair `(i, j)`
with `j = i + 7 (mod 14)`: `lhs = flip(Zamo(j, i))` (the side that is REPLACED)
and `rhs = Zamo(i, j)` (the side that is inserted). Both are 7-node morphisms with
the same endpoints; the rule is `lhs → rhs`.

Other `(i, j)` are NOT rejected (the machinery is general), but only for the
antipodal pairs do both sides have 7 nodes.
"""
zamo_rule_fixtures(i::Int, j::Int) = (flip(Zamo(j, i)), Zamo(i, j))

# leaf k of the LHS fixture ↦ leaf of the RHS fixture, by morphism position
# (bottom against bottom, top against top). Throws if the boundary lengths differ.
function _zamo_leaf_alignment(lhs::MorphismGraph, rhs::MorphismGraph)
    lb, lt = _bottom_leaves(lhs), _top_leaves(lhs)
    rb, rt = _bottom_leaves(rhs), _top_leaves(rhs)
    (length(lb) == length(rb) && length(lt) == length(rt)) || error(
        "zamo rule: the two fixtures have different boundary lengths " *
        "(bottom $(length(lb))/$(length(rb)), top $(length(lt))/$(length(rt)))")
    align = Dict{Int,Int}()
    for (a, b) in zip(lb, rb); align[a] = b; end
    for (a, b) in zip(lt, rt); align[a] = b; end
    return align
end

# ---- the conditions at the anchor --------------------------------------------
#
# Go through all 23-braids. Conditions: all
# three 3-arms go into 13-nodes; at least two 2-arms go into 12-nodes; the
# MIDDLE 13-node (the one between the two 12-connections) is wired to
# everyone; plus a further 23-node, everything wired together as in the
# (1,8) fixture — 7 nodes and 12 edges total.
#
# What stands here is the PREFILTER (the enumerable conditions). The condition
# "everything wired together as in the fixture" is not an enumeration but exactly a
# structural isomorphism — `_zamo_cluster_match` below does that.

"Is `nd` a 23-braid (m=3, colours {2,3})?"
_is_23_braid(nd::Node) = nd.kind === :braid && nd.m == 3 && Set(nd.colours) == Set((2, 3))
"Is `nd` a 12-braid (m=3, colours {1,2})?"
_is_12_braid(nd::Node) = nd.kind === :braid && nd.m == 3 && Set(nd.colours) == Set((1, 2))
"Is `nd` a 13-crossing (m=2, colours {1,3})?"
_is_13_cross(nd::Node) = nd.kind === :braid && nd.m == 2 && Set(nd.colours) == Set((1, 3))

# slot ↦ far port at node `ni` (only if every slot is occupied exactly once).
function _slot_map(g::WordGraph, ni::Int)
    slot = Dict{Int,Port}()
    for (_, pa, po) in _edges_at_node(g, ni)
        pa isa NodePort || return nothing
        haskey(slot, pa.slot) && return nothing
        slot[pa.slot] = po
    end
    length(slot) == _node_degree(g.nodes[ni]) || return nothing
    return slot
end

"""
    zamo_anchor_roles(g::WordGraph, A::Int) -> Union{Nothing, NamedTuple}

The conditions at the 23-braid `A`, checked point by point. Returns
`nothing` if one of them is violated, otherwise the roles found

    (anchor, crosses, twelves, middle, other23, cluster)

* `crosses`  — the three 13-nodes at the three 3-arms of `A`,
* `twelves`  — the (at least two) 12-nodes at 2-arms of `A`,
* `middle`   — the 13-node at the 3-arm BETWEEN the two 12-connections,
* `other23`  — the further 23-braid that `middle` sees,
* `cluster`  — all seven nodes.

This is the PREFILTER; whether the wiring really matches the fixture is decided
only by `find_zamo_matches`.
"""
function zamo_anchor_roles(g::WordGraph, A::Int)
    (1 <= A <= length(g.nodes)) || return nothing
    _is_23_braid(g.nodes[A]) || return nothing
    slot = _slot_map(g, A)
    slot === nothing && return nothing
    nd = g.nodes[A]
    s3 = [k for k in 1:6 if _slot_colour(nd, k) == 3]
    s2 = [k for k in 1:6 if _slot_colour(nd, k) == 2]

    # (a) all three 3-arms go into 13-nodes
    crosses = Int[]
    for k in s3
        p = slot[k]
        p isa NodePort || return nothing
        _is_13_cross(g.nodes[p.node]) || return nothing
        push!(crosses, p.node)
    end
    length(unique(crosses)) == 3 || return nothing

    # (b) at least two 2-arms go into 12-nodes
    tw_slots = [k for k in s2 if slot[k] isa NodePort && _is_12_braid(g.nodes[slot[k].node])]
    length(tw_slots) >= 2 || return nothing
    twelves = unique([slot[k].node for k in tw_slots])
    length(twelves) >= 2 || return nothing

    # (c) the MIDDLE 13-node: its 3-arm lies cyclically between two
    #     12-connections. Between two 2-slots at distance 2 there is exactly one
    #     3-slot; we take the pair for which that works out.
    mid = 0; midpair = (0, 0)
    for a in tw_slots, b in tw_slots
        a == b && continue
        mod(b - a, 6) == 2 || continue
        m = mod1(a + 1, 6)
        p = slot[m]
        p isa NodePort || continue
        mid = p.node; midpair = (a, b); break
    end
    mid == 0 && return nothing

    # (d) the middle 13-node is "wired to everyone": to the anchor, to BOTH
    #     12-nodes of the pair and to a FURTHER 23-braid.
    mslot = _slot_map(g, mid)
    mslot === nothing && return nothing
    nb = [p.node for p in values(mslot) if p isa NodePort]
    A in nb || return nothing
    n_a, n_b = slot[midpair[1]].node, slot[midpair[2]].node
    (n_a in nb && n_b in nb) || return nothing
    others = [x for x in nb if x != A && _is_23_braid(g.nodes[x])]
    isempty(others) && return nothing
    other23 = others[1]

    cluster = unique(vcat([A], crosses, [n_a, n_b], [other23]))
    length(cluster) == 7 || return nothing
    return (anchor = A, crosses = crosses, twelves = [n_a, n_b],
            middle = mid, other23 = other23, cluster = cluster)
end

"""
    zamo_anchors(g::WordGraph) -> Vector{Int}

All nodes of `g` satisfying the anchor conditions (`zamo_anchor_roles`).
"""
zamo_anchors(g::WordGraph) = [ni for ni in 1:length(g.nodes)
                              if zamo_anchor_roles(g, ni) !== nothing]

"The weak prefilter for the pair (2,9): all 23-braids."
_zamo_braid23_anchors(g::WordGraph) =
    [ni for (ni, nd) in enumerate(g.nodes) if _is_23_braid(nd)]

"""
    zamo_prefilter_kind(pair::Tuple{Int,Int}) -> Symbol

Which prefilter applies for this pair: `:roles` (the strict conditions from
`zamo_anchor_roles`, as for the pair `(1,8)`) or `:braid23` (only "is a 23-braid",
the condition for the pair `(2,9)`).
"""
zamo_prefilter_kind(pair::Tuple{Int,Int}) =
    isempty(zamo_anchors(zamo_rule_fixtures(pair[1], pair[2])[1].graph)) ? :braid23 : :roles

# ---- the actual matcher: isomorphism to the fixture ---------------------------
#
# Starting from an anchor pair (fixture node `fa` ↦ graph node `ga`) with slot
# offset `d`, the assignment is propagated ALONG THE EDGES. The slot offset is
# needed because the slot NUMBERING of a node is only cyclically meaningful
# (diagram/Graph.jl §PLANAR SLOT CONVENTION) — the fixture may be numbered
# differently from the match in the diagram.

# do two nodes fit under the offset `d` (kind, m, degree, slot colours)?
function _zamo_node_fits(fnd::Node, gnd::Node, d::Int)
    fnd.kind === gnd.kind || return false
    fnd.m == gnd.m || return false
    deg = _node_degree(fnd)
    deg == _node_degree(gnd) || return false
    return all(k -> _slot_colour(gnd, mod1(k + d, deg)) == _slot_colour(fnd, k), 1:deg)
end

"""
    _zamo_cluster_match(F, g, fa, ga, d) -> Union{Nothing, NamedTuple}

Structural isomorphism of the fixture `F` (ALL of whose nodes belong to the
pattern) into `g`, anchored at `fa ↦ ga` with slot offset `d`. Propagation along
the inner edges of `F`; every one must recur in `g` with the same slot and colour.

Result `(phi, shift, ext)`: `phi[fnode] = gnode`, `shift[fnode] = d_fnode`, and
`ext[k]` the port of `g` at which leaf `k` of `F` hangs (the outer connection of
the cluster). `nothing` as soon as something does not fit.
"""
function _zamo_cluster_match(F::WordGraph, g::WordGraph, fa::Int, ga::Int, d::Int)
    _zamo_node_fits(F.nodes[fa], g.nodes[ga], d) || return nothing
    phi = Dict{Int,Int}(fa => ga)
    shift = Dict{Int,Int}(fa => d)
    ext = Dict{Int,Port}()
    gslots = Dict{Int,Dict{Int,Port}}()
    getslots(ni) = get!(gslots, ni) do
        s = _slot_map(g, ni)
        s === nothing ? Dict{Int,Port}() : s
    end
    queue = [fa]
    while !isempty(queue)
        a = popfirst!(queue)
        A, da = phi[a], shift[a]
        degA = _node_degree(F.nodes[a])
        fslots = _slot_map(F, a)
        fslots === nothing && return nothing
        gs = getslots(A)
        length(gs) == _node_degree(g.nodes[A]) || return nothing
        for (sa, po) in fslots
            qo = get(gs, mod1(sa + da, degA), nothing)
            qo === nothing && return nothing
            if po isa NodePort
                b, sb = po.node, po.slot
                qo isa NodePort || return nothing
                degB = _node_degree(F.nodes[b])
                db = mod(qo.slot - sb, degB)
                if haskey(phi, b)
                    (phi[b] == qo.node && shift[b] == db) || return nothing
                else
                    qo.node in values(phi) && return nothing      # injective
                    _zamo_node_fits(F.nodes[b], g.nodes[qo.node], db) || return nothing
                    phi[b] = qo.node; shift[b] = db; push!(queue, b)
                end
            elseif po isa Leaf
                ext[po.k] = qo
            else
                return nothing            # circle in the pattern: not provided for
            end
        end
    end
    length(phi) == length(F.nodes) || return nothing
    return (phi = phi, shift = shift, ext = ext)
end

"""
    find_zamo_matches(g::WordGraph; pair = (1, 8), prefilter = true)
        -> Vector{NamedTuple}

All places in `g` where the LHS fixture `flip(Zamo(j,i))` of the pair
`pair = (i, j)` sits as a subdiagram. Each hit is

    (anchor, phi, shift, ext, roles)

with `anchor` = the graph node corresponding to the anchor 23-braid of the
fixture, and `roles` = the result of `zamo_anchor_roles` there.

With `prefilter = true` (default) only anchors satisfying the conditions
(`zamo_anchor_roles`) are tried. `prefilter = false` tries EVERY node — that makes
it measurable whether the prefilter throws away anything the isomorphism would
have accepted.

**TWO PREFILTERS, depending on the fixture.** If the LHS fixture of the pair
satisfies the strict conditions from `zamo_anchor_roles` (as it does for the pair
`(1,8)`), those are used for filtering. If it does NOT, the filter falls back to
the weaker condition for `(2,9)`: *"only look at 23-braids, all 6
neighbours wired as in the fixture"* — i.e. every 23-braid is a candidate, and the
isomorphism checks the wiring. Which filter applies is in `zamo_prefilter_kind`.
"""
function find_zamo_matches(g::WordGraph; pair::Tuple{Int,Int} = (1, 8),
                           prefilter::Bool = true)
    lhs, _ = zamo_rule_fixtures(pair[1], pair[2])
    F = lhs.graph
    strict = zamo_anchors(F)
    fanchors = isempty(strict) ? _zamo_braid23_anchors(F) : strict
    isempty(fanchors) && error(
        "find_zamo_matches: the LHS fixture for $(pair) does not even contain a " *
        "23-braid — the rule is not formulated for this pair")
    cands = !prefilter ? collect(1:length(g.nodes)) :
            isempty(strict) ? _zamo_braid23_anchors(g) : zamo_anchors(g)
    out = NamedTuple[]
    seen = Set{Any}()
    for fa in fanchors, ga in cands, d in 0:(_node_degree(F.nodes[fa]) - 1)
        m = _zamo_cluster_match(F, g, fa, ga, d)
        m === nothing && continue
        key = sort([(k, v) for (k, v) in m.phi])
        key in seen && continue
        push!(seen, key)
        push!(out, (anchor = ga, phi = m.phi, shift = m.shift, ext = m.ext,
                    roles = zamo_anchor_roles(g, ga)))
    end
    return out
end

# ---- the replacement ----------------------------------------------------------

"""
    apply_zamo_rule(g::WordGraph, anchor::Int; pair = (1, 8), check = true)
        -> Union{Nothing, WordGraph}

Applies the Zamolodchikov rule `rev(Zamo(j,i)) → Zamo(i,j)` at the site
`anchor`: the seven nodes of the LHS fixture are deleted and replaced by the seven
nodes of the RHS fixture, at the same twelve outer connections.

`nothing` if there is no match at `anchor`. With `check = true` it is verified
afterwards that the new nodes are fully wired and the edge count is right.

**Do not build this into a driver** — see the file header.
"""
function apply_zamo_rule(g::WordGraph, anchor::Int; pair::Tuple{Int,Int} = (1, 8),
                         check::Bool = true)
    ms = [m for m in find_zamo_matches(g; pair = pair) if m.anchor == anchor]
    isempty(ms) && return nothing
    m = ms[1]
    lhs, rhs = zamo_rule_fixtures(pair[1], pair[2])
    align = _zamo_leaf_alignment(lhs, rhs)
    R = rhs.graph

    dead = Set{Int}(values(m.phi))
    drop = Set{Int}()
    for ni in dead, (ei, _, _) in _edges_at_node(g, ni)
        push!(drop, ei)
    end
    gplus = WordGraph(g.word, vcat(g.nodes, R.nodes), g.edges)
    off = length(g.nodes)
    keep = Edge[e for (j, e) in enumerate(gplus.edges) if !(j in drop)]

    # leaf k of the RHS fixture ↦ the outer connection the CORRESPONDING
    # leaf of the LHS fixture is attached to.
    back = Dict{Int,Int}(v => k for (k, v) in align)
    newport(p::NodePort) = NodePort(off + p.node, p.slot)
    function newport(p::Leaf)
        haskey(back, p.k) || error("apply_zamo_rule: RHS leaf $(p.k) has no partner")
        haskey(m.ext, back[p.k]) || error(
            "apply_zamo_rule: no outer connection for LHS leaf $(back[p.k])")
        return m.ext[back[p.k]]
    end
    newport(p::Circle) = p
    for e in R.edges
        push!(keep, Edge(e.colour, newport(e.a), newport(e.b)))
    end
    out = _delete_nodes(gplus, dead, keep)

    if check
        length(out.nodes) == length(g.nodes) || error(
            "apply_zamo_rule: node count changed ($(length(g.nodes)) → $(length(out.nodes)))")
        length(out.edges) == length(g.edges) || error(
            "apply_zamo_rule: edge count changed ($(length(g.edges)) → $(length(out.edges)))")
        v = check_wiring(out)
        isempty(v) || error("apply_zamo_rule: the result violates the wiring: $v")
    end
    return out
end

"""
    apply_zamo_rule(m::MorphismGraph, anchor::Int; kwargs...) -> Union{Nothing, MorphismGraph}

The same rule on a morphism; the cuts stay where they are (the rule does not
touch the boundary).
"""
function apply_zamo_rule(m::MorphismGraph, anchor::Int; kwargs...)
    out = apply_zamo_rule(m.graph, anchor; kwargs...)
    out === nothing && return nothing
    return MorphismGraph(out, m.cut1, m.cut2)
end
