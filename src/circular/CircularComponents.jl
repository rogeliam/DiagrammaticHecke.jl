# circular/CircularComponents.jl — detect and locally evaluate freely floating connected
# components.
#
# WHAT THIS IS ABOUT. `compose`/`join_graphs` (diagram/Join.jl) only checks
# whether a pure EDGE CHAIN closes into a borderless circle (the existing
# `Circle` marker that `:free_circle` fires on). A chain WITH NODES — e.g. a
# braid node with several dots, all connected only among themselves, touching
# no leaf and no other node — is NOT caught by that and runs unchanged into
# the circular reduction. The only existing place that handles anything like this
# (in miniature: exactly 2 dots) is `_fr_barbell` (circular/rules/CircularRules.jl) —
# its "region without a boundary gap" fallback is the template for step D
# below, but it only ever sees the ONE step it makes itself, not the whole
# component at once.
#
# When diagrams are joined, i.e. top1 with bottom2, a connected component may
# result that floats free. If so, it is extracted and evaluated locally to a
# polynomial; during the join, this component is dropped and the polynomial is
# written directly into the right region.
#
# The check lives here on the CIRCULAR side, as the first step of
# `reduce_to_circular_leave` (CircularDecoratedRules.jl) — that covers not just
# the original `compose` but every recursive re-entry point (fusion, D4, Zamo
# all build new `CircularDecoratedMorphism`s themselves and call
# `reduce_to_circular_leave` again).
#
# Needs `_circular_delete_nodes` (circular/rules/CircularRules.jl) and `regions`/
# `_circular_dart_region` (diagram/Faces.jl, circular/CircularRegion.jl) — the include
# point in DiagrammaticHecke.jl is accordingly AFTER CircularRules.jl.

"""
    circular_connected_components(g::CircularGraph) -> Vector{Vector{Int}}

The connected components of the NODES of `g` (adjacency only via node-node
edges, `NodePort`↔`NodePort`; leaf edges do not count as a connection between
nodes, but they are the test for whether a component touches the boundary —
see [`_circular_component_touches_boundary`](@ref)). Each component as a sorted
list of node indices, component order = order of the first node.
"""
function circular_connected_components(g::CircularGraph)
    n = length(g.nodes)
    adj = [Int[] for _ in 1:n]
    for e in g.edges
        e.a isa NodePort && e.b isa NodePort || continue
        push!(adj[e.a.node], e.b.node)
        push!(adj[e.b.node], e.a.node)
    end
    seen = falses(n)
    comps = Vector{Int}[]
    for start in 1:n
        seen[start] && continue
        comp = Int[]
        stack = [start]
        while !isempty(stack)
            v = pop!(stack)
            seen[v] && continue
            seen[v] = true
            push!(comp, v)
            append!(stack, adj[v])
        end
        push!(comps, sort!(comp))
    end
    return comps
end

"`true` if any node of `comp` hangs off the boundary of `g` via a leaf edge."
function _circular_component_touches_boundary(g::CircularGraph, comp::Vector{Int})
    compset = Set(comp)
    return any(e -> (e.a isa NodePort && e.a.node in compset && e.b isa Leaf) ||
                    (e.b isa NodePort && e.b.node in compset && e.a isa Leaf),
               g.edges)
end

"""
    _circular_evaluate_closed_component(g::CircularGraph, comp::Vector{Int}) -> SoergelPoly

Evaluates the connected component `comp` (node indices in `g`, assumed to hang
off neither a leaf nor a node OUTSIDE `comp`) as a standalone, borderless
diagram: builds an isolated `CircularGraph` from just these nodes/edges (empty
boundary word), reduces it recursively with [`reduce_to_circular_leave`](@ref), and
reads off the result as a pure scalar.

A borderless connected diagram reduces, by the Soergel calculus, to a single
polynomial — either as the coefficient of an empty sum (the reduction already
extracts all scalars itself) or as the one label of the one remaining region.
Throws if the result is something ELSE (several terms, a term with actual
nodes remaining) — that would be its own finding, not a case for a silent
fallback.
"""
function _circular_evaluate_closed_component(g::CircularGraph, comp::Vector{Int})
    remap = Dict{Int,Int}(v => i for (i, v) in enumerate(comp))
    compset = Set(comp)
    newnodes = [g.nodes[v] for v in comp]
    newedges = Edge[]
    for e in g.edges
        e.a isa NodePort && e.a.node in compset || continue
        e.b isa NodePort && e.b.node in compset || continue
        push!(newedges, Edge(e.colour, NodePort(remap[e.a.node], e.a.slot),
                                        NodePort(remap[e.b.node], e.b.slot)))
    end
    sub = CircularGraph(CircularWord(Int[]), newnodes, newedges)
    fdm = CircularDecoratedMorphism(CircularMorphismGraph(sub, 0, 0),
                                fill(one(SoergelPoly), region_count(sub)))
    combo = reduce_to_circular_leave(fdm)

    terms = collect(pairs_of(combo))
    isempty(terms) && return zero(SoergelPoly)

    # SUM OVER ALL TERMS — do not insist on a single one. A borderless component
    # reduces to empty diagrams, and the scalar it is worth sits in the LABELS,
    # not in the coefficient: `_fr_barbell` puts its
    # `α` into `outer_label`, and `circular_extract_scalars` cannot pull it out here
    # because a borderless sub-diagram has no marking and hence no distances.
    # Two empty diagrams whose outer labels differ are therefore DIFFERENT keys
    # and stay apart in the combo — a component can legitimately be worth a
    # sum such as `α₁²α₂ + α₁α₂²`, coming back as several terms. Requiring a
    # single term would turn that sum into an error; summing them is the value
    # the caller wants.
    total = zero(SoergelPoly)
    for (dd, coeff) in terms
        isempty(dd.graph.nodes) || error(
            "_circular_evaluate_closed_component: borderless component does NOT reduce " *
            "to an empty diagram ($(length(dd.graph.nodes)) nodes remain) — " *
            "not a pure scalar.")
        # A completely empty diagram (0 nodes, 0 edges) has 0 regions; a single
        # leftover borderless cell (circle/lens/barbell whose label the fusion
        # has not pulled out) has 1. Anything else is unexpected.
        rc = region_count(dd.graph)
        rc in (0, 1) || error(
            "_circular_evaluate_closed_component: borderless remainder diagram has " *
            "$rc regions instead of 0 or 1 — unexpected.")
        # THE OUTER LABEL BELONGS TO THE SCALAR: the outer
        # region of this sub-diagram is exactly the host region the polynomial is
        # about to be written into — it is part of the value, not a second place.
        total += rc == 0 ? coeff * dd.outer_label :
                           coeff * dd.region_labels[1] * dd.outer_label
    end
    return total
end

"""
    circular_extract_floating_components(fdm::CircularDecoratedMorphism)
        -> (CircularDecoratedMorphism, SoergelPoly)

Detects connected components of `fdm.m.graph` that hang off NEITHER a leaf NOR
another component ("floating free" — can arise under `compose` when a piece of
diagram, after gluing, is only connected to itself), evaluates each one to a
polynomial ([`_circular_evaluate_closed_component`](@ref)), and writes the
polynomial into the region the component GEOMETRICALLY sat in (via the dart
tracer, generalising `circular_region_of_dot`'s mechanism — not via the incremental
rule-by-rule bookkeeping of the individual `CIRCULAR_RULES`).

Returns: the cleaned-up `CircularDecoratedMorphism` (components removed, their
polynomials written into the affected regions, or — where the region was
borderless — into the OUTER region) and a scalar factor (to be multiplied into
`fdm`'s overall value — analogous to `circular_extract_scalars`'s `factor`).
`(fdm, one(SoergelPoly))` if no component floats free (the normal case — then
this step only costs the connectivity check).

**The factor is always `1`**: the polynomial
of a borderless component region stays in place in the outer region, instead
of becoming the coefficient. It only gets pulled out by
`circular_extract_scalars`, and only at **distance 0**: a scalar may only be
pulled out of the diagram when it sits in the distance-0 region.

**The outer region is always at distance 0** whenever
there is a marking — so whatever lands here ends up as a coefficient
immediately, via `circular_extract_scalars`, in the marked driver. The path
(park it outside first, then pull it) is the same either way.
"""
function circular_extract_floating_components(fdm::CircularDecoratedMorphism)
    g = fdm.m.graph
    # A borderless word (empty `word`, no leaves) by definition has NO boundary
    # a component could touch — "floating free" is then no longer a meaningful
    # question, but ALWAYS true, for EVERY component. Without this guard,
    # `_circular_evaluate_closed_component` (which builds exactly such a borderless
    # graph) would call itself again via `reduce_to_circular_leave` and try to
    # "detect" its own, already isolated component again — infinite recursion.
    isempty(g.word) && return fdm, one(SoergelPoly)

    # `factor` is ALWAYS `1`: the polynomial of a borderless component region
    # stays in place in the outer region (`outer_cur`) instead of becoming the
    # scalar. The return value stays around anyway — it is the interface to
    # `reduce_to_circular_leave` and the place a future genuine scalar would arrive.
    factor = one(SoergelPoly)
    g_cur = g
    labels_cur = fdm.region_labels
    outer_cur = fdm.outer_label
    changed = false
    # Recompute components EVERY round: `_circular_delete_nodes` renumbers the
    # nodes, so a list of "all floating components" from BEFORE the first
    # deletion would be indexed wrongly for the second deletion and beyond.
    while true
        comps = circular_connected_components(g_cur)
        floating = [c for c in comps if !_circular_component_touches_boundary(g_cur, c)]
        isempty(floating) && break
        changed = true
        comp = floating[1]
        poly = _circular_evaluate_closed_component(g_cur, comp)
        dart_region, t = _circular_dart_region(g_cur)
        v0 = comp[1]
        slot0 = findfirst(s -> haskey(t.slot_dart, (v0, s)), 1:arm_count(g_cur.nodes[v0]))
        slot0 === nothing && error(
            "circular_extract_floating_components: component node $v0 has no wired slot")
        R = dart_region[t.slot_dart[(v0, slot0)]]

        dead = Set(comp)
        keep_edges = [e for e in g_cur.edges if
                      !((e.a isa NodePort && e.a.node in dead) ||
                        (e.b isa NodePort && e.b.node in dead))]
        g_new = _circular_delete_nodes(g_cur, dead, keep_edges)

        ro = regions(g_cur)
        # Labels on regions that die WITH the component (the phantom region
        # around it, the interior of a free circle) have nowhere left to go
        # after the deletion. They go into the OUTER region — the same
        # treatment `poly` already gets in the borderless branch below.
        outer = Ref(outer_cur)
        if isempty(ro[R].gaps)
            outer[] = outer[] * poly
            labels_cur = _circular_transfer_labels(g_cur, g_new, labels_cur; outer = outer)
        else
            rn = regions(g_new)
            kg = findfirst(S -> ro[R].gaps[1] in S.gaps, rn)
            kg === nothing && error(
                "circular_extract_floating_components: boundary gap $(ro[R].gaps[1]) of " *
                "the component region is in no region of the cleaned-up graph.")
            labels_cur = _circular_transfer_labels(g_cur, g_new, labels_cur;
                overrides = Dict(kg => labels_cur[R] * poly), outer = outer)
        end
        outer_cur = outer[]
        g_cur = g_new
    end

    changed || return fdm, one(SoergelPoly)
    m_new = CircularMorphismGraph(g_cur, fdm.m.cut1, fdm.m.cut2)
    return CircularDecoratedMorphism(m_new, labels_cur, outer_cur), factor
end
