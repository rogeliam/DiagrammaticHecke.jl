# circular/rules/CircularDotSlide.jl — the "dot-slide" rule as a driver step.
#
# THE IDEA: this equation is part of the D4 probe; there, start with the dots
# at the greatest distance first.
# Worked out for the diagram "two stacked 121-braids + 2-dot": it reduces
# directly to
# `id₁₂₁ ⊗ 2-dot + rest`; if D4 is forced FIRST at the dot (the middle
# 2-edge between the braids), the SAME diagram instead splits into terms
# WITHOUT the identity summand. The difference:
#
#     id-window (3 strands + dot)
#         =  [braid²-window, D4 forced at the dot immediately → 3 terms]
#         −  [rest-window]
#
# LOCAL WINDOW. A dot of color t whose distance path to the mark crosses, as
# its first three edges, the colors (a, t, a), where {a,t} is a braid pair
# (|a-t| == 1 — the commuting pair 1/3 is EXCLUDED). The three crossed edges
# are three nested strands around the dot.
#
# * Braid²-window: the three strands are replaced by TWO 6-valent braid
#   nodes `circular_node([a,t,a,t,a,t])` with three middle edges (colors t, a, t);
#   the dot stays. The forced D4 fires on one of the two t-middle edges.
#   Wiring template (color map 1→a, 2→t): b1 slots 1/6/5 to one side of the
#   strands (a,t,a), b2 slots 1/6/5 to the other side, middle edges
#   (b1,4)-(b2,2) color t, (b1,3)-(b2,3) color a, (b1,2)-(b2,4) color t.
#   WHICH strand end is "side 1"/"side 2" and which slot assignment is
#   planar is NOT guessed: every assignment is tried and the first one taken
#   that is planar AND on which the forced D4 fires at the dot (house
#   style — `circular_2parallel_apply` and the preimage stages search the same
#   way).
# * Rest-window: the dot and the three strands are replaced by ONE 4-valent
#   a-node `circular_node([a,a,a,a])` (takes the four a-strand ends) plus THREE
#   t-dots `circular_node([t])` (cap the two t-strand ends and the dot's old cap
#   end). Slot order again by planarity search.
#
# ┌────────────────────────────────────────────────────────────────────────────┐
# │ TERMINATION: NOT PROVEN. The step replaces one term by FOUR (3 D4 terms + │
# │ rest) and has no falling measure, so the driver does not call it; the     │
# │ `:dot` transitions go through `circular_dot_fusion_step` instead.         │
# └────────────────────────────────────────────────────────────────────────────┘
#
# THIS IS LABEL-FREE ONLY: `circular_dot_slide_step` requires
# `all(isone, region_labels) && isone(outer_label)`. Transferring labels
# through the double surgery (three edges out, two braids in, resp. a node
# out, four nodes in) is left out; it could go through
# `_circular_transfer_labels`.

"""
    CIRCULAR_DOT_SLIDE_ENABLED

Default **`false`**. `reduce_to_circular_leave` does not call this step: the
`:dot` transitions are handled by `circular_dot_fusion_step`, whose derived
relation gives the same result at the `(a, t, a)` window. The hand-built
surgery stays as the body of `circular_regionword_dihedral_step` and can be
called independently of the switch.

The slide is the BASE CASE. The general statement is the region word: a dot
sitting in a region whose word is `w`, where `ws` is unreduced, is the shape
this surgery handles at its shortest — and that is the same criterion
`circular_unreduced_transition` names, which is why the region-word step is
where the surgery lives.
"""
const CIRCULAR_DOT_SLIDE_ENABLED = Ref(false)

# ---- distance paths -------------------------------------------------------

"""
    _fds_paths(adj, dist, R) -> Vector{NTuple{3,Int}}

All edge triples `(e1, e2, e3)` for the first three steps of a distance path
down from region `R`: at each step, a neighbor at distance `d - 1` via
`circular_region_adjacency` edges. If a stage has SEVERAL `d - 1` neighbors, all
combinations are returned (the graphs are small).
"""
function _fds_paths(adj, dist, R::Int)
    out = NTuple{3,Int}[]
    dist[R] >= 3 || return out
    for (R1, e1) in adj[R]
        dist[R1] == dist[R] - 1 || continue
        for (R2, e2) in adj[R1]
            dist[R2] == dist[R1] - 1 || continue
            for (R3, e3) in adj[R2]
                dist[R3] == dist[R2] - 1 || continue
                push!(out, (e1, e2, e3))
            end
        end
    end
    return out
end

"The two ports of an edge."
_fds_ends(e::Edge) = (e.a, e.b)

"Does edge `e` touch node `ni`?"
_fds_touches(e::Edge, ni::Int) =
    (e.a isa NodePort && e.a.node == ni) || (e.b isa NodePort && e.b.node == ni)

# ---- surgery T': the braid²-window -----------------------------------------

"""
    _fds_braid2(g, dot, cap_ei, e1, e2, e3, a, t)
        -> Union{Nothing, Tuple{CircularGraph, Int, Int}}

Builds the braid²-window: removes the three strand edges `e1` (a, at the
dot), `e2` (t), `e3` (a, at the mark), attaches two 6-valent braid nodes
`circular_node([a,t,a,t,a,t])`, and connects the six free strand ends per the
derived template. Which side choice/slot assignment is planar is found by
searching ALL assignments (file header); the first one accepted has
`euler == 2`, empty `check_wiring`, `is_wired`, AND a firing `apply_circular_d4`
at the dot on a t-middle edge. Returns `(new graph, dot index, index of the
D4 middle edge)` or `nothing`.
"""
function _fds_braid2(g::CircularGraph, dot::Int, cap_ei::Int,
                     e1::Int, e2::Int, e3::Int, a::Int, t::Int)
    b1 = length(g.nodes) + 1
    b2 = length(g.nodes) + 2
    nodes2 = vcat(g.nodes, [circular_node([a, t, a, t, a, t]),
                            circular_node([a, t, a, t, a, t])])
    removed = (e1, e2, e3)
    kept = Edge[ee for (j, ee) in enumerate(g.edges) if !(j in removed)]
    strand_ends = [_fds_ends(g.edges[e1]), _fds_ends(g.edges[e2]),
                   _fds_ends(g.edges[e3])]

    # Slot triples (a-slot A, t-slot, a-slot B) — both mirror forms of the
    # template; on `[a,t,a,t,a,t]`, slots 1/3/5 carry color a, 2/4/6 color t.
    triples = ((1, 6, 5), (1, 2, 3))
    for tr1 in triples, tr2 in triples,
        aswap1 in (false, true), aswap2 in (false, true),
        sides in 0:7, tmswap in (false, true)

        # Side per strand: bit k of `sides` says which end of strand k goes
        # to b1 (the other to b2).
        side1 = [strand_ends[k][((sides >> (k - 1)) & 1) + 1] for k in 1:3]
        side2 = [strand_ends[k][2 - ((sides >> (k - 1)) & 1)] for k in 1:3]
        # a-strand assignment per braid: (strand for slot-triple position 1, 3).
        aA1, aB1 = aswap1 ? (1, 3) : (3, 1)
        aA2, aB2 = aswap2 ? (1, 3) : (3, 1)

        # Attachment slot per strand (1=inner a, 2=t, 3=outer a) at b1/b2, and
        # the middle-edge assignment b1-slot -> b2-slot:
        attach1 = Dict(aA1 => tr1[1], 2 => tr1[2], aB1 => tr1[3])
        attach2 = Dict(aA2 => tr2[1], 2 => tr2[2], aB2 => tr2[3])
        m1 = sort(setdiff(1:6, collect(tr1)))     # free b1 slots: [t, a, t]
        m2 = sort(setdiff(1:6, collect(tr2)))
        t1s = [s for s in m1 if iseven(s)]
        a1  = only(s for s in m1 if isodd(s))
        t2s = [s for s in m2 if iseven(s)]
        a2  = only(s for s in m2 if isodd(s))
        tm  = tmswap ? reverse(t2s) : t2s
        mid = Dict(t1s[1] => tm[1], a1 => a2, t1s[2] => tm[2])

        # IDENTITY CRITERION (the property of the derived template, checked
        # there: 1→4→2→5, 6→3→3→6, 5→2→4→1): every strand must arrive back
        # at its OWN other end after passing through b1, a middle edge, and
        # b2 — `opposite(attach1) --middle edge--> s2, opposite(s2) ==
        # attach2`. Only then is the window really braid² = id (the braid
        # relation applied twice); planar but differently wired candidates
        # represent a DIFFERENT morphism — without this criterion, the search
        # can accept a window whose reduction does not match the direct
        # reduction of the id-window (measured).
        opp(s) = mod1(s + 3, 6)
        all(opp(mid[opp(attach1[k])]) == attach2[k] for k in 1:3) || continue

        edges2 = copy(kept)
        push!(edges2, Edge(a, side1[aA1], NodePort(b1, tr1[1])))
        push!(edges2, Edge(t, side1[2],   NodePort(b1, tr1[2])))
        push!(edges2, Edge(a, side1[aB1], NodePort(b1, tr1[3])))
        push!(edges2, Edge(a, side2[aA2], NodePort(b2, tr2[1])))
        push!(edges2, Edge(t, side2[2],   NodePort(b2, tr2[2])))
        push!(edges2, Edge(a, side2[aB2], NodePort(b2, tr2[3])))
        # Middle edges: the three free slots each (colors t, a, t). The
        # a-middle is color-forced; the two t-middles have two pairings
        # (`tmswap`, already fixed above in `mid`).
        tmid1 = length(edges2) + 1
        push!(edges2, Edge(t, NodePort(b1, t1s[1]), NodePort(b2, tm[1])))
        push!(edges2, Edge(a, NodePort(b1, a1),     NodePort(b2, a2)))
        tmid2 = length(edges2) + 1
        push!(edges2, Edge(t, NodePort(b1, t1s[2]), NodePort(b2, tm[2])))

        gneu = try
            cand = CircularGraph(g.word, nodes2, edges2)
            (euler(cand) == 2 && isempty(check_wiring(cand)) && is_wired(cand)) ||
                continue
            cand
        catch
            continue                              # non-planar ⇒ next assignment
        end
        # The forced D4 must fire on a t-middle edge at the dot — ONLY a
        # firing assignment is accepted (exactly the middle t-edge between
        # the two braids is the move).
        for mei in (tmid1, tmid2)
            hit = try
                apply_circular_d4(circular_decorated(gneu), dot, mei)
            catch
                nothing
            end
            hit === nothing || return (gneu, dot, mei)
        end
    end
    return nothing
end

# ---- surgery T'': the rest-window -------------------------------------------

"""
    _fds_rest(g, dot, cap_ei, dse, e1, e2, e3, a, t) -> Union{Nothing, CircularGraph}

Builds the rest-window: removes `e1, e2, e3`, the cap edge, and the dot node;
attaches ONE 4-valent a-node to the four a-strand ends (slot order by
planarity search over all cyclic orders) and THREE t-dots to the two
t-strand ends and the existing cap end `dse`. Node removal runs through
`_circular_delete_nodes` (clean renumbering). `nothing` if no assignment is
planar.
"""
function _fds_rest(g::CircularGraph, dot::Int, cap_ei::Int, dse::Port,
                   e1::Int, e2::Int, e3::Int, a::Int, t::Int)
    a_ends = Port[_fds_ends(g.edges[e1])..., _fds_ends(g.edges[e3])...]
    t_ends = Port[_fds_ends(g.edges[e2])..., dse]
    X = length(g.nodes) + 1                       # the 4-valent a-node
    nodes2 = vcat(g.nodes, [circular_node([a, a, a, a]),
                            circular_node([t]), circular_node([t]), circular_node([t])])
    removed = (e1, e2, e3, cap_ei)
    kept = Edge[ee for (j, ee) in enumerate(g.edges) if !(j in removed)]

    # All cyclic orders of the four a-ends: first end fixed at slot 1, the
    # remaining three permuted (3! = 6 — this also covers the mirror images).
    for perm in ((2, 3, 4), (2, 4, 3), (3, 2, 4), (3, 4, 2), (4, 2, 3), (4, 3, 2))
        edges2 = copy(kept)
        push!(edges2, Edge(a, a_ends[1], NodePort(X, 1)))
        for (slot, k) in enumerate(perm)
            push!(edges2, Edge(a, a_ends[k], NodePort(X, slot + 1)))
        end
        for k in 1:3
            push!(edges2, Edge(t, t_ends[k], NodePort(X + k, 1)))
        end
        g2 = try
            # as in `apply_circular_d4`: attach the nodes first (edges still the
            # existing ones — the tracer tolerates unwired slots), then delete the
            # dot and let the NEW edge list (old node indices) be rewritten.
            _circular_delete_nodes(CircularGraph(g.word, nodes2, g.edges), Set([dot]), edges2)
        catch
            continue
        end
        ok = try
            euler(g2) == 2 && isempty(check_wiring(g2)) && is_wired(g2)
        catch
            false
        end
        ok && return g2
    end
    return nothing
end

# ---- the step ---------------------------------------------------------------

"""
    circular_dot_slide_step(fdm::CircularDecoratedMorphism) -> Union{Nothing, CircularComboR}

The dot-slide rule (file header): looks for a dot of color `t` whose distance
path to the mark crosses, as its first three edges, the colors `(a, t, a)`
with `|a - t| == 1` (dots at the GREATEST distance first), and replaces the
id-window with the 4-term sum

    [D4 at the dot on the braid²-window: 3 terms]  −  [rest-window].

THIS IS LABEL-FREE ONLY: requires `all(isone, region_labels)` and
`isone(outer_label)` — transferring labels through the double surgery is
left out (file header). `nothing` if no dot has the pattern, or
no planar wiring is found. Callable independently of
`CIRCULAR_DOT_SLIDE_ENABLED` (the switch only controls the driver hook).
"""
function circular_dot_slide_step(fdm::CircularDecoratedMorphism)
    all(isone, fdm.region_labels) && isone(fdm.outer_label) || return nothing
    g = fdm.m.graph
    dist = circular_region_distances(fdm.m)
    all(==(-1), dist) && return nothing
    adj, _ = circular_region_adjacency(g)

    # Dots by DESCENDING distance of their region ("start
    # with the dots at the greatest distance first"); skip distance -1.
    dots = Tuple{Int,Int}[]                       # (node, region)
    for (v, nd) in enumerate(g.nodes)
        arm_count(nd) == 1 || continue
        R = circular_region_of_dot(g, v)
        dist[R] >= 3 || continue                  # three downward steps needed
        push!(dots, (v, R))
    end
    sort!(dots; by = p -> (-dist[p[2]], p[1]))

    for (v, R) in dots
        t = arm_colour(g.nodes[v], 1)
        at = _circular_edges_at_node(g, v)
        length(at) == 1 || continue
        (cap_ei, _, dse) = at[1]
        for (e1, e2, e3) in _fds_paths(adj, dist, R)
            # pattern (a, t, a), braid pair (1/3 commutes and is excluded):
            a = g.edges[e1].colour
            g.edges[e3].colour == a || continue
            g.edges[e2].colour == t || continue
            abs(a - t) == 1 || continue
            # three DIFFERENT edges, none at the dot, no two sharing a port
            # (genuinely nested strands):
            length(unique((e1, e2, e3))) == 3 || continue
            any(_fds_touches(g.edges[e], v) for e in (e1, e2, e3)) && continue
            ports = Port[_fds_ends(g.edges[e1])..., _fds_ends(g.edges[e2])...,
                         _fds_ends(g.edges[e3])...]
            length(unique(ports)) == 6 || continue

            b2res = _fds_braid2(g, v, cap_ei, e1, e2, e3, a, t)
            b2res === nothing && continue
            (gneu, dotidx, mid) = b2res
            rest = _fds_rest(g, v, cap_ei, dse, e1, e2, e3, a, t)
            rest === nothing && continue

            d4 = apply_circular_d4(circular_decorated(gneu), dotidx, mid)
            d4 === nothing && continue            # can't happen (checked above)
            return d4 + (-1) * CircularComboR(circular_decorated(rest))
        end
    end
    return nothing
end
