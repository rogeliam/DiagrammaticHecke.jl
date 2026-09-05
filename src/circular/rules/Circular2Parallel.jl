# circular/rules/Circular2Parallel.jl — 2parallel: two parallel edges of the SAME colour.
#
# THE RULE: two non-touching edges of colour `s` lying on the same region `R`
# merge into one 4-armed `s`-node. `R` splits into the two sectors flanked by
# arms of DIFFERENT original edges; the result is a sum of two terms
# (coefficient 1 each) with `α_s/2` once in each new region — the
# decomposition of the identity on `B_s B_s` into its two idempotents (EMTW).
# Degree checks out on both sides (the 4-node = merge+split costs −2, `α_s/2`
# gives +2). Derived and measured.
#
# WHY THE RULE EXISTS: an inner region of boundary length >= 4 is a sign of
# not being in normal form (DL never produces one), and 2parallel is what
# breaks it into triangles: normal forms carry only length-2/3 inner regions,
# never length 4, while raw double leaves do have such spots.
#
# CONDITIONS:
#  1. same colour `s`, edges do not touch (`_circular_edges_touch`);
#  2. both lie on the SAME region `R` — any region, interior or boundary;
#  3. neither boundary arc between them carries a further `s`-edge
#     (`_circular_parallel_connection_clean`, same guard as P1, colour as parameter);
#  4. the pair is the (tree parent edge, transition edge) of the region carrying
#     the shortest unreduced transition (`circular_unreduced_transition`,
#     CircularRegionWord.jl, case `:direct`, kind `:edge`);
#  5. `R` carries a trivial label (else region-splitting can't decide which
#     half inherits the polynomial; not a real restriction in the driver,
#     which only runs this rule after the polynomials are cleared);
#  6. `strict` (default `true`): neither candidate edge already hangs off an
#     all-`s`-coloured node with >= 3 arms — the termination measure below.
#
# WHY THE REGION WORD DECIDES CONDITION 4. It is the same pair rex fusion fuses
# at `case == :direct`, so the two rules agree on the site instead of each
# picking their own. The rule therefore has nothing left to decide: the
# transition names the region and the pair, and what remains here is the search
# for that pair plus the guards.
#
# TERMINATION: the measure is the number of `s`-edges hanging off no all-`s`
# node with >= 3 arms (`_circular_edge_at_general_2parallel`). Every firing removes
# two such edges and creates none. Condition 4 does not enter the argument — it
# is a FILTER on candidates, never touching the measure or the condition 6 that
# feeds it, so it only changes how fast the measure falls. Runs in `reduce_to_circular_leave`, not `CIRCULAR_RULES`
# (it INCREASES `circular_weight`), in the `:after_d4` slot (`_CIRCULAR_2PARALLEL_SLOT`,
# CircularDecoratedRules.jl) — the last resort before a term counts as done.
# Driver-wide termination is separately guarded by `_circular_leave_growth_guard`.
#
# Selection among several hits: region number, then edge index — deterministic.
"""
    _circular_edge_at_general_2parallel(g::CircularGraph, e::Edge, s::Int) -> Bool

`true` if one end of `e` hangs off a node with ≥ 3 arms whose colours are ALL
`s` — exactly the kind of node [`circular_2parallel_apply`](@ref) creates. Such
edges are not candidates in strict mode; this is the rule's termination
measure (file header).
"""
function _circular_edge_at_general_2parallel(g::CircularGraph, e::Edge, s::Int)
    for p in (e.a, e.b)
        p isa NodePort || continue
        nd = g.nodes[p.node]
        arm_count(nd) >= 3 && all(==(s), nd.arms) && return true
    end
    return false
end

"""
    find_circular_2parallel(m::CircularMorphismGraph; strict = true)
        -> Union{Nothing, NamedTuple}

Searches for the candidate pair for **2parallel**: two non-touching edges
**of the same colour** on **the same region** — any region, interior or
boundary —, with a clean (= differently-coloured) connection on both
boundary arcs.

**The deciding condition** is the region word: the pair must be the (parent,
transition) edge pair of the region carrying the shortest unreduced transition
(`circular_unreduced_transition`, case `:direct`) — the same pair rex fusion
fuses there.

A marking must be set (`cut1` from `m`); without one there is no transition and
the rule never fires.

Returns `(ei, ej, region, colour, far_lo, far_hi)` or `nothing`. `far_lo`/
`far_hi` are the two neighbour regions beyond the two edges; they're just a
fixed order now (by distance when a marking is set, else by region number).
No caller in the driver reads them; they exist for notebooks and the log.

Deterministic: regions in order `1:region_count`, within a region the smaller
edge indices first.
"""
function find_circular_2parallel(m::CircularMorphismGraph; strict::Bool = true)
    g = m.graph
    dist = circular_region_distances(m)        # INFORMATIONAL only (log/return)
    adj, _ = circular_region_adjacency(g)

    # For the boundary-arc guard: dart -> edge index (all colours, since the
    # differently-coloured pieces on the arc must be recognized).
    _, t = _circular_dart_region(g)
    dart_edge = Dict{Int, Int}()
    for (ei, e) in enumerate(g.edges)
        (haskey(t.port_dart, e.a) && haskey(t.port_dart, e.b)) || continue
        dart_edge[t.port_dart[e.a]] = ei
        dart_edge[t.port_dart[e.b]] = ei
    end

    # ORDER: simply `1:region_count`, within a region the smaller edge indices.
    # Only one region can match anyway — the one the transition names — so the
    # order only fixes which pair is reported inside it.
    # The shortest unreduced transition, once per call (CircularRegionWord.jl).
    # It names the region and the pair; without one there is nothing to fuse.
    transition = circular_unreduced_transition(m)
    (transition !== nothing &&
     transition.case === :direct && transition.kind === :edge &&
     transition.prev_edge != 0) || return nothing
    pair = minmax(transition.prev_edge, transition.edge)
    for R in 1:region_count(g)
        lst = sort(adj[R]; by = x -> x[2])
        for a in 1:length(lst), b in (a + 1):length(lst)
            (Ra, ei) = lst[a]
            (Rb, ej) = lst[b]
            s = g.edges[ei].colour
            s == g.edges[ej].colour || continue
            _circular_edges_touch(g.edges[ei], g.edges[ej]) && continue
            # Condition 4: the pair is EXACTLY (tree parent edge, transition
            # edge) of the region the transition names. In case `:direct` the
            # region word ends on `s`, so the parent tree edge IS s-coloured and
            # the pair is well formed; the other cases are handled by the
            # rex-cycle derivation, not here. No distance test is needed.
            (R == transition.region && (ei, ej) == pair) || continue

            if strict
                (_circular_edge_at_general_2parallel(g, g.edges[ei], s) ||
                 _circular_edge_at_general_2parallel(g, g.edges[ej], s)) && continue
            end
            _circular_parallel_connection_clean(g, t, dart_edge, ei, ej, R;
                                           forbidden = Set([s])) || continue
            # `far_lo`/`far_hi`: by distance when a marking is set, else by
            # region number — just a fixed order, not a condition.
            key(X) = (dist[X] >= 0 ? dist[X] : typemax(Int), X)
            lo, hi = key(Ra) <= key(Rb) ? (Ra, Rb) : (Rb, Ra)

            return (ei = ei, ej = ej, region = R, colour = s,
                    far_lo = lo, far_hi = hi)
        end
    end
    return nothing
end
