# test/circularbraidchannels.jl — the CHANNEL normalisation at braid pairs
# (src/circular/rules/CircularBraidChannels.jl).
#
# The fixtures below are built with the planar constructor and frozen here as
# literals — each is wiring-clean with `euler == 2`.
#
# Everything is checked against the two invariants EVERY rule must preserve: the
# boundary word and `circular_degree`. Plus the termination evidence (`circular_weight` falls) —
# without which the rule could not stand in `CIRCULAR_RULES`.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_braid_back, _fr_braid_relation,
                     _fr_braid_back_gen, _fr_braid_relation_gen,
                     _circular_braid_channels, _circular_pull_trivalent, _circular_dissolve_through,
                     _circular_first_rule_match

# A1: two braids, THREE direct channels (classic C7).
_bc_classic_back() = CircularGraph(CircularWord([2, 2, 3, 2, 2, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(1, 2)),
    ])

# B1: two direct channels + a pure trivalent node (classic C8).
_bc_classic_rel() = CircularGraph(CircularWord([2, 2, 3, 2, 3, 2, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]), circular_node([3, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 2), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(3, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(2, 4)),
    ])

# A2: third channel through a pass-through [1,1,3,3] (nc = 2).
_bc_passthrough() = CircularGraph(CircularWord([1, 1, 2, 3, 2, 2, 3, 2]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]), circular_node([1, 1, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(3, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 4), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(1, 1)),
    ])

# B3a: general-13 node [1,1,3,3,3], connection slots (3,4) — the arc wraps around.
_bc_gen13_a() = CircularGraph(CircularWord([1, 1, 3, 2, 3, 2, 2, 3, 2]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]), circular_node([1, 1, 3, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(3, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 4), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(3, 5)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(1, 1)),
    ])

# B3b: general-13 node [1,1,3,3,3], connection slots (4,5).
_bc_gen13_b() = CircularGraph(CircularWord([1, 1, 2, 3, 2, 2, 3, 2, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]), circular_node([1, 1, 3, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(3, 4)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 5), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(3, 3)),
    ])

# B3c: general-13 node [1,3,1,3,3] — the 1-line CROSSES a leg.
_bc_gen13_crossing() = CircularGraph(CircularWord([1, 2, 3, 2, 2, 3, 2, 1, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]), circular_node([1, 3, 1, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(3, 4)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 5), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(3, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(3, 2)),
    ])

# Both outer channels occupied: pass-through [1,1,3,3] + general-13 [1,3,1,3,3].
_bc_two_sided() = CircularGraph(CircularWord([1, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]), circular_node([1, 1, 3, 3]), circular_node([1, 3, 1, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(3, 4)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 3), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(4, 4)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(4, 5), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(4, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(4, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(4, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(10), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(11), DiagrammaticHecke.NodePort(2, 3)),
    ])


# A3: a third channel across a pass-through with THREE foreign arms [1,1,1,3,3]
# (nc = 2, d = 5). Line by line B3b, except that its third 3-arm is here a
# 1-arm — as a result, dissolving leaves a genuine node [1,1,1] standing instead
# of an edge.
_bc_passthrough_three() = CircularGraph(CircularWord([1, 1, 2, 3, 2, 2, 3, 2, 1]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]), circular_node([1, 1, 1, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(3, 4)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 5), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(3, 3)),
    ])

# The A3 diagram automorphism 1 <-> 3: a 23-braid becomes a 12-braid.
# This tests the same fixtures for the 12-case: the rule never names concrete
# colours anywhere, so this must hold too.
_bc_sw(c) = c == 1 ? 3 : c == 3 ? 1 : 2
_bc_swap13(g::CircularGraph) = CircularGraph(
    CircularWord([_bc_sw(c) for c in letters(g.word)]),
    CircularNode[circular_node([_bc_sw(c) for c in arms(nd)]) for nd in g.nodes],
    DiagrammaticHecke.Edge[DiagrammaticHecke.Edge(_bc_sw(e.colour), e.a, e.b) for e in g.edges])

# The rule output of a diagram: (rule name, terms as CircularGraph).
function _bc_fires(g::CircularGraph)
    h = _circular_first_rule_match(circular_decorated(g), CIRCULAR_RULES)
    h === nothing && return (nothing, CircularGraph[])
    terms = CircularGraph[(t isa CircularGraph ? t : t.graph) for (t, _) in pairs_of(h[2])]
    return (h[1], terms)
end

# Every term must preserve the boundary word and the degree, be planar and
# wiring-clean, and `circular_weight` must fall strictly (the driver's assertion).
function _bc_pruefe(g::CircularGraph, terms::Vector{CircularGraph})
    for t in terms
        @test letters(t.word) == letters(g.word)
        @test circular_degree(t) == circular_degree(g)
        @test euler(t) == 2
        @test isempty(check_wiring(t))
        @test circular_weight(t) < circular_weight(g)
    end
end

@testset "channel normalisation at braid pairs (C7/C8 with general-13 nodes)" begin

    @testset "the fixtures are planar and wiring-clean" begin
        for g in (_bc_classic_back(), _bc_classic_rel(), _bc_passthrough(),
                  _bc_gen13_a(), _bc_gen13_b(), _bc_gen13_crossing(), _bc_two_sided())
            @test euler(g) == 2
            @test isempty(check_wiring(g))
        end
    end

    @testset "the classic case stays UNCHANGED" begin
        # The channel version calls the base rule first — with three direct channels
        # resp. a pure trivalent nothing may change.
        gA = _bc_classic_back()
        @test _fr_braid_back(gA) !== nothing
        @test isempty(_circular_braid_channels(gA, 1, 2))            # no channel node at all
        @test circular_canonical_key.(sort([(t isa CircularGraph ? t : t.graph)
                                       for (t, _) in pairs_of(_fr_braid_back_gen(gA))],
                                      by = x -> length((x isa CircularGraph ? x : x.graph).nodes))) ==
              circular_canonical_key.(sort([(t isa CircularGraph ? t : t.graph)
                                       for (t, _) in pairs_of(_fr_braid_back(gA))],
                                      by = x -> length((x isa CircularGraph ? x : x.graph).nodes)))

        gB = _bc_classic_rel()
        @test _fr_braid_relation(gB) !== nothing
        @test [circular_canonical_key(t isa CircularGraph ? t : t.graph)
               for (t, _) in pairs_of(_fr_braid_relation_gen(gB))] ==
              [circular_canonical_key(t isa CircularGraph ? t : t.graph)
               for (t, _) in pairs_of(_fr_braid_relation(gB))]
    end

    @testset "nc = 2 (pass-through) ⇒ dissolve, then C7" begin
        g = _bc_passthrough()
        chs = _circular_braid_channels(g, 1, 2)
        @test length(chs) == 1
        @test chs[1].nc == 2
        @test chs[1].colour == 3
        @test mod1(chs[1].k1 + 1, chs[1].d) == chs[1].k2      # cyclically adjacent

        @test _fr_braid_back(g) === nothing                    # the base rule does not apply
        @test _circular_dissolve_through(g, 1, 2) !== nothing
        @test _circular_pull_trivalent(g, 1, 2) === nothing         # nc = 2 ⇒ not responsible

        # The pass-through node (`[1,1,3,3]`, i.e. four-armed) remains standing —
        # the channel version of C7 fires directly on it.
        name, terms = _bc_fires(g)
        @test name === :braid_back
        @test length(terms) == 2                               # id + pinch
        @test sort([length(t.nodes) for t in terms]) == [0, 4]
        _bc_pruefe(g, terms)
    end

    @testset "nc = 2 with THREE foreign arms ⇒ edge + monochrome node" begin
        # Beside the channel edge a genuine node [1,1,1]
        # remains standing. Without the generalisation this diagram would be a
        # CIRCULAR_RULES fixed point.
        g = _bc_passthrough_three()
        chs = _circular_braid_channels(g, 1, 2)
        @test length(chs) == 1
        @test chs[1].nc == 2 && chs[1].d == 5
        @test mod1(chs[1].k1 + 1, chs[1].d) == chs[1].k2

        @test _fr_braid_back(g) === nothing                  # the base rule does not apply
        @test _circular_pull_trivalent(g, 1, 2) === nothing       # nc = 2 ⇒ not responsible

        h = _circular_dissolve_through(g, 1, 2)
        @test h !== nothing
        @test letters(h.word) == letters(g.word)
        @test circular_degree(h) == circular_degree(g)
        @test euler(h) == 2 && isempty(check_wiring(h))
        @test any(nd -> arms(nd) == [1, 1, 1], h.nodes)      # the node left standing
        @test count(nd -> nd.kind === :braid, h.nodes) == 2  # the braids are untouched

        name, terms = _bc_fires(g)
        @test name === :braid_back
        @test length(terms) == 2                             # id + pinch
        @test sort([length(t.nodes) for t in terms]) == [1, 5]
        _bc_pruefe(g, terms)
    end

    @testset "nc ≥ 3 (general-13) ⇒ split off a trivalent, then C8" begin
        for (nm, g, rest) in (("[1,1,3,3,3] (3,4)", _bc_gen13_a(), [3, 3, 1, 1]),
                              ("[1,1,3,3,3] (4,5)", _bc_gen13_b(), [3, 1, 1, 3]),
                              ("[1,3,1,3,3] crossing", _bc_gen13_crossing(), [3, 1, 3, 1]))
            chs = _circular_braid_channels(g, 1, 2)
            @test length(chs) == 1
            @test chs[1].nc == 3

            @test _fr_braid_relation(g) === nothing            # the base rule does not apply
            @test _circular_dissolve_through(g, 1, 2) === nothing   # nc ≠ 2

            # The split gives the PURE trivalent [3,3,3] plus the remainder.
            h = _circular_pull_trivalent(g, 1, 2)
            @test h !== nothing
            @test letters(h.word) == letters(g.word)
            @test circular_degree(h) == circular_degree(g)
            @test euler(h) == 2 && isempty(check_wiring(h))
            @test any(nd -> arms(nd) == [3, 3, 3], h.nodes)

            name, terms = _bc_fires(g)
            @test name === :braid_relation
            @test length(terms) == 1
            # braid + trivalent(2) + the trimmed node
            kinds = sort([arms(nd) for nd in terms[1].nodes])
            @test sort(rest) in [sort(a) for a in kinds]
            @test any(nd -> nd.kind === :braid, terms[1].nodes)
            @test any(nd -> arms(nd) == [2, 2, 2], terms[1].nodes)
            _bc_pruefe(g, terms)
        end
    end

    @testset "13-nodes in BOTH outer channels" begin
        g = _bc_two_sided()
        chs = _circular_braid_channels(g, 1, 2)
        @test length(chs) == 2
        @test sort([c.nc for c in chs]) == [2, 3]

        @test _fr_braid_relation(g) === nothing
        # Here too the four-armed pass-through
        # `[1,1,3,3]` remains standing and C8 fires directly.
        name, terms = _bc_fires(g)
        @test name === :braid_relation
        @test length(terms) == 1
        @test any(nd -> arms(nd) == [3, 1, 3, 1], terms[1].nodes)
        _bc_pruefe(g, terms)
    end

    @testset "12-braids: the rule is colour-independent (1 ↔ 3)" begin
        # Everything above was computed in the {2,3} case. Under the A3 automorphism
        # 1↔3 every 23-braid becomes a 12-braid; the rule
        # name, the term count and the terms themselves must follow.
        for (nm, g) in (("A1", _bc_classic_back()), ("B1", _bc_classic_rel()),
                        ("A2", _bc_passthrough()),      ("A3", _bc_passthrough_three()),
                        ("B3a", _bc_gen13_a()),       ("B3b", _bc_gen13_b()),
                        ("B3c", _bc_gen13_crossing()),   ("two-sided", _bc_two_sided()))
            gs = _bc_swap13(g)
            @test euler(gs) == 2
            @test isempty(check_wiring(gs))

            n23, t23 = _bc_fires(g)
            n12, t12 = _bc_fires(gs)
            @test n12 === n23
            @test length(t12) == length(t23)
            @test sort([circular_canonical_key(_bc_swap13(t)) for t in t23]) ==
                  sort([circular_canonical_key(t) for t in t12])
            _bc_pruefe(gs, t12)
        end
    end

    @testset "the connecting slots must be cyclically adjacent" begin
        # No channel node ⇒ no normalisation, and the rule behaves exactly like the
        # base one. (A braid itself is barred as a channel node.)
        g = _bc_classic_back()
        @test isempty(_circular_braid_channels(g, 1, 2))
        @test _circular_pull_trivalent(g, 1, 2) === nothing
        @test _circular_dissolve_through(g, 1, 2) === nothing
    end
end

