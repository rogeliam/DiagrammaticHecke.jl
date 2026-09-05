# test/circularregionrules.jl — boundary words (circular/CircularBoundaryWord.jl) and the
# region-based rules (circular/rules/CircularRegionRules.jl).
#
# ⚠️ PLANARITY: region rules read the face tracer, so their fixtures need a PLANAR
# embedding (`euler == 2`). The C7/C8 fixtures in test/circular.jl are not planar
# (gr5 −4, gr9 −4 — irrelevant for the algebraic rules, but not for faces). The
# back-to-back fixture here is therefore rewired (verified by brute force over the
# slot assignments, euler == 2); the planar fixture serves as the C8 oracle.

# Two back-to-back m=3 braids, wired PLANAR (euler == 2): the a-legs run
# the leaf ring backwards (Faces.jl §consequences), internal edges
# a4–b4/a5–b3/a6–b2, b-legs to leaves 4–6.
function _planar_braid_back_fixture()
    NP = DiagrammaticHecke.NodePort
    N = DiagrammaticHecke.Node
    nodes = [N(:braid, [1, 2], 3), N(:braid, [1, 2], 3)]
    edges = Edge[
        Edge(1, Leaf(1), NP(1, 3)), Edge(2, Leaf(2), NP(1, 2)), Edge(1, Leaf(3), NP(1, 1)),
        Edge(2, NP(1, 4), NP(2, 4)), Edge(1, NP(1, 5), NP(2, 3)), Edge(2, NP(1, 6), NP(2, 2)),
        Edge(1, Leaf(4), NP(2, 1)), Edge(2, Leaf(5), NP(2, 6)), Edge(1, Leaf(6), NP(2, 5)),
    ]
    return circular(WordGraph(CircularWord([1, 2, 1, 1, 2, 1]), nodes, edges))
end

@testset "circular_boundary_words — extraction and rotation invariance" begin
    NP = DiagrammaticHecke.NodePort

    # Needle fixture (test/circular.jl C1): self-loop at [1,1,1], slot 3 to the leaf.
    nd = circular_node([1, 1, 1])
    g_self = CircularGraph(CircularWord([1]), [nd],
        Edge[Edge(1, NP(1, 1), NP(1, 2)), Edge(1, Leaf(1), NP(1, 3))])
    ws = circular_boundary_words(g_self)
    @test length(ws) == region_count(g_self)
    inner = [w for w in ws if DiagrammaticHecke.is_interior(w)]
    # exactly one inner region: the monogon enclosed by the self-loop.
    @test length(inner) == 1 && length(only(inner)) == 1
    l = only(only(inner).letters)
    @test l.colour == 1 && l.corner.vertex == (:node, 1)

    # Back-to-back braids: 3 internal edges ⇒ two inner BIGON regions with
    # corners at both braids.
    fgr5 = _planar_braid_back_fixture()
    @test euler(fgr5) == 2
    bigons = [w for w in circular_boundary_words(fgr5)
              if DiagrammaticHecke.is_interior(w) && length(w) == 2]
    @test length(bigons) == 2
    @test all(w -> Set(l.corner.vertex for l in w.letters) ==
                   Set([(:node, 1), (:node, 2)]), bigons)

    # rotation invariance: rotate_arms changes only the slot numbering, not the diagram
    # — the multiset (length, corner-node set) of the inner boundary words stays the
    # same.
    fingerprint(g) = sort([(length(w), sort(collect(Set(l.corner.vertex[2] for l in w.letters))))
                           for w in circular_boundary_words(g) if DiagrammaticHecke.is_interior(w)])
    fp0 = fingerprint(fgr5)
    @test fingerprint(rotate_arms(fgr5, 1, 2)) == fp0
    @test fingerprint(rotate_arms(fgr5, 2, 4)) == fp0
end

@testset "region rules — oracle equality with C1/C7/C8" begin
    NP = DiagrammaticHecke.NodePort
    N = DiagrammaticHecke.Node

    # C1 oracle: free circle and self-loop ⇒ 0, dot/braid self-loop not.
    fgloop = CircularGraph(CircularWord(Int[]), CircularNode[], Edge[Edge(1, Circle(1), Circle(1))])
    r = DiagrammaticHecke._fr_needle(fgloop)
    @test r !== nothing && isempty(r)
    nd = circular_node([1, 1, 1])
    g_self = CircularGraph(CircularWord([1]), [nd],
        Edge[Edge(1, NP(1, 1), NP(1, 2)), Edge(1, Leaf(1), NP(1, 3))])
    r2 = DiagrammaticHecke._fr_needle(g_self)
    @test r2 !== nothing && isempty(r2)
    @test DiagrammaticHecke._fr_needle(circular(dot(1))) === nothing
    fgb = circular(braid(1, 2; m = 3))
    ndb = fgb.nodes[1]
    g_braid_self = CircularGraph(CircularWord(Int[]), [ndb, circular_node([1, 1, 1])],
        Edge[Edge(1, NP(1, 1), NP(1, 4)), Edge(2, NP(1, 2), NP(1, 5)),
             Edge(1, NP(1, 3), NP(2, 1)), Edge(1, NP(1, 6), NP(2, 2))])
    # Because of planarity, 2 arms of a node can only be self-connected if 1. they are
    # adjacent, 2. everything between them goes to dots. If you traced the dots
    # through, you'd have other nodes whose neighbours are self-connected, which is
    # then mathematically always 0. The loop 1↔4 encloses arms 2,3, and C1 says 0
    # there too — the edge pattern does not care what the loop encloses.
    rbs = DiagrammaticHecke._fr_needle(g_braid_self)
    @test rbs !== nothing && isempty(rbs)

    # C7 oracle: the identical result to _fr_braid_back, also under slot rotation (the
    # shape variant at issue).
    fgr5 = _planar_braid_back_fixture()
    @test DiagrammaticHecke._frr_braid_back(fgr5) == DiagrammaticHecke._fr_braid_back(fgr5)
    @test DiagrammaticHecke._frr_braid_back(fgr5) !== nothing
    fgr5r = rotate_arms(fgr5, 1, 2)
    @test DiagrammaticHecke._frr_braid_back(fgr5r) == DiagrammaticHecke._fr_braid_back(fgr5r)

    # C8 oracle on the PLANAR fixture (a route-A cluster in the middle of the
    # diagram); negative case K1 (two braids, no trivalent).
    fbig = circular(double_leaf([1, 2, 1, 3, 2, 1], [0, 1, 0, 1, 1, 1],
                           [2, 1, 3, 2, 3, 1], [1, 0, 1, 1, 0, 1]))
    @test euler(fbig) == 2
    @test DiagrammaticHecke._frr_braid_relation(fbig) == DiagrammaticHecke._fr_braid_relation(fbig)
    @test DiagrammaticHecke._frr_braid_relation(fbig) !== nothing
    k1 = WordGraph(CircularWord([1, 2, 2, 1]),
        [N(:braid, [1, 2], 3), N(:braid, [1, 2], 3)],
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(2, Leaf(2), NP(1, 2)),
             Edge(2, Leaf(3), NP(2, 2)), Edge(1, Leaf(4), NP(2, 1)),
             Edge(1, NP(1, 3), NP(2, 3)), Edge(1, NP(1, 5), NP(2, 5)),
             Edge(2, NP(1, 4), NP(2, 4)), Edge(2, NP(1, 6), NP(2, 6))])
    @test DiagrammaticHecke._frr_braid_relation(circular(k1)) === nothing
end

@testset "the gap (dot at a braid, neighbouring arms joined) — needle covers it" begin
    NP = DiagrammaticHecke.NodePort
    # Like the C6 self-loop fixture in test/circular.jl but wired PLANAR (the leaves run
    # BACKWARDS through the slots, cw against ccw; the original fixture has euler = 0):
    # braid
    # [2,3,2,3,2,3], dot(2) on slot 1, self-loop(3) slot 2 ↔ 6 —
    # the loop encloses EXACTLY the dot-capped arm.
    #
    # C6 (`_fr_dot_into_braid`) does not fire here; under the convention that a braid
    # self-loop counts as 0, the needle rule alone covers the case completely.
    gsl = CircularGraph(CircularWord([2, 3, 2]),
        [circular_node([2, 3, 2, 3, 2, 3]), circular_node([2])],
        Edge[Edge(2, NP(1, 1), NP(2, 1)), Edge(3, NP(1, 2), NP(1, 6)),
             Edge(2, Leaf(1), NP(1, 5)), Edge(3, Leaf(2), NP(1, 4)),
             Edge(2, Leaf(3), NP(1, 3))])
    @test euler(gsl) == 2
    @test DiagrammaticHecke._fr_dot_into_braid(gsl) === nothing
    rn = DiagrammaticHecke._fr_needle(gsl)
    @test rn !== nothing && isempty(rn)
    @test isempty(reduce_circular_full(gsl))
    @test isempty(reduce_circular_full(gsl; rules = CIRCULAR_RULES_REGION))
end

@testset "C21 pitchfork — the figure C1 needle does not see" begin
    NP = DiagrammaticHecke.NodePort
    # The reduction `121 → 212 → 22 → ε` drawn as a diagram: a 6-armed braid
    # `[1,2,1,2,1,2]` (the 121/212 node), its slot 2 (colour 2) capped by a dot, and
    # its two neighbouring 1-arms — slots 1 and 3 — running to ONE trivalent
    # `[1,1,1]`. That is the pitchfork: a dot between two same-coloured arms that
    # meet again.
    #
    # `_fr_needle` requires both ends of ONE edge on the same node and declines —
    # the two arms end on a DIFFERENT node. C21 is what makes the term 0.
    g = CircularGraph(CircularWord([2, 1, 2, 1]),
        CircularNode[circular_node([1, 2, 1, 2, 1, 2]),
                     circular_node([2]),
                     circular_node([1, 1, 1])],
        Edge[Edge(2, NP(1, 2), NP(2, 1)),
             Edge(1, NP(1, 1), NP(3, 1)),
             Edge(1, NP(1, 3), NP(3, 3)),
             Edge(2, Leaf(1), NP(1, 4)),
             Edge(1, Leaf(2), NP(3, 2)),
             Edge(2, Leaf(3), NP(1, 6)),
             Edge(1, Leaf(4), NP(1, 5))])
    @test euler(g) == 2 && isempty(check_wiring(g))          # a legal diagram
    @test DiagrammaticHecke._fr_needle(g) === nothing        # C1 does not fire
    r = _fr_pitchfork(g)
    @test r !== nothing && isempty(r)                        # C21 gives 0
    @test isempty(reduce_circular(g)[1])
    @test reduce_circular(g)[2] == [:pitchfork]              # and it is C21 that does it
end

@testset "CIRCULAR_RULES_REGION — registry construction and full comparison" begin
    names = [r.name for r in CIRCULAR_RULES_REGION]
    # C1 has no region version: its pattern is an edge pattern, so the region
    # registry takes `_fr_needle` itself, under its own name.
    @test :needle in names
    # C21 pitchfork IS registered, in both registries. It is not covered by C1:
    # `_fr_needle` wants a genuine self-edge (`e.a.node == e.b.node`), while the
    # pitchfork's two same-coloured neighbour arms may end on ANOTHER node —
    # then needle declines and C21 does not (the fixture below).
    @test :pitchfork in names
    @test findfirst(==(:pitchfork), names) < findfirst(==(:dot_walks_on), names)
    @test findfirst(==(:braid_back_region), names) < findfirst(==(:braid_back), names)
    @test findfirst(==(:braid_relation_region), names) < findfirst(==(:braid_relation), names)
    # The "big braid node" rules stay LAST (ordering trap, CircularDriver.jl):
    # building a large braid node is the "only when nothing cheaper is left"
    # move. `gen12_merge` (C13) is followed by `trivalent_into_gen12` (C19) and
    # `mixed_into_gen12` (C23), which sit next to it for the same reason.
    @test names[end - 2] === :gen12_merge
    @test names[end - 1] === :trivalent_into_gen12
    @test names[end] === :mixed_into_gen12
    fnames = [r.name for r in CIRCULAR_RULES]
    @test :pitchfork in fnames

    # Full comparison on both fixtures: both registries give the same
    # normal form.
    fgr5 = _planar_braid_back_fixture()
    @test reduce_circular_full(fgr5) == reduce_circular_full(fgr5; rules = CIRCULAR_RULES_REGION)
    fbig = circular(double_leaf([1, 2, 1, 3, 2, 1], [0, 1, 0, 1, 1, 1],
                           [2, 1, 3, 2, 3, 1], [1, 0, 1, 1, 0, 1]))
    @test reduce_circular_full(fbig) == reduce_circular_full(fbig; rules = CIRCULAR_RULES_REGION)
end

