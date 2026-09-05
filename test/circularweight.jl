# test/circularweight.jl — circular_weight deltas per rule C1–C9; the delta
# table is in the file header of src/circular/CircularWeight.jl.

@testset "circular_weight" begin
    NP = DiagrammaticHecke.NodePort
    N = DiagrammaticHecke.Node
    sc(nd, s) = DiagrammaticHecke._slot_colour(nd, s)

    # Mini-check 1: component ordering — a braid node weighs more than a
    # dot, although the dot has fewer nodes/arms (the braid count dominates lexicographically).
    @test circular_weight(circular(dot(1))) < circular_weight(circular(braid(1, 2; m = 3)))

    # Mini-check 2: circular_weight(CircularDecorated) == circular_weight(underlying graph).
    g_dec = circular(trivalent(1))
    @test circular_weight(circular_decorated(g_dec)) == circular_weight(g_dec)

    # C1 _fr_needle: self-connection ↦ 0 (empty CircularComboR) — trivially falling,
    # no after-term to compare; only the before-fixture is anchored here.
    nd_self = circular_node([1, 1, 1])
    g_self = CircularGraph(CircularWord([1]), [nd_self],
        Edge[Edge(1, NP(1, 1), NP(1, 2)), Edge(1, Leaf(1), NP(1, 3))])
    r1 = DiagrammaticHecke._fr_needle(g_self)
    @test r1 !== nothing && isempty(r1)

    # C2 _fr_dot_walks_on: [1,3,3,1,3,3]+dot(1) ↦ [3,3,3,3]+fresh dot.
    nd2 = circular_node([1, 3, 3, 1, 3, 3])
    before2 = CircularGraph(CircularWord([1, 3, 3, 3, 3]), [nd2, circular_node([1])],
        Edge[
            Edge(1, NP(1, 1), NP(2, 1)), Edge(3, Leaf(1), NP(1, 2)),
            Edge(3, Leaf(2), NP(1, 3)), Edge(1, Leaf(3), NP(1, 4)),
            Edge(3, Leaf(4), NP(1, 5)), Edge(3, Leaf(5), NP(1, 6)),
        ])
    ((after2, _),) = DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_dot_walks_on(before2))
    @test circular_weight(after2) < circular_weight(before2)

    # C3 _fr_barbell: two degree-1 nodes of the same colour ↦ empty graph + α_i.
    g3 = CircularGraph(CircularWord(Int[]), [circular_node([1]), circular_node([1])],
        Edge[Edge(1, NP(1, 1), NP(2, 1))])
    ((d3, _),) = DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_barbell(circular_decorated(g3)))
    @test circular_weight(d3) < circular_weight(g3)

    # C4 _fr_bead: the bead [1,1] disappears, its two edges splice.
    g4 = CircularGraph(CircularWord([1, 1]), [DiagrammaticHecke._unchecked_circularnode(:mixed, [1, 1])],
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2))])
    ((after4, _),) = DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_bead(g4))
    @test circular_weight(after4) < circular_weight(g4)

    # C5 _fr_merge: [1,1,1]+dot ↦ bead collapses immediately (splice fallback), 0 nodes.
    g5 = CircularGraph(CircularWord([1, 1]), [circular_node([1, 1, 1]), circular_node([1])],
        Edge[Edge(1, NP(1, 3), NP(2, 1)), Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2))])
    ((after5, _),) = DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_merge(g5))
    @test circular_weight(after5) < circular_weight(g5)

    # C6 _fr_dot_into_braid: braid(2,3;m=3)+dot(2) ↦ sum of two terms, BOTH
    # strictly smaller (term A even raises the node count / arm sum, but the braid count goes 1↦0).
    gr4 = WordGraph(CircularWord([3, 2, 3, 2, 3]),
        [N(:braid, [2, 3], 3), N(:dot, [2], 0)],
        Edge[Edge(2, NP(1, 1), NP(2, 1)), Edge(3, Leaf(1), NP(1, 2)),
             Edge(2, Leaf(2), NP(1, 3)), Edge(3, Leaf(3), NP(1, 4)),
             Edge(2, Leaf(4), NP(1, 5)), Edge(3, Leaf(5), NP(1, 6))])
    fg6 = circular(gr4)
    r6 = DiagrammaticHecke._fr_dot_into_braid(fg6)
    @test all(circular_weight(gg) < circular_weight(fg6) for (gg, _) in DiagrammaticHecke.pairs_of(r6))

    # C7 _fr_braid_back: two back-to-back m=3 braids(1,2) ↦ sum of two terms.
    nodesr5 = [N(:braid, [1, 2], 3), N(:braid, [1, 2], 3)]
    edgesr5 = Edge[]
    for (s, lf) in [(1, 1), (2, 2), (3, 3)]
        push!(edgesr5, Edge(sc(nodesr5[1], s), Leaf(lf), NP(1, s)))
    end
    for s in [4, 5, 6]
        push!(edgesr5, Edge(sc(nodesr5[1], s), NP(1, s), NP(2, s)))
    end
    for (s, lf) in [(1, 4), (2, 5), (3, 6)]
        push!(edgesr5, Edge(sc(nodesr5[2], s), Leaf(lf), NP(2, s)))
    end
    fg7 = circular(WordGraph(CircularWord([1, 2, 1, 1, 2, 1]), nodesr5, edgesr5))
    r7 = DiagrammaticHecke._fr_braid_back(fg7)
    @test all(circular_weight(gg) < circular_weight(fg7) for (gg, _) in DiagrammaticHecke.pairs_of(r7))

    # C8 _fr_braid_relation: two m=3 braids + one pure trivalent ↦ 1 braid + 1 trivalent.
    fg8 = circular(WordGraph(CircularWord([2, 1, 2, 1, 2, 1, 1]),
        [N(:braid, [2, 1], 3), N(:trivalent, [1], 0), N(:braid, [1, 2], 3)],
        Edge[Edge(2, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2)), Edge(2, Leaf(3), NP(1, 3)),
             Edge(1, NP(1, 4), NP(2, 1)), Edge(1, Leaf(4), NP(2, 2)), Edge(1, NP(1, 6), NP(3, 1)),
             Edge(2, NP(1, 5), NP(3, 2)), Edge(1, NP(2, 3), NP(3, 3)), Edge(2, Leaf(7), NP(3, 6)),
             Edge(1, Leaf(6), NP(3, 5)), Edge(2, Leaf(5), NP(3, 4))]))
    ((after8, _),) = DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_braid_relation(fg8))
    @test circular_weight(after8) < circular_weight(fg8)

end

