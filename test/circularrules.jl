# test/circularrules.jl — the reduction rules C1–C5, C10, C12 (needle, dot-walk,
# barbell, bead, merge, mono-double, commutation-merge, two-adjacent-merge) and
# their R11 subsumption.

@testset "_fr_needle (C1) — self-loop / free circle ↦ 0" begin
    NP = DiagrammaticHecke.NodePort
    # free circle.
    fgloop = CircularGraph(CircularWord(Int[]), CircularNode[], Edge[Edge(1, Circle(1), Circle(1))])
    r = DiagrammaticHecke._fr_needle(fgloop)
    @test r !== nothing && isempty(r)

    # self-loop at a non-:braid node.
    nd = circular_node([1, 1, 1])
    g_self = CircularGraph(CircularWord([1]), [nd],
        Edge[Edge(1, NP(1, 1), NP(1, 2)), Edge(1, Leaf(1), NP(1, 3))])
    r2 = DiagrammaticHecke._fr_needle(g_self)
    @test r2 !== nothing && isempty(r2)

    # colour symmetry (1↔3, swapc).
    swapc(c) = c == 1 ? 3 : c == 3 ? 1 : c
    nd3 = circular_node([swapc(1), swapc(1), swapc(1)])
    g_self3 = CircularGraph(CircularWord([swapc(1)]), [nd3],
        Edge[Edge(swapc(1), NP(1, 1), NP(1, 2)), Edge(swapc(1), Leaf(1), NP(1, 3))])
    r3 = DiagrammaticHecke._fr_needle(g_self3)
    @test r3 !== nothing && isempty(r3)

    # negative case: no match on an ordinary dot.
    @test DiagrammaticHecke._fr_needle(circular(dot(1))) === nothing

    # A braid self-loop IS 0 — a conservative exception for it would be
    # gone.
    fgb = circular(braid(1, 2; m = 3))
    ndb = fgb.nodes[1]
    g_braid_self = CircularGraph(CircularWord(Int[]), [ndb, circular_node([1, 1, 1])],
        Edge[Edge(1, NP(1, 1), NP(1, 4)), Edge(2, NP(1, 2), NP(1, 5)),
             Edge(1, NP(1, 3), NP(2, 1)), Edge(1, NP(1, 6), NP(2, 2))])
    # the only relevant point: the (1,4) edge is a self-loop at the :braid — it
    # matches and gives 0.
    rbs = DiagrammaticHecke._fr_needle(g_braid_self)
    @test rbs !== nothing && isempty(rbs)
end

@testset "_fr_dot_walks_on (C2) — the case split and the ordering trap" begin
    NP = DiagrammaticHecke.NodePort

    # main case: [1,3,3,1,3,3] + 1-dot ↦ [3,3,3,3] + dot, is_wired on both.
    nd = circular_node([1, 3, 3, 1, 3, 3])
    nd_dot = circular_node([1])
    before = CircularGraph(CircularWord([1, 3, 3, 3, 3]), [nd, nd_dot],
        Edge[
            Edge(1, NP(1, 1), NP(2, 1)),
            Edge(3, Leaf(1), NP(1, 2)),
            Edge(3, Leaf(2), NP(1, 3)),
            Edge(1, Leaf(3), NP(1, 4)),
            Edge(3, Leaf(4), NP(1, 5)),
            Edge(3, Leaf(5), NP(1, 6)),
        ])
    @test is_wired(before)
    r = DiagrammaticHecke._fr_dot_walks_on(before)
    @test r !== nothing
    @test length(r) == 1
    ((after, coeff),) = DiagrammaticHecke.pairs_of(r)
    @test coeff == one(SoergelPoly)
    @test is_wired(after)
    kinds_after = sort([sort(arms(n)) for n in after.nodes])
    @test kinds_after == sort([[3, 3, 3, 3], [1]])
    # the new dot hangs on a fresh degree-1 node with exactly one edge.
    freshdots = [ni for (ni, n) in enumerate(after.nodes) if arms(n) == [1]]
    @test length(freshdots) == 1
    @test length(DiagrammaticHecke._circular_edges_at_node(after, only(freshdots))) == 1

    # THE ORDERING TRAP: on the R2-1 form [1,3,1,3]+dot(1), _fr_dot_walks_on gives a
    # clean result. The [3,3] remainder is itself not a valid bead (CircularNode
    # forbids degree 2), so _fr_dot_walks_on collapses it IMMEDIATELY (as _fr_bead
    # would): the two remaining neighbours are joined directly, no intermediate node
    # is left standing — only the fresh dot.
    #
    # `_fr_merge` THROWS here: `[3,1,3]` is a two-coloured node with only ONE 1,
    # i.e. an invalid `CircularNode`. The order C2 before C5 is therefore the
    # condition for the driver to run through at all.
    nd21 = circular_node([1, 3, 1, 3])
    g21 = CircularGraph(CircularWord([3, 1, 3]), [nd21, circular_node([1])],
        Edge[Edge(1, NP(1, 1), NP(2, 1)), Edge(3, Leaf(1), NP(1, 2)),
             Edge(1, Leaf(2), NP(1, 3)), Edge(3, Leaf(3), NP(1, 4))])
    r_walk = DiagrammaticHecke._fr_dot_walks_on(g21)
    @test r_walk !== nothing
    @test_throws ArgumentError DiagrammaticHecke._fr_merge(g21)
    gw = only(DiagrammaticHecke.pairs_of(r_walk))[1]
    @test is_wired(gw)
    walk_kinds = sort([sort(arms(n)) for n in gw.nodes])
    @test walk_kinds == [[1]]                          # clean: only the fresh dot remains
    @test DiagrammaticHecke._circular_leaf_colour(gw, 1) == 3 && DiagrammaticHecke._circular_leaf_colour(gw, 3) == 3

    # THE [1,1,3,3] CASE: here the two slots of colour c are NEIGHBOURS, not
    # opposite — one remainder block is empty, the other carries both remaining arms.
    # A bead branch assuming "one remaining arm per side" would build the
    # invariant-violating `[1,1,3]` node. Required behaviour: one continuous line,
    # and the dot carries on — a 1-edge Leaf(1)–Leaf(2) plus a 3-dot at Leaf(3).
    g1133 = CircularGraph(CircularWord([1, 1, 3]), [circular_node([1, 3, 3, 1]), circular_node([3])],
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 4)),
             Edge(3, Leaf(3), NP(1, 3)), Edge(3, NP(2, 1), NP(1, 2))])
    @test is_wired(g1133)
    expected = CircularGraph(CircularWord([1, 1, 3]), [circular_node([3])],
        Edge[Edge(1, Leaf(1), Leaf(2)), Edge(3, Leaf(3), NP(1, 1))])
    c1133, h1133 = reduce_circular(g1133)
    @test h1133 == [:dot_walks_on]
    g1133_out = only(DiagrammaticHecke.pairs_of(c1133))[1].graph
    @test circular_canonical_key(g1133_out) == circular_canonical_key(expected)
    @test is_wired(g1133_out) && euler(g1133_out) == 2

    # and the same at each of the four arms, for both rotations of the arm sequence:
    # no path lands at `merge`, each ends at the clean dot.
    for armseq in ([1, 3, 3, 1], [1, 1, 3, 3]), slot in 1:4
        c = armseq[slot]
        rest = [k for k in 1:4 if k != slot]
        es = Edge[Edge(c, NP(2, 1), NP(1, slot))]
        for (l, k) in enumerate(rest)
            push!(es, Edge(armseq[k], Leaf(l), NP(1, k)))
        end
        gk = CircularGraph(CircularWord(armseq[rest]), [circular_node(armseq), circular_node([c])], es)
        @test is_wired(gk)
        ck, hk = reduce_circular(gk)
        @test hk == [:dot_walks_on]
        @test [arms(n) for n in only(DiagrammaticHecke.pairs_of(ck))[1].graph.nodes] == [[c]]
    end

    # negative case: count(==(c), arms) == 3 (not 2) must NOT go through this rule.
    nd_neg = circular_node([1, 1, 1, 3, 3])
    g_neg = CircularGraph(CircularWord([1, 1, 3, 3]), [nd_neg, circular_node([1])],
        Edge[Edge(1, NP(1, 1), NP(2, 1)), Edge(1, Leaf(1), NP(1, 2)),
             Edge(1, Leaf(2), NP(1, 3)), Edge(3, Leaf(3), NP(1, 4)),
             Edge(3, Leaf(4), NP(1, 5))])
    @test DiagrammaticHecke._fr_dot_walks_on(g_neg) === nothing

    # colour symmetry: the same shape with 1↔3 swapped.
    nd_sw = circular_node([3, 1, 1, 3, 1, 1])
    g_sw = CircularGraph(CircularWord([3, 1, 1, 1, 1]), [nd_sw, circular_node([3])],
        Edge[
            Edge(3, NP(1, 1), NP(2, 1)),
            Edge(1, Leaf(1), NP(1, 2)),
            Edge(1, Leaf(2), NP(1, 3)),
            Edge(3, Leaf(3), NP(1, 4)),
            Edge(1, Leaf(4), NP(1, 5)),
            Edge(1, Leaf(5), NP(1, 6)),
        ])
    r_sw = DiagrammaticHecke._fr_dot_walks_on(g_sw)
    @test r_sw !== nothing
    ((after_sw, _),) = DiagrammaticHecke.pairs_of(r_sw)
    @test is_wired(after_sw)
    @test sort([sort(arms(n)) for n in after_sw.nodes]) == sort([[1, 1, 1, 1], [3]])
end

@testset "_fr_barbell (C3) — factor α_i, MUST fire before merge" begin
    NP = DiagrammaticHecke.NodePort
    nd1 = circular_node([1]); nd2 = circular_node([1])
    g = CircularGraph(CircularWord(Int[]), [nd1, nd2], Edge[Edge(1, NP(1, 1), NP(2, 1))])
    fd = circular_decorated(g)
    r = DiagrammaticHecke._fr_barbell(fd)
    @test r !== nothing
    @test length(r) == 1
    ((d, coeff),) = DiagrammaticHecke.pairs_of(r)
    # the barbell lies in the OUTER face here (empty boundary word). `α_i` is placed
    # IN PLACE into the `outer_label` and the coefficient stays `1`. The same is
    # checked in `test/circulardriver.jl` at the driver end.
    @test coeff == one(SoergelPoly)
    @test d.outer_label == alpha(1)
    @test isempty(d.graph.nodes)
    @test isempty(d.graph.edges)

    # colour symmetry: colour 3.
    nd1c = circular_node([3]); nd2c = circular_node([3])
    gc = CircularGraph(CircularWord(Int[]), [nd1c, nd2c], Edge[Edge(3, NP(1, 1), NP(2, 1))])
    rc = DiagrammaticHecke._fr_barbell(circular_decorated(gc))
    @test rc !== nothing

    # negative case: no barbell if one node is NOT of degree 1.
    nd_triv = circular_node([1, 1, 1])
    g_no = CircularGraph(CircularWord([1, 1]), [nd1, nd_triv],
        Edge[Edge(1, NP(1, 1), NP(2, 1)), Edge(1, Leaf(1), NP(2, 2)), Edge(1, Leaf(2), NP(2, 3))])
    @test DiagrammaticHecke._fr_barbell(circular_decorated(g_no)) === nothing

    # _fr_merge must NOT fire on a genuine barbell (degree1+degree1 is _fr_barbell's
    # business) — exactly the ordering guarantee.
    @test DiagrammaticHecke._fr_merge(g) === nothing
end

@testset "_fr_bead (C4) — a degree-2 bead disappears" begin
    NP = DiagrammaticHecke.NodePort
    # `circular_node([1,1])` throws under the invariant (beads are not valid
    # CircularNodes, see CircularGraph.jl `circular_node`), but this test checks the
    # collapse path `_fr_bead` that removes such beads, so the bead is built with
    # the unchecked inner constructor as an intermediate state.
    nd = CircularNode(:mixed, [1, 1])
    g = CircularGraph(CircularWord([1, 1]), [nd],
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2))])
    r = DiagrammaticHecke._fr_bead(g)
    @test r !== nothing
    ((after, coeff),) = DiagrammaticHecke.pairs_of(r)
    @test coeff == one(SoergelPoly)
    @test isempty(after.nodes)
    @test only(after.edges) == Edge(1, Leaf(1), Leaf(2))

    # Farbsymmetrie.
    nd3 = CircularNode(:mixed, [3, 3])
    g3 = CircularGraph(CircularWord([3, 3]), [nd3],
        Edge[Edge(3, Leaf(1), NP(1, 1)), Edge(3, Leaf(2), NP(1, 2))])
    r3 = DiagrammaticHecke._fr_bead(g3)
    @test r3 !== nothing

    # negative case: arm_count == 3 must not match.
    @test DiagrammaticHecke._fr_bead(circular(trivalent(1))) === nothing

    # embedded in context: a bead in the middle of a longer chain, again built with
    # the unchecked constructor (see above).
    nd_ctx = CircularNode(:mixed, [1, 1])
    nd_far = circular_node([1, 1, 1])
    g_ctx = CircularGraph(CircularWord([1, 1, 1]), [nd_ctx, nd_far],
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, NP(1, 2), NP(2, 1)),
             Edge(1, Leaf(2), NP(2, 2)), Edge(1, Leaf(3), NP(2, 3))])
    r_ctx = DiagrammaticHecke._fr_bead(g_ctx)
    @test r_ctx !== nothing
    ((after_ctx, _),) = DiagrammaticHecke.pairs_of(r_ctx)
    @test length(after_ctx.nodes) == 1
    @test is_wired(after_ctx)
end

@testset "_fr_merge (C5) — arm deletion in the general branch" begin
    NP = DiagrammaticHecke.NodePort
    # [1,1,1] + dot ↦ [1,1] (the unit rule as a special case) — but under the
    # "no degree-2 node" invariant the resulting bead [1,1] collapses IMMEDIATELY in
    # the splice fallback of `_merge_at` (CircularGraph.jl): `circular_node([1,1])` throws, so
    # `_merge_at` takes the same path as `_fr_bead` and wires the two neighbours
    # straight through — 0 nodes remain instead of a single [1,1] node.
    nd_triv = circular_node([1, 1, 1]); nd_dot = circular_node([1])
    g = CircularGraph(CircularWord([1, 1]), [nd_triv, nd_dot],
        Edge[Edge(1, NP(1, 3), NP(2, 1)), Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2))])
    r = DiagrammaticHecke._fr_merge(g)
    @test r !== nothing
    ((after, coeff),) = DiagrammaticHecke.pairs_of(r)
    @test coeff == one(SoergelPoly)
    @test isempty(after.nodes)
    @test only(after.edges) == Edge(1, Leaf(1), Leaf(2))
    @test is_wired(after)

    # On [1,3,3,1,3,3]+dot the dot colour occurs TWICE; the GENERAL branch needs
    # count != 2 (see the C2 test for the two-arm case), so the fixture is
    # [1,1,1,3,3] with count(==(1)) == 3.
    nd_gen = circular_node([1, 1, 1, 3, 3])
    g_gen = CircularGraph(CircularWord([1, 1, 3, 3]), [nd_gen, circular_node([1])],
        Edge[Edge(1, NP(1, 1), NP(2, 1)), Edge(1, Leaf(1), NP(1, 2)),
             Edge(1, Leaf(2), NP(1, 3)), Edge(3, Leaf(3), NP(1, 4)),
             Edge(3, Leaf(4), NP(1, 5))])
    @test DiagrammaticHecke._fr_dot_walks_on(g_gen) === nothing   # not the two-arm case
    r_gen = DiagrammaticHecke._fr_merge(g_gen)
    ((after_gen, _),) = DiagrammaticHecke.pairs_of(r_gen)
    @test arms(after_gen.nodes[1]) == [1, 1, 3, 3]   # the dot arm (slot 1) is deleted

    # negative case: no edge between two distinct non-:braid nodes.
    @test DiagrammaticHecke._fr_merge(circular(dot(1))) === nothing
    @test DiagrammaticHecke._fr_merge(circular(braid(1, 2; m = 3))) === nothing
end

@testset "_fr_merge (C5) — MERGE LOCK for mixed 1/3 nodes" begin
    NP = DiagrammaticHecke.NodePort
    # Two nodes joined by EXACTLY ONE edge; the remaining arms go to leaves.
    # `a`/`b` are the arm sequences, slot 1 is the connection in each.
    function joined(a, b)
        na, nb = circular_node(a), circular_node(b)
        a[1] == b[1] || error("connection colours must match")
        w = vcat(a[2:end], b[2:end])
        es = Edge[Edge(a[1], NP(1, 1), NP(2, 1))]
        for s in 2:length(a); push!(es, Edge(a[s], Leaf(s - 1), NP(1, s))); end
        for s in 2:length(b); push!(es, Edge(b[s], Leaf(length(a) - 1 + s - 1), NP(2, s))); end
        return CircularGraph(CircularWord(w), [na, nb], es)
    end

    # ---- LOCKED: both nodes carry 1 AND 3, a single connection ---------------
    # This is the case that must NOT be merged: merging makes it indistinguishable
    # whether the 1-strands really connect through.
    g_block = joined([1, 3, 1, 3], [1, 3, 1, 3])
    @test DiagrammaticHecke._circular_is_bicoloured(g_block.nodes[1])
    @test DiagrammaticHecke._circular_is_bicoloured(g_block.nodes[2])
    @test DiagrammaticHecke._circular_conn_count(g_block, 1, 2) == 1
    @test DiagrammaticHecke._fr_merge(g_block) === nothing
    # and the full driver leaves it standing too (no detour via another rule onto
    # the same merged node)
    ((fix, _),) = DiagrammaticHecke.pairs_of(reduce_circular_full(g_block))
    @test length(fix.graph.nodes) == 2

    # "carries 1 AND 3" and "pure 1/3 node with both colours" are the same thing:
    # `circular_node` does not let colour 2 mix with 1/3 at all. The case where the two
    # formulations could differ is thus not constructible — recorded here so the lock
    # is not later extended to "mixed colour-2 nodes", which do not exist.
    @test_throws ArgumentError circular_node([1, 3, 2, 3])

    # ---- ALLOWED: at least one of the two is monochrome ----------------------
    for (a, b) in ([[1, 1, 1], [1, 3, 1, 3]],     # trivalent against crossing
                   [[1, 3, 1, 3], [1, 1, 1]],     # different order
                   [[3, 3, 3], [3, 1, 3, 1]],
                   [[2, 2, 2], [2, 2, 2]])        # no 1/3 world at all
        g_ok = joined(a, b)
        @test DiagrammaticHecke._fr_merge(g_ok) !== nothing
    end

    # ---- ALLOWED: a double connection stays untouched ------------------------
    # Two mixed nodes with TWO connections: the lock does not apply because it only
    # concerns the SINGLE connection. (In the driver C12/C10 take this case before C5
    # anyway; here C5 is checked directly.)
    n1, n2 = circular_node([1, 3, 1, 3]), circular_node([3, 1, 3, 1])
    g_two = CircularGraph(CircularWord([1, 3, 3, 1]), [n1, n2],
        Edge[Edge(1, NP(1, 1), NP(2, 2)), Edge(3, NP(1, 2), NP(2, 1)),
             Edge(1, Leaf(1), NP(1, 3)), Edge(3, Leaf(2), NP(1, 4)),
             Edge(3, Leaf(3), NP(2, 3)), Edge(1, Leaf(4), NP(2, 4))])
    @test DiagrammaticHecke._circular_conn_count(g_two, 1, 2) == 2
    @test DiagrammaticHecke._fr_merge(g_two) !== nothing
end

@testset "_fr_merge (C5) — the co-oriented {1,3} double connection" begin
    NP = DiagrammaticHecke.NodePort
    # THE MERGE CONDITION: one connection merges only
    # if at least one node contains just a single colour. Two connections: a)
    # different colours ⇒ merge · b) same colours ⇒ needle = 0.
    # Sentence 1 = the merge lock in the testset above. Sentence 2 belongs to
    # `_fr_mono_double` (b) and C12/C10 (a) — all three BEFORE C5 in the registry.
    #
    # THE GAP: two mixed 1/3 nodes, EXACTLY TWO connections of DIFFERENT colours,
    # neighbouring at both nodes, wired CO-ORIENTED. C12 and C10 require
    # counter-orientation, `_fr_mono_double` equal colours: NONE of the three
    # matches, and C5 clears the situation alone. Locking C5 wholesale at
    # `_circular_conn_count >= 2` would leave exactly this diagram standing as a
    # fixed point, so it is NOT locked (see the block comment before `_fr_merge`).
    n1, n2 = circular_node([1, 3, 1, 3]), circular_node([1, 3, 1, 3])
    g_same = CircularGraph(CircularWord([1, 3, 1, 3]), [n1, n2],
        Edge[Edge(1, NP(1, 1), NP(2, 1)), Edge(3, NP(1, 2), NP(2, 2)),
             Edge(1, Leaf(1), NP(1, 3)), Edge(3, Leaf(2), NP(1, 4)),
             Edge(1, Leaf(3), NP(2, 3)), Edge(3, Leaf(4), NP(2, 4))])
    @test DiagrammaticHecke._circular_conn_count(g_same, 1, 2) == 2
    @test DiagrammaticHecke._circular_is_bicoloured(g_same.nodes[1])
    @test DiagrammaticHecke._circular_is_bicoloured(g_same.nodes[2])
    # co-oriented: A slot 1 ↔ B slot 1 and A slot 2 ↔ B slot 2 (counter-oriented
    # would be 1↔2 / 2↔1, the C12 case in the testset below).
    @test DiagrammaticHecke._fr_mono_double(g_same) === nothing        # colours differ
    @test DiagrammaticHecke._fr_two_adjacent_merge(g_same) === nothing # not counter-oriented
    @test DiagrammaticHecke._fr_commutation_merge(g_same) === nothing  # ditto
    @test DiagrammaticHecke._fr_merge(g_same) !== nothing              # C5 applies on its own

    # ... and NO other registry rule applies to this diagram — the fixed-point claim,
    # pinned down mechanically here.
    fd_same = DiagrammaticHecke.circular_decorated(g_same)
    for r in DiagrammaticHecke.CIRCULAR_RULES
        r.name === :merge && continue
        @test r.apply(fd_same) === nothing
    end

    # ---- second class: a MONOCHROME double connection at a 3-ARMED node.
    # This is sentence 2b ("equal colours ⇒ needle = 0"). `_fr_mono_double` bounds
    # the arm count at `>= 3` (any valid non-dot arm count), so it is responsible
    # itself and the history is `[:mono_double]`. Result: 0.
    a6, b3 = circular_node(fill(2, 6)), circular_node([2, 2, 2])
    g_mono3 = CircularGraph(CircularWord(fill(2, 5)), [a6, b3],
        Edge[Edge(2, NP(1, 1), NP(2, 2)), Edge(2, NP(1, 2), NP(2, 1)),
             Edge(2, Leaf(1), NP(1, 3)), Edge(2, Leaf(2), NP(1, 4)),
             Edge(2, Leaf(3), NP(1, 5)), Edge(2, Leaf(4), NP(1, 6)),
             Edge(2, Leaf(5), NP(2, 3))])
    @test DiagrammaticHecke._circular_conn_count(g_mono3, 1, 2) == 2
    @test DiagrammaticHecke._fr_mono_double(g_mono3) !== nothing   # 3-armed node — a hit
    @test isempty(DiagrammaticHecke._fr_mono_double(g_mono3))      # ⇒ 0
    @test DiagrammaticHecke._fr_two_adjacent_merge(g_mono3) === nothing
    @test DiagrammaticHecke._fr_commutation_merge(g_mono3) === nothing
    @test isempty(DiagrammaticHecke.pairs_of(reduce_circular_full(g_mono3)))   # ⇒ 0, as in sentence 2b
    (c_mono3, hist_mono3) = reduce_circular(g_mono3)
    @test isempty(collect(pairs_of(c_mono3)))
    @test hist_mono3 == [:mono_double]
end

@testset "_fr_two_adjacent_merge (C12) — EXACTLY TWO adjacent connections" begin
    NP = DiagrammaticHecke.NodePort

    function build_two_adjacent(armsA, k, armsB, l)
        n, m = length(armsA), length(armsB)
        ndA = DiagrammaticHecke._unchecked_circularnode(:mixed, armsA)
        ndB = DiagrammaticHecke._unchecked_circularnode(:mixed, armsB)
        k2 = mod1(k + 1, n); l2 = mod1(l + 1, m)
        edges = Edge[
            Edge(armsA[k], NP(1, k), NP(2, l2)),
            Edge(armsA[k2], NP(1, k2), NP(2, l)),
        ]
        leaf = 0
        for s in 1:n
            (s == k || s == k2) && continue
            leaf += 1
            push!(edges, Edge(armsA[s], Leaf(leaf), NP(1, s)))
        end
        for s in 1:m
            (s == l || s == l2) && continue
            leaf += 1
            push!(edges, Edge(armsB[s], Leaf(leaf), NP(2, s)))
        end
        return CircularGraph(CircularWord(fill(1, max(leaf, 1))), [ndA, ndB], edges)
    end

    # ---- MONOCHROME case (1,1): not C12's business — C12 matches DIFFERENT
    # colours only. The monochrome case is `_fr_mono_double` (⇒ 0, the monochrome
    # bigon). C10 does not see the case either (Set != (1,3)). It is example 1 in
    # the `_fr_mono_double` test block below.
    armsA = [1, 1, 3, 3, 3]; armsB = [1, 1, 3, 3, 3]
    g = build_two_adjacent(armsA, 1, armsB, 1)
    @test DiagrammaticHecke._fr_commutation_merge(g) === nothing     # C10 does not see this case
    @test DiagrammaticHecke._fr_two_adjacent_merge(g) === nothing    # C12: monochrome, not a match
    r0 = DiagrammaticHecke._fr_mono_double(g)
    @test r0 !== nothing
    @test isempty(pairs_of(r0))                                  # the diagram is 0

    # through the full driver: `_fr_mono_double` fires BEFORE C12 (the ordering trap,
    # CircularDriver.jl header) and pulls the diagram to 0.
    names = [r.name for r in DiagrammaticHecke.CIRCULAR_RULES]
    @test !(:adjacent_pair_out in names)
    @test findfirst(==(:mono_double), names) < findfirst(==(:two_adjacent_merge), names)
    @test findfirst(==(:two_adjacent_merge), names) < findfirst(==(:commutation_merge), names)

    fd = circular_decorated(g)
    c, steps = reduce_circular(fd; maxsteps = 1)
    @test steps == [:mono_double]
    @test isempty(pairs_of(c))                                  # 0

    # ---- different colours (overlaps C10, same result — pure overlap, no
    # competition).
    armsA2 = [1, 3, 1, 1]; armsB2 = [3, 1, 1, 1]
    g2 = build_two_adjacent(armsA2, 1, armsB2, 1)
    r10_2 = DiagrammaticHecke._fr_commutation_merge(g2)
    r12_2 = DiagrammaticHecke._fr_two_adjacent_merge(g2)
    @test r10_2 !== nothing && r12_2 !== nothing
    ((a10_2, _),) = pairs_of(r10_2)
    ((a12_2, _),) = pairs_of(r12_2)
    minrot(v) = minimum([circshift(v, -k) for k in 0:(length(v) - 1)])
    @test minrot(arms(a10_2.nodes[1])) == minrot(arms(a12_2.nodes[1]))

    # ---- negative case: only ONE connection ⇒ no match (that is the
    # single-connection case, _fr_merge/_merge_at's business, not C12's). [1,1,3]
    # itself violates the "every colour at least twice" invariant (only one 3) —
    # invalid as a STANDALONE node, but exactly the intermediate state merge_at_edge
    # is about to splice (the same trick as in the merge_nodes splice testset
    # above).
    nd_single_i = DiagrammaticHecke._unchecked_circularnode(:mixed, [1, 1, 3])
    nd_single_j = DiagrammaticHecke._unchecked_circularnode(:mixed, [3, 1, 3])
    g_single = CircularGraph(CircularWord([1, 1, 1, 1, 3, 3]), [nd_single_i, nd_single_j],
        Edge[Edge(3, NP(1, 3), NP(2, 1)), Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2)),
             Edge(1, Leaf(3), NP(2, 2)), Edge(3, Leaf(4), NP(2, 3))])
    @test DiagrammaticHecke._fr_two_adjacent_merge(g_single) === nothing

    # ---- negative case: THREE connections ⇒ no match (stays needle/merge's
    # business, see the C10 test block below for the three-strand case).
    n3 = [Node(:braid, [1, 3], 2), Node(:braid, [3, 1], 2)]
    sc(nd, s) = DiagrammaticHecke._slot_colour(nd, s)
    e3 = Edge[Edge(sc(n3[1], 3), NP(1, 3), NP(2, 4)), Edge(sc(n3[1], 4), NP(1, 4), NP(2, 3)),
              Edge(sc(n3[1], 2), NP(1, 2), NP(2, 2)),
              Edge(sc(n3[1], 1), Leaf(1), NP(1, 1)), Edge(sc(n3[2], 1), Leaf(2), NP(2, 1))]
    g3 = circular(WordGraph(CircularWord([1, 3]), n3, e3))
    @test DiagrammaticHecke._fr_two_adjacent_merge(g3) === nothing

    # ---- negative case: wired co-oriented ⇒ no match (merge_at_edges' ArgumentError
    # is caught internally, no error escapes — the rule simply gives `nothing`).
    function build_two_same(armsA, k, armsB, l)
        n, m = length(armsA), length(armsB)
        ndA = DiagrammaticHecke._unchecked_circularnode(:mixed, armsA)
        ndB = DiagrammaticHecke._unchecked_circularnode(:mixed, armsB)
        k2 = mod1(k + 1, n); l2 = mod1(l + 1, m)
        edges = Edge[Edge(armsA[k], NP(1, k), NP(2, l)), Edge(armsA[k2], NP(1, k2), NP(2, l2))]
        leaf = 0
        for s in 1:n
            (s == k || s == k2) && continue
            leaf += 1
            push!(edges, Edge(armsA[s], Leaf(leaf), NP(1, s)))
        end
        for s in 1:m
            (s == l || s == l2) && continue
            leaf += 1
            push!(edges, Edge(armsB[s], Leaf(leaf), NP(2, s)))
        end
        return CircularGraph(CircularWord(fill(1, max(leaf, 1))), [ndA, ndB], edges)
    end
    g_same = build_two_same(copy(armsA), 1, copy(armsB), 1)
    @test DiagrammaticHecke._fr_two_adjacent_merge(g_same) === nothing

    # ---- negative case: no match on a bare dot / :braid.
    @test DiagrammaticHecke._fr_two_adjacent_merge(circular(dot(1))) === nothing
    @test DiagrammaticHecke._fr_two_adjacent_merge(circular(braid(1, 2; m = 3))) === nothing
end

@testset "_fr_mono_double — a MONOCHROME double connection ↦ 0" begin
    NP = DiagrammaticHecke.NodePort

    function build_two_adjacent(armsA, k, armsB, l)
        n, m = length(armsA), length(armsB)
        ndA = DiagrammaticHecke._unchecked_circularnode(:mixed, armsA)
        ndB = DiagrammaticHecke._unchecked_circularnode(:mixed, armsB)
        k2 = mod1(k + 1, n); l2 = mod1(l + 1, m)
        edges = Edge[
            Edge(armsA[k], NP(1, k), NP(2, l2)),
            Edge(armsA[k2], NP(1, k2), NP(2, l)),
        ]
        leaf = 0
        for s in 1:n
            (s == k || s == k2) && continue
            leaf += 1
            push!(edges, Edge(armsA[s], Leaf(leaf), NP(1, s)))
        end
        for s in 1:m
            (s == l || s == l2) && continue
            leaf += 1
            push!(edges, Edge(armsB[s], Leaf(leaf), NP(2, s)))
        end
        return CircularGraph(CircularWord(fill(1, max(leaf, 1))), [ndA, ndB], edges)
    end

    # ---- example 1 from the C12 block comment (CircularRules.jl): 5+5 arms, monochrome
    # (1,1) ⇒ 0. C12 must NOT give a merged [3,3,3,3,3,3] node here — the trap this
    # rule fixes.
    g1 = build_two_adjacent([1, 1, 3, 3, 3], 1, [1, 1, 3, 3, 3], 1)
    @test DiagrammaticHecke._fr_two_adjacent_merge(g1) === nothing
    r1 = DiagrammaticHecke._fr_mono_double(g1)
    @test r1 !== nothing
    @test isempty(pairs_of(r1))

    # ---- example 3 from the C12 block comment: 4+4 arms, monochrome (1,1) ⇒ 0.
    # C12 must not give [3,3,3,3] here.
    g3 = build_two_adjacent([1, 1, 3, 3], 1, [1, 1, 3, 3], 1)
    @test DiagrammaticHecke._fr_two_adjacent_merge(g3) === nothing
    r3 = DiagrammaticHecke._fr_mono_double(g3)
    @test r3 !== nothing
    @test isempty(pairs_of(r3))

    # ---- example 2 from the C12 block comment: DIFFERENT colours ⇒
    # `_fr_mono_double` does NOT match, C12 stays responsible (merge, unchanged).
    g2 = build_two_adjacent([1, 3, 1, 1, 3], 1, [3, 1, 1, 3, 1], 1)
    @test DiagrammaticHecke._fr_mono_double(g2) === nothing
    r2 = DiagrammaticHecke._fr_two_adjacent_merge(g2)
    @test r2 !== nothing
    ((a2, c2),) = pairs_of(r2)
    @test c2 == one(SoergelPoly)
    @test length(a2.nodes) == 1
    @test is_wired(a2)

    # ---- through the full driver: the example-1 diagram reduces to 0.
    fd = circular_decorated(g1)
    c, steps = reduce_circular(fd; maxsteps = 1)
    @test steps == [:mono_double]
    @test isempty(pairs_of(c))
    cfull = reduce_circular_full(g1)
    @test isempty(pairs_of(cfull))

    # ---- negative case: only ONE connection ⇒ no match.
    nd_single_i = DiagrammaticHecke._unchecked_circularnode(:mixed, [1, 1, 3])
    nd_single_j = DiagrammaticHecke._unchecked_circularnode(:mixed, [3, 1, 3])
    g_single = CircularGraph(CircularWord([1, 1, 1, 1, 3, 3]), [nd_single_i, nd_single_j],
        Edge[Edge(3, NP(1, 3), NP(2, 1)), Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2)),
             Edge(1, Leaf(3), NP(2, 2)), Edge(3, Leaf(4), NP(2, 3))])
    @test DiagrammaticHecke._fr_mono_double(g_single) === nothing

    # ---- negative case: THREE connections ⇒ no match.
    n3 = [Node(:braid, [1, 3], 2), Node(:braid, [3, 1], 2)]
    sc(nd, s) = DiagrammaticHecke._slot_colour(nd, s)
    e3 = Edge[Edge(sc(n3[1], 3), NP(1, 3), NP(2, 4)), Edge(sc(n3[1], 4), NP(1, 4), NP(2, 3)),
              Edge(sc(n3[1], 2), NP(1, 2), NP(2, 2)),
              Edge(sc(n3[1], 1), Leaf(1), NP(1, 1)), Edge(sc(n3[2], 1), Leaf(2), NP(2, 1))]
    g3n = circular(WordGraph(CircularWord([1, 3]), n3, e3))
    @test DiagrammaticHecke._fr_mono_double(g3n) === nothing

    # ---- negative case: no match on a bare dot / :braid.
    @test DiagrammaticHecke._fr_mono_double(circular(dot(1))) === nothing
    @test DiagrammaticHecke._fr_mono_double(circular(braid(1, 2; m = 3))) === nothing
end

@testset "_fr_commutation_merge (C10)" begin
    NP = DiagrammaticHecke.NodePort
    minrot(v) = minimum([circshift(v, -k) for k in 0:(length(v) - 1)])

    # Specification: A = (i₁…iₙ a b), B = (b a j₁…jₘ), with a and b connected
    # COUNTER-ORIENTED (clockwise at one node, counter-clockwise at the other).
    # Result: ONE node (i₁…iₙ j₁…jₘ), i and j keep their wiring. No case split by
    # arm count.
    #
    # Builds exactly that situation: the pair sits on A at slots 1,2 and on B at
    # slots 2,1 (i.e. the other way round), the rest goes to leaves.
    function pairjoined(aA, aB)
        A = circular_node(aA); B = circular_node(aB)
        edges = Edge[Edge(1, NP(1, 1), NP(2, 2)), Edge(3, NP(1, 2), NP(2, 1))]
        lf = 0; w = Int[]
        for s in 3:length(aA)
            lf += 1; push!(w, arm_colour(A, s)); push!(edges, Edge(arm_colour(A, s), Leaf(lf), NP(1, s)))
        end
        for s in 3:length(aB)
            lf += 1; push!(w, arm_colour(B, s)); push!(edges, Edge(arm_colour(B, s), Leaf(lf), NP(2, s)))
        end
        CircularGraph(CircularWord(w), [A, B], edges)
    end

    # slot → leaf, to check the WIRING (not just the arm colours: two nodes can have
    # the same arm sequence and still hang differently).
    function wiring(g)
        out = Tuple{Int,Int,Int}[]
        for e in g.edges
            if e.a isa Leaf && e.b isa NP
                push!(out, (e.b.slot, e.a.k, e.colour))
            elseif e.b isa Leaf && e.a isa NP
                push!(out, (e.a.slot, e.b.k, e.colour))
            end
        end
        sort(out)
    end

    # ---- (13)1313 + (31)3131 ⇒ 31311313 --------------------------------------
    g66 = pairjoined([1, 3, 1, 3, 1, 3], [3, 1, 3, 1, 3, 1])
    r66 = DiagrammaticHecke._fr_commutation_merge(g66)
    @test r66 !== nothing
    ((a66, c66),) = DiagrammaticHecke.pairs_of(r66)
    @test c66 == one(SoergelPoly)
    @test length(a66.nodes) == 1
    @test arm_count(a66.nodes[1]) == 8
    @test minrot(arms(a66.nodes[1])) == minrot([3, 1, 3, 1, 1, 3, 1, 3])
    # the i and j still hang on the same leaves as before, in order
    @test wiring(a66) == [(k, k, arms(a66.nodes[1])[k]) for k in 1:8]

    # ---- 4 + >4: one node with the arm count of the larger one ---------------
    g46 = pairjoined([1, 3, 1, 3], [3, 1, 3, 1, 3, 1])
    ((a46, _),) = DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_commutation_merge(g46))
    @test length(a46.nodes) == 1
    @test arm_count(a46.nodes[1]) == 6
    @test wiring(a46) == [(k, k, arms(a46.nodes[1])[k]) for k in 1:6]

    # ---- 4 + 4: C10 gives a 1133 node -----------------------------------------
    g44 = pairjoined([1, 3, 1, 3], [3, 1, 3, 1])
    ((a44, _),) = DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_commutation_merge(g44))
    @test length(a44.nodes) == 1
    @test minrot(arms(a44.nodes[1])) == [1, 1, 3, 3]      # NOT the crossing
    # THE 1133 NODE STAYS. It is not dissolved into two pass-through strands:
    # 13-nodes may remain together even when all colours are consecutive.
    #
    # ⚠️ The normal form here is the NODE, not the identity. Two parallel edges of
    # colour 1 and 3 that do not touch and are connected with matching colours are
    # merged — that is P1's direction (`circular_parallel_merge_step`), and it is why
    # a pass-through rule dissolving the node has no place in `CIRCULAR_RULES`.
    r44 = reduce_circular_full(g44)
    @test length(DiagrammaticHecke.pairs_of(r44)) == 1
    ((a44f, _),) = DiagrammaticHecke.pairs_of(r44)
    @test length(a44f.graph.nodes) == 1
    @test minrot(arms(a44f.graph.nodes[1])) == [1, 1, 3, 3]
    @test is_wired(a44f.graph)

    # ---- orientation: connected CO-ORIENTED ⇒ C10 does not apply -------------
    # (slots 3↔3 and 4↔4 instead of 3↔4/4↔3 — then it is not (ab)/(ba).)
    sc(nd, s) = DiagrammaticHecke._slot_colour(nd, s)
    nsame = [Node(:braid, [1, 3], 2), Node(:braid, [1, 3], 2)]
    esame = Edge[Edge(sc(nsame[1], 3), NP(1, 3), NP(2, 3)),
                 Edge(sc(nsame[1], 4), NP(1, 4), NP(2, 4)),
                 Edge(sc(nsame[1], 1), Leaf(1), NP(1, 1)), Edge(sc(nsame[1], 2), Leaf(2), NP(1, 2)),
                 Edge(sc(nsame[2], 1), Leaf(3), NP(2, 1)), Edge(sc(nsame[2], 2), Leaf(4), NP(2, 2))]
    gsame = circular(WordGraph(CircularWord([1, 3, 1, 3]), nsame, esame))
    @test DiagrammaticHecke._fr_commutation_merge(gsame) === nothing

    # ---- THREE connections ⇒ one remains as a self-loop ⇒ needle ⇒ 0 ---------
    n3 = [Node(:braid, [1, 3], 2), Node(:braid, [3, 1], 2)]
    e3 = Edge[Edge(sc(n3[1], 3), NP(1, 3), NP(2, 4)), Edge(sc(n3[1], 4), NP(1, 4), NP(2, 3)),
              Edge(sc(n3[1], 2), NP(1, 2), NP(2, 2)),
              Edge(sc(n3[1], 1), Leaf(1), NP(1, 1)), Edge(sc(n3[2], 1), Leaf(2), NP(2, 1))]
    g3 = circular(WordGraph(CircularWord([1, 3]), n3, e3))
    @test isempty(DiagrammaticHecke.pairs_of(reduce_circular_full(g3)))

    # ---- rule order ----------------------------------------------------------
    # `:commutation_merge` (C10) must come before `:merge`, and there is no
    # `:adjacent_pair_out` rule.
    names = [r.name for r in DiagrammaticHecke.CIRCULAR_RULES]
    @test !(:adjacent_pair_out in names)
    @test findfirst(==(:commutation_merge), names) < findfirst(==(:merge), names)
end

@testset "R11 subsumption — bigon → merge → self-loop → needle" begin
    NP = DiagrammaticHecke.NodePort
    c = 1
    nd1 = circular_node([c, c, c]); nd2 = circular_node([c, c, c])
    edges_bigon = Edge[
        Edge(c, Leaf(1), NP(1, 1)),
        Edge(c, Leaf(2), NP(2, 1)),
        Edge(c, NP(1, 2), NP(2, 2)),
        Edge(c, NP(1, 3), NP(2, 3)),
    ]
    gb = CircularGraph(CircularWord([c, c]), [nd1, nd2], edges_bigon)
    m = merge_at_edge(gb, 3)
    @test length(m.nodes) == 1
    @test all(==(c), arms(m.nodes[1]))
    @test arm_count(m.nodes[1]) == 4
    e_surv = only(e for e in m.edges if e.a isa NodePort && e.b isa NodePort)
    @test e_surv.a.node == e_surv.b.node
    r = DiagrammaticHecke._fr_needle(m)
    @test r !== nothing
    @test isempty(r)
end

