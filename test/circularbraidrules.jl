# test/circularbraidrules.jl — the two-colour braid rules C6–C8, the arm_count == 6
# guard, the :braid node kind, and the C7/R5 pinch term port accuracy.

@testset "C6–C8 — the two-colour rules" begin
    N = DiagrammaticHecke.Node
    NP = DiagrammaticHecke.NodePort

    # ---- C6 _fr_dot_into_braid (R4) — colour symmetry, the pair {2,3}. ---------
    nodesr4 = [N(:braid, [2, 3], 3), N(:dot, [2], 0)]
    edgesr4 = [
        Edge(2, NP(1, 1), NP(2, 1)),
        Edge(3, Leaf(1), NP(1, 2)),
        Edge(2, Leaf(2), NP(1, 3)),
        Edge(3, Leaf(3), NP(1, 4)),
        Edge(2, Leaf(4), NP(1, 5)),
        Edge(3, Leaf(5), NP(1, 6)),
    ]
    gr4 = WordGraph(CircularWord([3, 2, 3, 2, 3]), nodesr4, edgesr4)
    r4_new = DiagrammaticHecke._fr_dot_into_braid(circular(gr4))
    @test r4_new !== nothing
    @test length(r4_new) == 2

    # ---- slot counting at a braid SELF-LOOP.
    # Braid [2,3,2,3,2,3]: dot(2) on slot 1, self-loop(3) slot 2 <-> 6, slots 3/4/5
    # to leaves. All 6 slots are occupied, but only 5 EDGES hang on the braid — which
    # is why `_circular_edges_at_node` gives 5 while `_circular_slot_ports` gives 6.
    # Wired PLANAR — shared fixture test/fixtures_circular.jl, identical to `gsl` in
    # test/circularregionrules.jl.
    gsl = _planar_braid_selfloop_fixture()
    @test euler(gsl) == 2
    @test isempty(check_wiring(gsl))
    @test gsl.nodes[1].kind === :braid
    @test length(DiagrammaticHecke._circular_edges_at_node(gsl, 1)) == 5      # edges
    @test length(DiagrammaticHecke._circular_slot_ports(gsl, 1)) == 6         # slots
    @test sort(collect(keys(DiagrammaticHecke._circular_slot_ports(gsl, 1)))) == collect(1:6)
    # C6 leaves the case alone by an explicit self-loop condition; C1 (`_fr_needle`)
    # is responsible.
    @test DiagrammaticHecke._fr_dot_into_braid(gsl) === nothing

    # A self-loop at a braid is 0: it is only possible with a dot in between, and
    # then it is a pitchfork. The only planar situation is exactly this one (dot on
    # slot 1, self-loop slot2<->slot6) — C1 `_fr_needle` matches directly (braid-like
    # nodes included) and the driver reduces to 0.
    @test isempty(DiagrammaticHecke._fr_needle(gsl))
    (c_sl, hist_sl) = reduce_circular(gsl)
    @test isempty(collect(pairs_of(c_sl)))            # 0 (empty sum)
    @test hist_sl == [:needle]

    # ---- C7 _fr_braid_back (R5) — two back-to-back m=3 braids. ------------------
    # Wired PLANAR — shared fixture test/fixtures_circular.jl, the same as the
    # `circular(...)` version in test/circularregionrules.jl.
    gr5 = _planar_braid_back_wordgraph()
    @test euler(gr5) == 2
    @test isempty(check_wiring(gr5))
    r5_new = DiagrammaticHecke._fr_braid_back(circular(gr5))
    @test r5_new !== nothing
    @test length(r5_new) == 2

    # ---- C8 _fr_braid_relation (R9) — two m=3 braids + one pure trivalent node.
    #      Same node/edge topology as the R9 fixture in test/diagrams.jl, but wired
    #      PLANAR — shared fixture test/fixtures_circular.jl. ---
    gr9 = _planar_braid_relation_wordgraph()
    @test euler(gr9) == 2
    @test isempty(check_wiring(gr9))
    r9_new = DiagrammaticHecke._fr_braid_relation(circular(gr9))
    @test r9_new !== nothing
    gnew9 = only(DiagrammaticHecke.pairs_of(r9_new))[1]
    @test length(gnew9.nodes) == 2   # 3 ↦ 2, a genuine simplification

    # negative case: the K1 pre-form (two braids, no trivalent) must NOT fire.
    k1 = WordGraph(CircularWord([1, 2, 2, 1]),
        [N(:braid, [1, 2], 3), N(:braid, [1, 2], 3)],
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(2, Leaf(2), NP(1, 2)),
             Edge(2, Leaf(3), NP(2, 2)), Edge(1, Leaf(4), NP(2, 1)),
             Edge(1, NP(1, 3), NP(2, 3)), Edge(1, NP(1, 5), NP(2, 5)),
             Edge(2, NP(1, 4), NP(2, 4)), Edge(2, NP(1, 6), NP(2, 6))])
    @test DiagrammaticHecke._fr_braid_relation(circular(k1)) === nothing

    # Regression for the same finding in the circular stack: the cluster IN THE MIDDLE
    # in a larger diagram — a collision double leaf — must leave no
    # loose leaves and must keep the rest of the diagram.
    loose_leaves(gg) = begin
        n = length(DiagrammaticHecke.letters(gg.word)); belegt = falses(n)
        for e in gg.edges, pp in (e.a, e.b)
            pp isa DiagrammaticHecke.Leaf && 1 <= pp.k <= n && (belegt[pp.k] = true)
        end
        [k for k in 1:n if !belegt[k]]
    end
    fbig = circular(double_leaf([1, 2, 1, 3, 2, 1], [0, 1, 0, 1, 1, 1],
                           [2, 1, 3, 2, 3, 1], [1, 0, 1, 1, 0, 1]))
    rfbig = DiagrammaticHecke._fr_braid_relation(fbig)
    @test rfbig !== nothing
    hfbig = first(DiagrammaticHecke.pairs_of(rfbig))[1]
    @test length(hfbig.nodes) == length(fbig.nodes) - 1
    @test isempty(loose_leaves(hfbig))
    @test DiagrammaticHecke.is_wired(hfbig)
    @test all(DiagrammaticHecke._circular_leaf_colour(hfbig, k) != 0
              for k in 1:length(DiagrammaticHecke.letters(hfbig.word)))

    # ---- C7b _fr_braid_back --------------------------------------------------
    # A wrong triC wiring in the pinch
    # term when the inner sector sits in the middle of the slot numbers. This
    # end-to-end test runs the concrete case (x,y,e,f) = (132312,323123,011010,111101)
    # through reduce_circular + apply_circular_d4 and must not throw.
    x = [1, 3, 2, 3, 1, 2]
    y = [3, 2, 3, 1, 2, 3]
    e = [0, 1, 1, 0, 1, 0]
    f = [1, 1, 1, 1, 0, 1]
    t = only(tt for tt in double_leaves(x, y) if tt.e == e && tt.f == f)
    fm = circular_morphism(t.morphism)
    fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
    faktor, fdm0 = DiagrammaticHecke.circular_extract_scalars(fdm)
    combo, _ = reduce_circular(DiagrammaticHecke._circular_decorated(fdm0))
    combo = faktor * combo
    terms = collect(DiagrammaticHecke.pairs_of(combo))
    term2 = terms[findfirst(((fd, _),) -> find_circular_d4_match(CircularMorphismGraph(fd.graph, fm.cut1, fm.cut2)) !== nothing, terms)]
    fd2, _ = term2
    mt = find_circular_d4_match(CircularMorphismGraph(fd2.graph, fm.cut1, fm.cut2))
    @test mt !== nothing
    # the call must not throw a KeyError on the region index.
    @test begin apply_circular_d4(fd2, mt...); true end
end

# ---- C6–C9 do NOT fire on 8-armed braid nodes --------------------------------
#
# The explicit guard `arm_count == 6` keeps the eight-armed nodes away from C6–C9.
# This block pins the behaviour down without reference to the node kind.

"The bare `2k`-armed braid node with colour pair `(s,t)` (leaves backwards)."
function _p11_bare(n::Int, s::Int = 1, t::Int = 2)
    a = [isodd(i) ? s : t for i in 1:n]
    lf(i) = i == 1 ? 1 : n + 2 - i
    word = [a[i] for i in 1:n][invperm([lf(i) for i in 1:n])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:n])
end

"Dot on boundary leaf `k` (like `_c14_attach_dot` in test/circulardotongen12.jl)."
function _p11_dot_at(g::CircularGraph, k::Int)
    lets = DiagrammaticHecke.letters(g.word); c = lets[k]
    nodes = vcat(copy(g.nodes), CircularNode[circular_node([c])]); di = length(nodes)
    shift(p) = p isa Leaf ? (p.k == k ? NodePort(di, 1) : Leaf(p.k > k ? p.k - 1 : p.k)) : p
    CircularGraph(CircularWord(vcat(lets[1:k-1], lets[k+1:end])), nodes,
             Edge[Edge(e.colour, shift(e.a), shift(e.b)) for e in g.edges])
end

@testset "C6–C9 with the explicit arm_count == 6 guard" begin
    # ⚠ C14 depends on the dot policy; here we test the
    # surgery itself, not the gate (that's tested in circulardotmerge.jl).
    _p11_old_policy = CIRCULAR_DOT_POLICY[]
    CIRCULAR_DOT_POLICY[] = false
    g6 = _p11_bare(6)
    g8 = _p11_bare(8)
    @test arm_count(g6.nodes[1]) == 6 && arm_count(g8.nodes[1]) == 8
    @test isempty(check_wiring(g8)) && euler(g8) == 2

    # C6: a dot at each of the eight arms — no match. Responsible instead is C14
    # (`dot_on_gen12`; shut down in the driver, where the generic preimage stage
    # takes over), and it gives a 6-armed braid plus a trivalent.
    for k in 1:8
        gd = _p11_dot_at(g8, k)
        @test DiagrammaticHecke._fr_dot_into_braid(gd) === nothing
        r = DiagrammaticHecke._fr_dot_on_gen12(gd)
        @test r !== nothing
        h = first(DiagrammaticHecke.pairs_of(r))[1]
        @test sort([arm_count(nd) for nd in h.nodes]) == [3, 6]
    end
    # Control check on the 6-armed node: C6 does indeed fire there.
    @test DiagrammaticHecke._fr_dot_into_braid(_p11_dot_at(g6, 1)) !== nothing
    CIRCULAR_DOT_POLICY[] = _p11_old_policy
    # and with the policy on: ONE dot on the 8-armed node is a normal form
    @test DiagrammaticHecke._fr_dot_on_gen12(_p11_dot_at(g8, 1)) === nothing

    # C7/C8: the candidate lists must not take up the 8-armed node.
    # The two-node case (two 8-armed nodes, three shared edges — the
    # C16 pattern) is in test/circulargen12ongen12.jl, where the planar fixture lives.
    @test DiagrammaticHecke._fr_braid_back(g8) === nothing
    @test DiagrammaticHecke._fr_braid_relation(g8) === nothing
end

@testset "ONE kind :braid, carrying the arm count as an attribute" begin
    for n in (6, 8, 10, 12, 14)
        nd = circular_node([isodd(i) ? 1 : 2 for i in 1:n])
        @test nd.kind === :braid && arm_count(nd) == n
        @test DiagrammaticHecke._circular_braidlike(nd)
    end
    # the invariants hold for EVERY arm count; odd / too small / wrong colours or
    # non-alternating still throw.
    @test_throws ArgumentError CircularNode(:braid, [1, 2, 1])              # 3 arms
    @test_throws ArgumentError CircularNode(:braid, [1, 2, 1, 2])           # 4 arms
    @test_throws ArgumentError CircularNode(:braid, [1, 2, 1, 2, 1, 2, 1])  # odd
    @test_throws ArgumentError CircularNode(:braid, [1, 3, 1, 3, 1, 3, 1, 3])   # colour pair
    @test_throws ArgumentError CircularNode(:braid, [1, 2, 1, 2, 2, 1, 1, 2])   # not alternating

    # `_circular_braidlike` is literally `kind === :braid`.
    @test !DiagrammaticHecke._circular_braidlike(circular_node([1, 1, 1]))
    @test !DiagrammaticHecke._circular_braidlike(circular_node([2, 2, 2, 2]))
    @test !DiagrammaticHecke._circular_braidlike(circular_node([1, 3, 1, 3]))

    @test DiagrammaticHecke._ck_kind_id(:braid) == 3

    # `.fwg` round trip: the kind is written as `braid`, the arms follow it.
    g8io = CircularGraph(CircularWord([isodd(i) ? 1 : 2 for i in 1:8]),
                    CircularNode[circular_node([isodd(i) ? 1 : 2 for i in 1:8])],
                    Edge[Edge(isodd(i) ? 1 : 2, Leaf(i == 1 ? 1 : 10 - i), NodePort(1, i))
                         for i in 1:8])
    txt = circulargraph_to_string(g8io)
    @test occursin("n braid 1,2,1,2,1,2,1,2", txt)
    @test circulargraph_from_string(txt) == g8io
end

@testset "port accuracy of the C7/R5 pinch term" begin
    # `legsB_tri_nested` must take B's slot order, not A's traversal order — a
    # mistake the node-kind/colour comparison of the C6–C9 tests does not catch.
    # Here: the exact port wiring of the pinch term, on a fixture whose inner sector
    # does NOT sit at the slot-number boundary — `_planar_braid_back_wordgraph`
    # (test/fixtures_circular.jl) has that shape: node 2's internal slots are
    # {2,3,4}, NOT {1,..} or {..,6}.
    gr5b = _planar_braid_back_wordgraph()
    fgr5b = circular(gr5b)
    NP = DiagrammaticHecke.NodePort
    internal2 = Set(e.b.slot for e in fgr5b.edges
                     if e.a isa NP && e.a.node == 1 && e.b isa NP && e.b.node == 2)
    @test internal2 == Set([2, 3, 4])   # inner sector at node 2, not on the boundary

    pinch_of(pairs) = pairs[argmax([length(t.nodes) for (t, _) in pairs])][1]

    # ---- C7 (circular): the pinch term is planar and does NOT change under rotation of
    # node 2's slot numbering (rotate_arms represents the SAME diagram — a wrong
    # sector order gives different `circular_canonical_key`s per rotation here).
    pinch_circular = pinch_of(DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_braid_back(fgr5b)))
    @test euler(pinch_circular) == 2 && isempty(check_wiring(pinch_circular))
    fgr5b_rot = rotate_arms(fgr5b, 2, 1)   # shifts the inner sector at node 2
    pinch_circular_rot = pinch_of(DiagrammaticHecke.pairs_of(DiagrammaticHecke._fr_braid_back(fgr5b_rot)))
    @test euler(pinch_circular_rot) == 2 && isempty(check_wiring(pinch_circular_rot))
    @test circular_canonical_key(pinch_circular) == circular_canonical_key(pinch_circular_rot)


end

