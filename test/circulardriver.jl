# test/circulardriver.jl — the circular reduction driver. circular_weight deltas
# per rule live in test/circularweight.jl, not here. Fixtures are slot-exact
# (diagram/Graph.jl §PLANAR SLOT CONVENTION).

@testset "reduce_circular — driver basics" begin
    NP = DiagrammaticHecke.NodePort

    # Needle in context: a self-loop at a [1,1,1] trivalent, WITH a genuine boundary
    # leaf beside it (context, not the isolated case from the C1 rule fixture in
    # test/circular.jl) -> the whole combination is 0.
    nd = circular_node([1, 1, 1])
    g_self = CircularGraph(CircularWord([1]), [nd],
        Edge[Edge(1, NP(1, 1), NP(1, 2)), Edge(1, Leaf(1), NP(1, 3))])
    combo0, hist0 = reduce_circular(circular_decorated(g_self))
    @test isempty(combo0)
    @test hist0 == [:needle]

    # Barbell in context: two degree-1 nodes of colour 1 (a barbell) BESIDE a separate
    # trivalent [2,2,2] with an α₂ label on its only cell. A bare barbell has NO
    # boundary gap of its own (it is always an isolated two-node component away from
    # the boundary, checked structurally) — so the α_i always lands in the COEFFICIENT,
    # never as a cell label; "the right cell" means here: the existing α₂ label on the
    # surviving trivalent cell stays UNTOUCHED while the barbell
    # vanishes cleanly.
    nd1 = circular_node([1]); nd2 = circular_node([1])
    nd3 = circular_node([2, 2, 2])
    g_barb = CircularGraph(CircularWord([2, 2, 2]), [nd1, nd2, nd3],
        # ORIENTATION (src/diagram/Faces.jl): boundary ccw, arms cw — a
        # a node at several leaves runs the leaf ring BACKWARDS. Written forwards
        # (Leaf(k) at slot k) the embedding is not planar, the regions collapse from 4
        # to 2 and the construction throws right here.
        Edge[Edge(1, NP(1, 1), NP(2, 1)),
             Edge(2, Leaf(1), NP(3, 1)), Edge(2, Leaf(3), NP(3, 2)), Edge(2, Leaf(2), NP(3, 3))])
    fd_barb = circular_decorated(g_barb, [one(SoergelPoly), alpha(2), one(SoergelPoly), one(SoergelPoly)])
    combo1, hist1 = reduce_circular(fd_barb)
    @test hist1 == [:barbell]
    @test length(combo1) == 1
    ((d1, coeff1),) = pairs_of(combo1)
    # The barbell factor α₁ sits in the `outer_label`, not in the coefficient.
    # The phantom region around the barbell has no boundary gap and dies
    # with it; after the deletion it is nowhere addressable any more (the inner anchor
    # does not help either — `_circular_delete_nodes` renumbers), so it takes the route
    # outwards. The VALUE is the same, only the place it is stored differs.
    @test coeff1 == one(SoergelPoly)
    @test d1.outer_label == alpha(1)                  # the barbell itself
    @test length(d1.graph.nodes) == 1                  # only the trivalent remains
    @test face_count(d1) == 3                           # 3 boundary leaves, 3 sector cells
    # The labels are REGION-indexed: the α₂ sits on the region with gap 2, and exactly
    # that one survives the barbell's disappearance unchanged.
    #
    # The region is addressed via its BOUNDARY GAP, not via its number:
    # the numbering follows the cell order and shifts with the orientation — the order
    # of the regions really doesn't matter, what matters is that distances are
    # correct. What actually counts is that EXACTLY the region that had the label
    # before keeps it.
    r_in = 2                                          # that is how fd_barb above is built
    gaps_in = DiagrammaticHecke.regions(g_barb)[r_in].gaps
    @test fd_barb.region_labels[r_in] == alpha(2)     # that's how it was before
    r_out = findfirst(R -> !isempty(intersect(R.gaps, gaps_in)),
                          DiagrammaticHecke.regions(d1.graph))
    @test region_label(d1, r_out) == alpha(2)          # the pre-existing label survives
    @test count(!isone, d1.region_labels) == 1             # and nothing else

    # THE CORE PROMISE: the two bracketings of 111 -> 1 — two trivalents [1,1,1], wired
    # differently (a different slot choice for the connecting edge, with
    # correspondingly rotated boundary legs) but with the same 4 outer leaves. As
    # CircularGraphs they ARE structurally different (`==` false), but after reduce_circular they
    # must arrive at THE SAME circular_canonical_key. The leaf-rotation fix in
    # _circular_leaf_rotations (CircularGraph.jl) is needed exactly for that.
    c = 1
    bracketA = CircularGraph(CircularWord([c, c, c, c]), [circular_node([c, c, c]), circular_node([c, c, c])],
        Edge[Edge(c, Leaf(1), NP(1, 1)), Edge(c, Leaf(2), NP(1, 2)),
             Edge(c, Leaf(3), NP(2, 1)), Edge(c, Leaf(4), NP(2, 2)),
             Edge(c, NP(1, 3), NP(2, 3))])
    bracketB = CircularGraph(CircularWord([c, c, c, c]), [circular_node([c, c, c]), circular_node([c, c, c])],
        Edge[Edge(c, Leaf(1), NP(1, 2)), Edge(c, Leaf(2), NP(1, 3)),
             Edge(c, Leaf(3), NP(2, 2)), Edge(c, Leaf(4), NP(2, 3)),
             Edge(c, NP(1, 1), NP(2, 1))])
    @test is_wired(bracketA) && is_wired(bracketB)
    # `bracketA == bracketB`: the cyclic neighbour sequences are identical up to
    # rotation (n1: `L1→L2→n2` against `n2→L1→L2`, n2 likewise), the leaf numbers
    # agree, and both nodes are `[1,1,1]`, i.e. fully rotation-symmetric. They are
    # the SAME diagram.
    #
    # A and B are NOT "the two bracketings" of 1111 — in
    # both, leaves 1,2 hang on one trivalent and 3,4 on the other; it is the SAME
    # bracketing, numbered two different ways. The genuinely other bracketing groups
    # (2,3) and (4,1); it stays distinct, because leaf numbers enter the key absolutely
    # and there is NO minimisation over leaf rotations (the left marker, see
    # `circular_canonical_key`). That is exactly what the test right below is for.
    @test bracketA == bracketB
    # the GENUINELY other bracketing (leaves 2,3 | 4,1) stays distinguishable —
    # otherwise too many diagrams would merge
    bracketC = CircularGraph(CircularWord([c, c, c, c]), [circular_node([c, c, c]), circular_node([c, c, c])],
        Edge[Edge(c, Leaf(2), NP(1, 1)), Edge(c, Leaf(3), NP(1, 2)),
             Edge(c, Leaf(4), NP(2, 1)), Edge(c, Leaf(1), NP(2, 2)),
             Edge(c, NP(1, 3), NP(2, 3))])
    @test is_wired(bracketC)
    @test bracketA != bracketC
    @test bracketB != bracketC
    comboA, histA = reduce_circular(circular_decorated(bracketA))
    comboB, histB = reduce_circular(circular_decorated(bracketB))
    @test histA == [:merge] && histB == [:merge]
    ((termA, _),) = pairs_of(comboA)
    ((termB, _),) = pairs_of(comboB)
    @test circular_canonical_key(termA) == circular_canonical_key(termB)   # SAME result
    @test arms(termA.graph.nodes[1]) == [1, 1, 1, 1]
end

@testset "reduce_circular_full — all_matches confluence and the shared memo" begin
    NP = DiagrammaticHecke.NodePort
    c = 1
    # the same bracketing-A fixture: several rules could in principle compete for the
    # first match (here only _fr_merge); all_matches=true additionally checks that
    # EVERY applicable rule gives the same result.
    bracketA = CircularGraph(CircularWord([c, c, c, c]), [circular_node([c, c, c]), circular_node([c, c, c])],
        Edge[Edge(c, Leaf(1), NP(1, 1)), Edge(c, Leaf(2), NP(1, 2)),
             Edge(c, Leaf(3), NP(2, 1)), Edge(c, Leaf(4), NP(2, 2)),
             Edge(c, NP(1, 3), NP(2, 3))])
    fdA = circular_decorated(bracketA)
    resA = reduce_circular_full(fdA; all_matches = true)          # does NOT throw ⇒ confluent
    @test length(resA) == 1

    # memoisation checked indirectly via equality of results: two calls with a SHARED
    # memo must give the same result as one call with a
    # fresh memo (no access to memo internals, only comparing results).
    memo = Dict{Any,CircularComboR}()
    r1 = reduce_circular_full(fdA; memo = memo)
    r2 = reduce_circular_full(fdA; memo = memo)                   # cache hit
    r_fresh = reduce_circular_full(fdA)
    @test r1 == r2
    @test r1 == r_fresh
end

@testset "_lift — label transfer for CircularGraph rules" begin
    NP = DiagrammaticHecke.NodePort

    # Positive: _fr_bead fires on a bead [1,1] (UNINVOLVED) while a FAR AWAY cell (the
    # outer cell of the distant trivalent) carries a non-trivial α₃ label — after the
    # bead collapse the label must land unchanged on the corresponding cell of the
    # result.
    nd_ctx = CircularNode(:mixed, [1, 1])
    nd_far = circular_node([1, 1, 1])
    g_ctx = CircularGraph(CircularWord([1, 1, 1]), [nd_ctx, nd_far],
        # ORIENTATION: arms cw against the ccw leaf ring — at nd_far the leaves run
        # backwards through the slots.
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, NP(1, 2), NP(2, 1)),
             Edge(1, Leaf(3), NP(2, 2)), Edge(1, Leaf(2), NP(2, 3))])
    @test face_count(g_ctx) == 3
    fd_ctx = circular_decorated(g_ctx, [one(SoergelPoly), alpha(3), one(SoergelPoly)])
    combo, hist = reduce_circular(fd_ctx)
    @test hist == [:bead]
    @test length(combo) == 1
    ((d, coeff),) = pairs_of(combo)
    @test coeff == one(SoergelPoly)
    @test length(d.graph.nodes) == 1                    # only the distant trivalent remains
    @test any(f -> f == alpha(3), d.region_labels)         # the label SURVIVES

    # A barbell whose only cell, STRUCTURALLY having no boundary gap (every barbell is
    # an isolated two-node component away from the boundary, see the test above),
    # carries a non-trivial label by hand — that cell disappears with the barbell, so
    # the label must be transferred somewhere.
    #
    # The rule: a scalar may only be pulled out of the diagram if it sits in the
    # region at distance 0. It puts the polynomial IN PLACE into the outer region
    # instead: `outer_label = α₁ (from the barbell) · α₂ (the
    # label)`, and the coefficient stays 1. The value is the same, only the place it
    # is stored differs; `circular_extract_scalars` only pulls it out at distance 0
    # (here there is no marker at all — a wordless boundary).
    nd1 = circular_node([1]); nd2 = circular_node([1])
    g_bad = CircularGraph(CircularWord(Int[]), [nd1, nd2], Edge[Edge(1, NP(1, 1), NP(2, 1))])
    @test face_count(g_bad) == 1
    fd_bad = circular_decorated(g_bad, [alpha(2)])
    combo_bad, hist_bad = reduce_circular(fd_bad)
    @test hist_bad == [:barbell]
    @test length(combo_bad) == 1
    ((d_bad, coeff_bad),) = pairs_of(combo_bad)
    @test coeff_bad == one(SoergelPoly)
    @test d_bad.outer_label == alpha(1) * alpha(2)
    @test isempty(d_bad.graph.nodes)
    # without a label: the pure barbell factor — likewise outside
    ((d_plain, c_plain),) = pairs_of(first(reduce_circular(circular_decorated(g_bad))))
    @test c_plain == one(SoergelPoly)
    @test d_plain.outer_label == alpha(1)
    # an already present outer label is not lost, it multiplies up (the target slot is
    # seeded by `_fr_barbell`).
    ((d_seed, _),) = pairs_of(first(reduce_circular(
        circular_decorated(g_bad, [alpha(2)], alpha(3)))))
    @test d_seed.outer_label == alpha(1) * alpha(2) * alpha(3)

    # The hardness sits one level deeper: `outer` is opt-in, and without it the
    # transfer throws (as all other call sites do).
    g_bad2 = DiagrammaticHecke._circular_delete_nodes(g_bad, Set([1, 2]), Edge[])
    @test_throws ErrorException DiagrammaticHecke._circular_transfer_labels(
        g_bad, g_bad2, [alpha(2)])
    # with `outer` it lands there, and the slot is seeded by the caller
    outer = Ref(alpha(1))
    neu = DiagrammaticHecke._circular_transfer_labels(g_bad, g_bad2, [alpha(2)]; outer = outer)
    @test outer[] == alpha(1) * alpha(2)
    @test all(isone, neu)
end

# ---- the outer label on the circular side: storage, key, distance ----------------
#
# Mirror of the plain testset "outer_label: storage, key, rendering"
# (test/decorated_rules.jl) for `CircularDecorated`/`CircularDecoratedMorphism`.
@testset "outer_label (circular): storage, key, the distance-0 rule" begin
    NP = DiagrammaticHecke.NodePort
    R  = SoergelPoly

    g = CircularGraph(CircularWord(Int[]), [circular_node([1]), circular_node([1])],
                 Edge[Edge(1, NP(1, 1), NP(2, 1))])
    labels = fill(one(R), region_count(g))

    # backward-compatible constructor: third argument optional, default 1
    D  = CircularDecorated(g, labels)
    Da = CircularDecorated(g, labels, alpha(1))
    @test outer_label(D) == one(R)
    @test outer_label(Da) == alpha(1)
    @test circular_decorated(g, labels, alpha(1)) == Da

    # `==`/`hash`/`circular_canonical_key` take the outer label into account
    @test D != Da
    @test hash(D) != hash(Da)
    @test circular_canonical_key(D) != circular_canonical_key(Da)
    @test circular_canonical_key(Da) == circular_canonical_key(CircularDecorated(g, labels, alpha(1)))

    # two terms differing only outside do NOT merge in a combo
    c = CircularComboR(D) + CircularComboR(Da)
    @test length(c) == 2

    # a 0 outside annihilates the term just as a 0 inside does
    @test isempty(CircularComboR(CircularDecorated(g, labels, zero(R))))

    # `show` names the outer region explicitly
    @test occursin("[outer]", sprint(show, Da))

    # ---- the distance-0 rule in `circular_extract_scalars` -----------------------
    # The boundary circle is a wall, and the outer region sits right in front of the
    # marking gap — so it touches the marking and its distance is 0 whenever there IS
    # a marking. A genuine polynomial outside therefore becomes a coefficient, exactly
    # like a constant.
    m = circular_morphism(light_leaf([1], [1]))
    lab = fill(one(R), region_count(m.graph))
    @test DiagrammaticHecke.circular_outer_region_distance(m) == 0
    f1, fdm1 = DiagrammaticHecke.circular_extract_scalars(
        CircularDecoratedMorphism(m, lab, alpha(1)))
    @test f1 == alpha(1)
    @test fdm1.outer_label == one(R)
    f2, fdm2 = DiagrammaticHecke.circular_extract_scalars(
        CircularDecoratedMorphism(m, lab, 3 * one(R)))
    @test f2 == 3 * one(R)
    @test fdm2.outer_label == one(R)
    # An UNMARKED `CircularDecorated` has no distance at all (`-1`), so the label
    # stays put there — and the way through the driver does not lose it.
    @test DiagrammaticHecke._circular_decorated(
        CircularDecoratedMorphism(m, lab, alpha(1))).outer_label == alpha(1)
end

@testset "reduce_circular — negative cases" begin
    NP = DiagrammaticHecke.NodePort

    # an irreducible diagram stays: a bare dot matches no rule.
    fd_dot = circular_decorated(circular(dot(1)))
    combo, hist = reduce_circular(fd_dot)
    @test isempty(hist)
    @test length(combo) == 1
    @test only(pairs_of(combo))[1] == fd_dot

    # a :braid self-loop ⇒ 0: the driver fires :needle and gives the empty sum.
    fgb = circular(braid(1, 2; m = 3))
    ndb = fgb.nodes[1]
    g_braid_self = CircularGraph(CircularWord(Int[]), [ndb, circular_node([1, 1, 1])],
        Edge[Edge(1, NP(1, 1), NP(1, 4)), Edge(2, NP(1, 2), NP(1, 5)),
             Edge(1, NP(1, 3), NP(2, 1)), Edge(1, NP(1, 6), NP(2, 2))])
    combo_bs, hist_bs = reduce_circular(circular_decorated(g_braid_self))
    @test :needle in hist_bs
    @test isempty(combo_bs)                             # the empty sum = 0
end

