# test/circularmisc.jl — CircularDecorated fusion (D1–D3), circular SVG with cell
# labels, canonical-key leaf rotation, tree_word labels, and the distance-word
# criterion.

@testset "CircularDecorated + fusion (step 3, parts D1-D3)" begin
    NP = DiagrammaticHecke.NodePort
    # D1: the length of region_labels must match region_count.
    nd = circular_node([1, 1, 1])
    g = CircularGraph(CircularWord([1, 1, 1]), [nd],
        # ORIENTATION (src/diagram/Faces.jl): boundary ccw, arms cw — the slots run
        # BACKWARDS along the leaf ring (slot 1 holds leaf 1 fixed).
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(3), NP(1, 2)), Edge(1, Leaf(2), NP(1, 3))])
    @test_throws ArgumentError CircularDecorated(g, [alpha(1)])   # only 1 instead of 3 labels
    fd_default = circular_decorated(g)
    @test all(isone, fd_default.region_labels)
    fd = circular_decorated(g, [alpha(1), one(SoergelPoly), one(SoergelPoly)])
    @test region_label(fd, 1) == alpha(1)
    @test face_count(fd) == 3
    @test length(fd) == length(g)

    # canonical_key includes the labels: same structure, different labels ⇒ different
    # key.
    fd_other = circular_decorated(g, [alpha(2), one(SoergelPoly), one(SoergelPoly)])
    @test circular_canonical_key(fd) != circular_canonical_key(fd_other)
    @test fd != fd_other
    @test fd == circular_decorated(g, [alpha(1), one(SoergelPoly), one(SoergelPoly)])

    # D3 fusion: s_1(α_1) = -α_1, Δ_1(α_1) = 2.
    r = DiagrammaticHecke._fr_fusion(fd, 1)
    @test r !== nothing
    @test length(r) == 2
    for (dd, coeff) in DiagrammaticHecke.pairs_of(r)
        @test coeff == one(SoergelPoly)
        if length(dd.graph.nodes) == 1
            # term 1: the strand stays, the label moves (act(1, α_1) = -α_1).
            @test any(f -> f == -alpha(1), dd.region_labels)
        else
            # term 2: broken open, two fresh dots, demazure(1, α_1) = 2.
            @test length(dd.graph.nodes) == 3
            @test any(f -> f == 2 * one(SoergelPoly), dd.region_labels)
        end
    end

    # nothing if the edge does not separate two distinct cells (degree-1 dot, its own
    # capping edge).
    @test DiagrammaticHecke._fr_fusion(circular_decorated(circular(dot(1))), 1) === nothing
end

@testset "circular SVG with cell labels (step 3, part D4)" begin
    NP = DiagrammaticHecke.NodePort
    nd = circular_node([1, 1, 1])
    g = CircularGraph(CircularWord([1, 1, 1]), [nd],
        # ORIENTATION (src/diagram/Faces.jl): boundary ccw, arms cw — the slots run
        # BACKWARDS along the leaf ring (slot 1 holds leaf 1 fixed).
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(3), NP(1, 2)), Edge(1, Leaf(2), NP(1, 3))])
    fd = circular_decorated(g, [alpha(1), one(SoergelPoly), one(SoergelPoly)])
    svg = circular_decorated_svg(fd)
    @test occursin("<svg", svg)
    @test occursin("α", svg)   # the α_1 label appears in the SVG

    # the label 1 is hidden: a completely unlabelled diagram has no <text> cell label
    # (the boundary-leaf captions from circular_tutte_svg do not count, they are digits,
    # not <text> cell labels — checking that NO 'α' appears suffices here).
    fd_triv = circular_decorated(g)
    svg_triv = circular_decorated_svg(fd_triv)
    @test !occursin("α", svg_triv)

    @test showable("image/svg+xml", fd)
    @test showable("image/png", fd)
    @test occursin("<svg", repr("image/svg+xml", fd))
    @test length(repr("image/png", fd)) > 100

    # show_circular_combo: must not throw on an empty combination.
    @test DiagrammaticHecke.show_circular_combo(CircularComboR()) === nothing
end

@testset "show_circular_combo(c, m) — the mark_between convenience overload" begin
    # Regression: result terms of reduce_to_circular_leave (bare CircularDecorated, i.e. "cuts
    # forgotten") drew NO boundary marker without this overload — noticed when
    # comparing against the correctly marked DL diagrams. Test: the
    # overload gives a DIFFERENT SVG from the default (mark_between=nothing), and
    # pulls `_circular_left_mark(m)` correctly out of `m`.
    ll = light_leaf([1, 1], [0, 1])
    hf = hflip(ll)
    fm = circular_morphism(hf)
    fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
    combo = reduce_to_circular_leave(fdm)
    @test !isempty(combo)
    for (d, _) in pairs_of(combo)
        svg_default = circular_decorated_svg(d)
        svg_marked = circular_decorated_svg(d; mark_between = DiagrammaticHecke._circular_left_mark(fm))
        @test svg_default != svg_marked
    end
    # show_circular_combo(c, m) must use the same mark_between derivation as checked by
    # hand above — no `display(...)` call here (no MIME backend is registered in the
    # headless test run, as with the other `display_*` calls in this file: only
    # `repr`/`showable` are checked, never `display` directly).
    @test DiagrammaticHecke._circular_left_mark(fm) !== nothing
end

@testset "the GREY boundary marker in the decorated circular path" begin
    # The term pictures of `show_circular_combo(c, m)` must not carry
    # the black but NO grey marker, although `display(circular_morphism(m))` has both. The
    # cause was not `show_circular_combo` itself — it draws via `CircularDecoratedMorphism`,
    # whose `circular_decorated_svg` did not pass `_circular_right_mark` through (only
    # `_circular_left_mark`). The marker is not mere decoration: it enters
    # `circular_tutte_positions` through `_add_nudge!`.
    m = vflip(circular_morphism(light_leaf([2, 3, 2, 3], [0, 1, 1, 1])))
    @test DiagrammaticHecke._circular_right_mark(m) !== nothing
    grey(s) = length(collect(eachmatch(r"fill=\"#999\"", s)))
    fdm = CircularDecoratedMorphism(m, fill(one(SoergelPoly), region_count(m.graph)))
    # the same number of grey markers as on the direct path `circular_tutte_svg(m)`.
    @test grey(circular_decorated_svg(fdm)) == grey(circular_tutte_svg(m)) == 1
end


@testset "the crossing shape is REACHABLE by merging" begin
    # That a `[1,3,1,3]` node is structurally indistinguishable from a converted
    # `braid(1,3;m=2)` is in the testset above. Here: REACHABILITY — the alternating
    # crossing shape does ARISE by merging, so the coverage of R8/R12 is not
    # hypothetical.
    cross = DiagrammaticHecke._min_rotation([1, 3, 1, 3])

    # all 4-arm shapes a splice of two 1/3 nodes can give — and only the VALID ones:
    # a splice whose result `circular_node` rejects throws. "Reachable" therefore
    # means "reachable by a LEGAL merge", which is what the `circular_node` filter
    # below expresses.
    reach = Set{Vector{Int}}()
    cands = ([1, 3, 1], [3, 1, 3], [1, 1, 3], [3, 3, 1], [1, 3], [1, 1], [3, 3])
    for ax in cands, ay in cands
        length(ax) + length(ay) - 2 == 4 || continue
        for si in 1:length(ax), sj in 1:length(ay)
            ax[si] == ay[sj] || continue
            na = vcat([ax[mod1(si + k, length(ax))] for k in 1:(length(ax) - 1)],
                      [ay[mod1(sj + k, length(ay))] for k in 1:(length(ay) - 1)])
            try
                circular_node(na)
            catch
                continue                       # not a legal merge — does not count
            end
            push!(reach, DiagrammaticHecke._min_rotation(na))
        end
    end
    @test cross in reach                       # the crossing shape IS reachable
    @test DiagrammaticHecke._min_rotation([1, 1, 3, 3]) in reach   # control: not everything is equal

    # BUT: reachability depends on the SLOT CHOICE — the same two nodes
    # `[1,3,1]`+`[1,3,1]` give the crossing `[3,1,3,1]` across the 1-edge at slot
    # 1/slot 1, but `[3,1,1,3]` (~`[1,1,3,3]`) across slot 1/slot 3, i.e. precisely
    # NOT the crossing. `[1,3,1]` itself is invalid as a STANDALONE node under the
    # "every colour at least twice" invariant (circular_node would throw) — exactly the
    # intermediate state `_unchecked_circularnode` exists for: built here to pass into a
    # valid shape immediately afterwards (via merge_at_edge). The pair
    # `[1,3,1]`+`[3,1,3]` → `[1,3,3,3]` gives a valid node for NO slot choice and
    # throws, so it does not serve as a contrast case.
    NP = DiagrammaticHecke.NodePort
    ax = [1, 3, 1]
    function mk(sj)
        nd = DiagrammaticHecke._unchecked_circularnode(:mixed, ax)
        es = Edge[Edge(1, NP(1, 1), NP(2, sj))]      # the 1-edge, slot 1 onto slot sj
        push!(es, Edge(3, Leaf(1), NP(1, 2)))
        push!(es, Edge(1, Leaf(2), NP(1, 3)))
        w = [3, 1]
        for k in 1:3
            k == sj && continue
            push!(es, Edge(ax[k], Leaf(length(w) + 1), NP(2, k)))
            push!(w, ax[k])
        end
        return CircularGraph(CircularWord(w), [nd, nd], es)
    end
    m_cross = DiagrammaticHecke.merge_at_edge(mk(1), 1)
    @test arms(m_cross.nodes[1]) == [3, 1, 3, 1]              # THE crossing
    @test DiagrammaticHecke._min_rotation(arms(m_cross.nodes[1])) == cross
    m_flat = DiagrammaticHecke.merge_at_edge(mk(3), 1)
    @test arms(m_flat.nodes[1]) == [3, 1, 1, 3]               # the same nodes, different slots
    @test DiagrammaticHecke._min_rotation(arms(m_flat.nodes[1])) != cross

    # and despite the same arm sequence the crossing is not the same object as a
    # converted braid(1,3; m=2).
    @test m_flat.nodes[1] != circular(braid(1, 3; m = 2)).nodes[1]
end

@testset "circular_canonical_key and leaf rotation" begin
    # A THIRD freedom: `g.word` is a CircularWord, defined only up to rotation — but
    # `Leaf(k)` is an ABSOLUTE position. `[2,2,2]` (rotation-symmetric) at a
    # trivalent: the same wiring, once started with leaf 1 at slot 1, once with leaf 1
    # at slot 2 — the same diagram apart from where the leaf count starts.
    NP = DiagrammaticHecke.NodePort
    nd = circular_node([2, 2, 2])
    gA = CircularGraph(CircularWord([2, 2, 2]), [nd],
        # ORIENTATION (src/diagram/Faces.jl): boundary ccw, arms cw — the slots run
        # BACKWARDS along the leaf ring (slot 1 holds leaf 1 fixed).
        Edge[Edge(2, Leaf(1), NP(1, 1)), Edge(2, Leaf(3), NP(1, 2)), Edge(2, Leaf(2), NP(1, 3))])
    gB = CircularGraph(CircularWord([2, 2, 2]), [nd],
        Edge[Edge(2, Leaf(1), NP(1, 2)), Edge(2, Leaf(2), NP(1, 3)), Edge(2, Leaf(3), NP(1, 1))])
    @test is_wired(gA) && is_wired(gB)
    # The LEFT MARKER is relevant to whether two diagrams are equal
    # — practically, the word is fixed (no rotation in the circular word), and only
    # then can diagrams be compared. `circular_canonical_key` therefore does NOT
    # minimise over `_circular_leaf_rotations`; gA and gB differ exactly by where the
    # leaf count starts and are therefore DIFFERENT.

    @test circular_canonical_key(gA) != circular_canonical_key(gB)
    @test gA != gB

    # NEGATIVE CONTROL against over-collapse: `[1,1,3,3]` is NOT rotation-symmetric
    # (the only valid rotation is r=0, checked) — an implementation minimising blindly
    # over ALL r instead of only the word-valid ones would wrongly identify two
    # different diagrams here.
    # gC arises from gA' by EXACTLY the leaf renaming Leaf(k) -> Leaf(mod1(k-1,4)),
    # i.e. the rotation r=1 — which is NO symmetry for this word. So gA' and gC are
    # genuinely different diagrams and MUST keep different keys.
    @test DiagrammaticHecke._circular_leaf_rotations(CircularWord([1, 1, 3, 3])) == [0]
    nd2 = circular_node([1, 1, 3, 3])
    gAp = CircularGraph(CircularWord([1, 1, 3, 3]), [nd2],
        Edge[Edge(1, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 2)),
             Edge(3, Leaf(3), NP(1, 3)), Edge(3, Leaf(4), NP(1, 4))])
    gC = CircularGraph(CircularWord([1, 1, 3, 3]), [nd2],
        Edge[Edge(1, Leaf(4), NP(1, 1)), Edge(1, Leaf(1), NP(1, 2)),
             Edge(3, Leaf(2), NP(1, 3)), Edge(3, Leaf(3), NP(1, 4))])
    @test is_wired(gAp) && is_wired(gC)
    @test circular_canonical_key(gAp) != circular_canonical_key(gC)
    @test gAp != gC
end

# ---- tree_word labels --------------------------------------------------------
# Alongside the region distance (circular_region_distances, the version with node
# adjacencies) each region gets the word of the crossed edge colours —
# `circular_region_tree_words` (src/circular/CircularRegion.jl). Invariant: word
# length == distance. CONJECTURE: on end terms of `reduce_to_circular_leave` (normal
# form) all tree_words are reduced, so a non-reduced tree_word indicates "not in
# normal form".
@testset "tree_word labels per region (circular_region_tree_words)" begin
    # the invariant on three Zamo fixtures: length(word) == region distance,
    # unreachable (-1) <-> nothing; the start region carries the empty word.
    for m0 in (zamo_lhs(), zamo_rhs(), Zamo(1, 7))
        fm = circular_morphism(m0)
        d = circular_region_distances(fm)
        w = circular_region_tree_words(fm)
        @test length(w) == length(d)
        @test all(i -> (d[i] == -1) == (w[i] === nothing), eachindex(d))
        @test all(i -> w[i] === nothing || length(w[i]) == d[i], eachindex(d))
        @test any(x -> x == Int[], w)                  # the start region
    end
    # the conjecture on the end terms of Zamo(1,7): everything is reduced.
    fm = circular_morphism(Zamo(1, 7))
    fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
    c = reduce_to_circular_leave(fdm; maxrounds = 100, node_growth_limit = 15)
    @test !isempty(collect(pairs_of(c)))
    for (dd, _) in pairs_of(c)
        m3 = CircularMorphismGraph(dd.graph, fm.cut1, fm.cut2)
        ws = circular_region_tree_words(m3)
        @test all(x -> x === nothing || is_reduced(x), ws)
    end
end

# ---- the distance-word criterion ---------------------------------------------
# The rule: for dots we look at the distance word; if it can end
# on the dot's colour we can reduce further. […] That's exactly what the 1212
# example was built for. The word comes from circular_region_words (the edge BFS,
# distance 0 -> outside) and is only rotated — NOT searched
# backwards from the dot.
@testset "the distance-word criterion (circular_dot_reducible)" begin
    # the braid class does not shorten: the length stays.
    @test circular_braid_class([1, 2, 1]) == Set([[1, 2, 1], [2, 1, 2]])
    @test circular_braid_class([1, 3])    == Set([[1, 3], [3, 1]])
    @test circular_braid_class([2])       == Set([[2]])
    @test all(w -> length(w) == 4, circular_braid_class([1, 2, 1, 2]))

    # The class is `reduced_words`, so no representative may carry a double letter.
    # Applying `[2,1,2]->[1,2,1]` at two overlapping positions of `21212` would give
    # `21121`/`12112`/`22122`, which are excluded.
    for w in ([1, 2, 1], [1, 2, 1, 2], [2, 1, 3, 2], [1, 2, 3, 1, 2])
        for v in circular_braid_class(w)
            @test all(i -> v[i] != v[i + 1], 1:(length(v) - 1))
            @test same_element(v, w)
        end
    end
    # a NON-reduced word has no braid class here: the singleton, no witnesses
    # (non-reduced region words belong to the separate region-word criterion, E2).
    @test circular_braid_class([2, 2]) == Set([[2, 2]])
    @test circular_braid_class([2, 1, 2, 1, 2]) == Set([[2, 1, 2, 1, 2]])

    # the 1212 example: 1212 at the bottom, 121 at
    # the top, a 2-dot on the right. Its tree_word is 121, the class {121, 212} — and
    # 212 ends in 2 = the dot colour, hence REDUCIBLE.
    NP, L, E = NodePort, Leaf, Edge
    nodes = CircularNode[circular_node([1,2,1,2,1,2]), circular_node([1,2,1,2,1,2]), circular_node([2])]
    edges = E[E(1, L(2), NP(1,1)), E(2, L(3), NP(1,6)), E(1, L(4), NP(1,5)),
              E(1, L(6), NP(2,1)), E(2, L(7), NP(2,6)), E(1, L(1), NP(2,5)),
              E(2, L(5), NP(3,1)),
              E(2, NP(1,4), NP(2,2)), E(1, NP(1,3), NP(2,3)), E(2, NP(1,2), NP(2,4))]
    g = CircularGraph(CircularWord([1,2,1,2,1,2,1]), nodes, edges)
    @test euler(g) == 2
    @test isempty(check_wiring(g))
    m = CircularMorphismGraph(g, 1, 5)

    @test circular_distance_word(m, 3) == [1, 2, 1]
    @test circular_dot_reducible(m, 3)
    @test circular_non_normal_dots(m) == [3]

    # no dot -> error (node 1 has 6 arms).
    @test_throws ErrorException circular_distance_word(m, 1)

    # counter-check: the end terms of Zamo(1,7) are normal forms, so NO dot may
    # violate the criterion there (the same conjecture as for the tree_word).
    fm = circular_morphism(Zamo(1, 7))
    fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
    c = reduce_to_circular_leave(fdm; maxrounds = 100, node_growth_limit = 15)
    @test !isempty(collect(pairs_of(c)))
    for (dd, _) in pairs_of(c)
        m3 = CircularMorphismGraph(dd.graph, fm.cut1, fm.cut2)
        @test isempty(circular_non_normal_dots(m3))
    end
end

