# Focused checks for D4 in the circular world, on
# CircularGraph/CircularDecorated/CircularMorphismGraph.

@testset "region adjacency: leaf edges separate the gaps k-1 and k" begin
    # The dart->region assignment must not be shifted by a boundary segment (a wrong
    # gap label or a wrong search direction does that, see
    # src/circular/CircularRegion.jl).
    #
    # A criterion independent of the tracer: an edge at leaf k that does not hang on a
    # dot separates exactly the regions of gaps k-1 and k.
    function leaf_edges_separate_own_gaps(g)
        n = length(g.word)
        g2r = Dict{Int,Int}()
        for R in DiagrammaticHecke.regions(g), k in R.gaps
            g2r[k] = R.id
        end
        dr, t = DiagrammaticHecke._circular_dart_region(g)
        checked = 0
        for e in g.edges
            lp = e.a isa Leaf ? e.a : (e.b isa Leaf ? e.b : nothing)
            lp === nothing && continue
            np = e.a isa Leaf ? e.b : e.a
            (np isa NodePort && arm_count(g.nodes[np.node]) == 1) && continue
            haskey(t.port_dart, e.a) || continue
            d = t.port_dart[e.a]
            R1 = get(dr, d, 0); R2 = get(dr, t.darts[d].rev, 0)
            (R1 == 0 || R2 == 0) && continue
            checked += 1
            want = Set([get(g2r, mod1(lp.k - 1, n), 0), get(g2r, lp.k, 0)])
            Set([R1, R2]) == want || return (checked, false)
        end
        (checked, true)
    end

    for (w, e) in [([3, 2, 1, 3, 2, 3], [1, 1, 1, 1, 0, 0]),
                   ([3, 2, 1, 3, 2], [1, 0, 1, 0, 0]),
                   ([3, 2, 1], [1, 1, 0]),
                   ([1, 2], [1, 0])]
        g = circular_morphism(light_leaf(w, e)).graph
        checked, ok = leaf_edges_separate_own_gaps(g)
        @test checked > 0
        @test ok
    end

    # The concrete case: the two regions in question are at distance 1, not 2. What is
    # checked is the DISTRIBUTION of the distances, not their assignment to region
    # numbers — that follows the cell order and shifts with the orientation, and the
    # order of the regions does not matter.
    fm = circular_morphism(light_leaf([3, 2, 1, 3, 2, 3], [1, 1, 1, 1, 0, 0]))
    d = circular_region_distances(fm)
    @test sort(d) == sort([1, 2, 3, 4, 4, 3, 2, 1, 0, 3, 1, 2, 2, 3, 1])
    @test count(==(0), d) == 1        # exactly one region carries the marker
    @test maximum(d) == 4
    @test all(>=(0), d)               # no region is unreachable
    # and the distances are a genuine BFS layering: every level occurs before
    @test sort(unique(d)) == collect(0:maximum(d))
end

@testset "D4 on Circular — apply_circular_d4 / find_circular_d4_match" begin
    # dot_A = dot_B + dot_C - BC_decorated_trivalent (cyclically), on CircularGraph.
    function circular_dotmap(k::Int)
        others = filter(!=(k), (1, 2, 3))
        a, b = others
        CircularGraph(CircularWord([1, 1, 1]), [circular_node([1])],
            [Edge(1, Leaf(k), NodePort(1, 1)), Edge(1, Leaf(a), Leaf(b))])
    end
    function circular_trivmap()
        # ORIENTATION (src/diagram/Faces.jl): boundary ccw, arms cw — the slots run
        # BACKWARDS along the leaf ring (slot 1 holds leaf 1 fixed). Written forwards,
        # the embedding would not be planar.
        CircularGraph(CircularWord([1, 1, 1]), [circular_node([1, 1, 1])],
            [Edge(1, Leaf(1), NodePort(1, 1)), Edge(1, Leaf(3), NodePort(1, 2)),
             Edge(1, Leaf(2), NodePort(1, 3))])
    end
    # Which region carries the α₁ in the trivalent term, by dot position k. The
    # numbers are region/cell numbers and therefore depend on the numbering; with the
    # orientation correction k = 2 and k = 3 swap.
    pair_cell = Dict(1 => 3, 2 => 1, 3 => 2)

    # ---- apply_circular_d4 on the direct circular_dotmap(k) shapes, all three k ---------
    for k in 1:3
        others = filter(!=(k), (1, 2, 3))
        d = circular_decorated(circular_dotmap(k))
        res = apply_circular_d4(d, 1, 2)
        @test res !== nothing
        # 3 structurally different terms. `circular_canonical_key` does NOT minimise
        # over leaf rotations (the left marker is part of the identity), and the two
        # "+1" dot terms differ exactly there — the dot sits on one or the other of the
        # two free leaves — so they stay apart: twice +1 dot, −1 trivalent.
        @test length(pairs_of(res)) == 3
        dotterms = 0
        trivterms = 0
        for (dd, c) in pairs_of(res)
            if arm_count(dd.graph.nodes[1]) == 1
                dotterms += 1
                @test c == one(R)
                @test any(o -> circular_canonical_key(dd.graph) == circular_canonical_key(circular_dotmap(o)), others)
            else
                trivterms += 1
                @test c == -one(R)
                @test circular_canonical_key(dd.graph) == circular_canonical_key(circular_trivmap())
                # The labels are REGION-indexed while `pair_cell` is phrased in FACE
                # numbers, so it is translated via `R.cell` here.
                want = findfirst(R -> R.cell == pair_cell[k], DiagrammaticHecke.regions(dd.graph))
                nontrivial = findall(!isone, dd.region_labels)
                @test nontrivial == [want]
                @test dd.region_labels[want] == alpha(1)
            end
        end
        @test dotterms == 2
        @test trivterms == 1
    end

    # ---- negative matches ------------------------------------------------------
    d1 = circular_decorated(circular_dotmap(1))
    @test apply_circular_d4(d1, 1, 1) === nothing            # ei = own capping edge
    @test apply_circular_d4(d1, 99, 2) === nothing            # invalid dot_node
    @test apply_circular_d4(d1, 1, 99) === nothing            # invalid ei

    # ---- find_circular_d4_match: needs a marking, matches the base case --------
    ll = light_leaf([1, 1], [0, 1])
    hf = hflip(ll)
    fm = circular_morphism(hf)
    @test find_circular_d4_match(fm) == (1, 2)
end

@testset "reduce_to_circular_leave — base case (hflip of light_leaf)" begin
    ll = light_leaf([1, 1], [0, 1])
    hf = hflip(ll)
    fm = circular_morphism(hf)

    fdm = CircularDecoratedMorphism(fm, fill(one(R), region_count(fm.graph)))
    result = reduce_to_circular_leave(fdm)
    @test result isa DiagrammaticHecke.CircularComboR
    @test !isempty(result)
    # every surviving term must have NO further D4 match (fixed point reached)
    for (d, c) in pairs_of(result)
        m2 = CircularMorphismGraph(d.graph, fm.cut1, fm.cut2)
        @test find_circular_d4_match(m2) === nothing
    end
    # the α_1 factor from the decorated trivalent term must survive somewhere in
    # the result (as a coefficient, extracted from a marking-distance-0 cell by
    # circular_extract_scalars)
    @test any(c -> occursin("α", string(c)), last.(pairs_of(result)))
    # the two dot terms do not coincide (see the D4 testset above): 2 dot terms +
    # 1 trivalent term.
    @test length(pairs_of(result)) == 3
end

@testset "reduce_to_circular_leave — termination on a multi-dot example" begin
    dls = double_leaves([1, 2, 1], [1, 2, 1])
    r = only(filter(x -> x.e == [0, 0, 1] && x.f == [1, 0, 0], dls))
    @test dot_count(r.morphism.graph) == 3   # genuine multi-dot starting point

    fm = circular_morphism(r.morphism)
    fdm = CircularDecoratedMorphism(fm, fill(one(R), region_count(fm.graph)))
    result = reduce_to_circular_leave(fdm)
    @test result isa DiagrammaticHecke.CircularComboR
    @test !isempty(result)
    for (d, c) in pairs_of(result)
        m2 = CircularMorphismGraph(d.graph, fm.cut1, fm.cut2)
        @test find_circular_d4_match(m2) === nothing   # fixed point: no term has a match left
    end
end

@testset "reduce_to_circular_leave — negative case: no D4 match is a fixed point" begin
    # identity_morphism has no dots at all — no D4 match possible.
    im = identity_morphism()
    fm = circular_morphism(im)
    @test find_circular_d4_match(fm) === nothing

    fdm = CircularDecoratedMorphism(fm, fill(one(R), region_count(fm.graph)))
    result = reduce_to_circular_leave(fdm)
    @test length(pairs_of(result)) == 1
    (d, c) = only(pairs_of(result))
    @test c == one(R)
    @test circular_canonical_key(d.graph) == circular_canonical_key(fm.graph)
    @test all(isone, d.region_labels)
end

@testset "circular D4 exported from the top-level module" begin
    @test DiagrammaticHecke.apply_circular_d4 === apply_circular_d4
    @test DiagrammaticHecke.find_circular_d4_match === find_circular_d4_match
    @test DiagrammaticHecke.reduce_to_circular_leave === reduce_to_circular_leave
    @test DiagrammaticHecke.CircularDecoratedMorphism === CircularDecoratedMorphism
end

# ---- THE D4 NODE VARIANT -----------------------------------------------------

@testset "_circular_d4_node_far_arc — the arc between the arms, away from the dot" begin
    far = DiagrammaticHecke._circular_d4_node_far_arc
    # A 4-armed node whose sectors after slots 1..4 are the regions 10, 20, 30, 40.
    # Arms 1 and 3 split it into {10,20} and {30,40}.
    sec = [10, 20, 30, 40]
    @test far(sec, 4, 1, 3, 20) == [30, 40]   # dot in 20 ⇒ the far arc is the other one
    @test far(sec, 4, 1, 3, 30) == [10, 20]
    # dot region in BOTH arcs (can happen when a dot lets sectors merge) ⇒ not
    # unique, so no match.
    @test far([10, 20, 10, 40], 4, 1, 3, 10) === nothing
    # dot region in NEITHER of the two ⇒ likewise not unique.
    @test far(sec, 4, 1, 3, 99) === nothing
    # neighbouring arms: the far arc is everything except the one sector
    @test far(sec, 4, 1, 2, 10) == [20, 30, 40]
end

@testset "the D4 node variant — matcher, terms, label" begin
    NP = DiagrammaticHecke.NodePort
    # Fixture: the smallest such case (x = y = 1313, e = 1000, f = 0111, after the
    # first reduce_circular round). A 3-dot, a 1-dot, another
    # 3-dot and a [3,1,3,1,1] node: the dot region is bounded by 1-arms there, and the
    # two 3-arms (slots 1 and 3) lie on the other side — exactly the case for which the
    # edge variant is NOT responsible.
    g = CircularGraph(CircularWord([1, 1, 3, 1, 3, 3, 1, 3]),
        [circular_node([3]), circular_node([3]), circular_node([1]), circular_node([3, 1, 3, 1, 1])],
        Edge[Edge(3, Leaf(2), NP(1, 1)), Edge(1, Leaf(1), NP(4, 5)),
             Edge(1, Leaf(3), NP(4, 4)), Edge(3, Leaf(4), NP(2, 1)),
             Edge(1, Leaf(8), NP(3, 1)), Edge(3, Leaf(7), NP(4, 1)),
             Edge(1, Leaf(6), NP(4, 2)), Edge(3, NP(4, 3), Leaf(5))])
    @test is_wired(g)
    m = CircularMorphismGraph(g, 0, 4)

    # ---- matcher ----------------------------------------------------------
    # THE MATCHER DOES NOT FIRE HERE, and the reason is the distance:
    #
    # `circular_region_distances` additionally measures node adjacencies across
    # general-2/general-13 nodes (see the docstring in `src/circular/CircularRegion.jl`). Node 4
    # (`[3,1,3,1,1]`) is ITSELF a general-13 node, so its region 3 (sector after slot
    # 3, the dot region of dot 2) gains an extra colour-1 edge to region 2 (distance 0,
    # near the marker). `dist` is `[1,0,1,1,2]`, so region 3 sits at distance 1. A
    # hit needs `dist[R_between] < dist[R]` with `R_between` = region 4 (distance 1)
    # and `R` = region 3; with `dist[R] == dist[R_between] == 1` the condition fails
    # and NEITHER the node NOR the edge variant fires here (both checked below).
    #
    # The same computation resolves the 2-cycle in which `reduce_to_circular_leave`
    # does not otherwise terminate at DL 132 of the pair `212321/231231`: the cycle
    # node (`[1,3,1,3,3]`, 5 arms) is structurally the same general-13 type as this
    # test node (`[3,1,3,1,1]`, also 5 arms).
    #
    # So what is tested here is THAT both matcher variants stay empty for this
    # fixture; the `apply_circular_d4_node` term structure below is checked with the
    # values (`dot_node=2, v=4, p=1, q=3, insert_after=3`) — it does not depend on the
    # matcher but on `apply_circular_d4_node` itself.
    dist = circular_region_distances(m)
    @test dist == [1, 0, 1, 1, 2]
    @test find_circular_d4_node_match(m) === nothing
    @test find_circular_d4_match(m) === nothing
    @test find_circular_d4_any(m) === nothing

    dot_node, v, p, q, insert_after = 2, 4, 1, 3, 3
    @test arm_colour(g.nodes[v], p) == arm_colour(g.nodes[dot_node], 1)  # dot colour
    @test arm_colour(g.nodes[v], q) == arm_colour(g.nodes[dot_node], 1)
    # the inserted arm points into the dot region (sector 3, slot insert_after=3)
    @test circular_node_sector_regions(g)[v][insert_after] == circular_region_of_dot(g, dot_node)

    # ---- application ----------------------------------------------------------
    fd = CircularDecorated(g, fill(one(R), region_count(g)))
    res = apply_circular_d4_node(fd, dot_node, v, p, q, insert_after)
    @test res !== nothing
    terms = collect(pairs_of(res))
    @test length(terms) == 3
    @test count(((_, c),) -> c == one(R), terms) == 2
    @test count(((_, c),) -> c == -one(R), terms) == 1

    s = arm_colour(g.nodes[dot_node], 1)
    triv = only([t for (t, c) in terms if c == -one(R)])
    dots = [t for (t, c) in terms if c == one(R)]
    @test length(dots) == 2

    # trivalent term: merely connected — the node gains ONE arm and the old dot is
    # gone (so one dot fewer than before).
    ndots(t) = count(n -> arm_count(n) == 1, t.graph.nodes)
    @test ndots(triv) == ndots(fd) - 1
    @test maximum(arm_count, triv.graph.nodes) == arm_count(g.nodes[v]) + 1
    @test is_wired(triv.graph)
    # alpha_s sits in EXACTLY ONE region
    marked = [i for (i, l) in enumerate(triv.region_labels) if !isone(l)]
    @test length(marked) == 1
    @test triv.region_labels[only(marked)] == alpha(s)

    # dot terms: one arm goes and one arrives ⇒ the arm count is unchanged; the old
    # dot is gone and a fresh one sits at the far end.
    for t in dots
        @test ndots(t) == ndots(fd)
        @test maximum(arm_count, t.graph.nodes) == arm_count(g.nodes[v])
        @test is_wired(t.graph)
        @test all(isone, t.region_labels)      # no alpha in the dot terms
    end
    # the two dot terms are DIFFERENT (p goes in one, q in the other)
    @test circular_canonical_key(dots[1].graph) != circular_canonical_key(dots[2].graph)

    # ---- negative cases -------------------------------------------------------
    @test apply_circular_d4_node(fd, dot_node, dot_node, 1, 1, 1) === nothing  # v == Dot
    @test apply_circular_d4_node(fd, v, v, p, q, insert_after) === nothing     # no dot
    # p/q must carry the DOT colour: slot 2 of [3,1,3,1,1] is colour 1
    @test apply_circular_d4_node(fd, dot_node, v, 2, q, insert_after) === nothing
end

@testset "the D4 node variant is exported" begin
    @test DiagrammaticHecke.apply_circular_d4_node === apply_circular_d4_node
    @test DiagrammaticHecke.find_circular_d4_node_match === find_circular_d4_node_match
    @test DiagrammaticHecke.find_circular_d4_any === find_circular_d4_any
    @test DiagrammaticHecke.circular_node_sector_regions === circular_node_sector_regions
end

@testset "D4 chord case — ei is a chord of the dot cell" begin
    # The D4 edge ei has the same cell on both sides (dot_cell == bc_cell_old), so the
    # chains pair up as a 2-cycle {dse <-> X} plus a fixed point {Y -> Y} instead of a
    # 3-cycle. Planar-wise that is fine: the new trivalent node sits in the pocket
    # containing the dot and splits it into 2 regions, so
    # `_circular_d4_cyclic_ends` must accept it and not complain that the chain
    # concatenation does not close.
    dls = double_leaves([3, 2, 3, 1, 2, 3], [3, 2, 3, 1, 2, 3])
    r = only(filter(x -> x.e == [1, 1, 0, 1, 0, 0] && x.f == [1, 1, 1, 1, 0, 1], dls))
    fm = circular_morphism(r.morphism)
    fdm = CircularDecoratedMorphism(fm, fill(one(R), region_count(fm.graph)))
    result = reduce_to_circular_leave(fdm)      # round 2 is the delicate one
    @test result isa DiagrammaticHecke.CircularComboR
    @test !isempty(result)
    for (d, c) in pairs_of(result)
        # every result graph must have a consistent cell structure (planarity of the
        # new wiring) and no D4 match left open
        @test is_wired(d.graph)
        m2 = CircularMorphismGraph(d.graph, fm.cut1, fm.cut2)
        @test find_circular_d4_any(m2) === nothing
    end
end

@testset "a label in a gapless inner region — the FUSION resolves it" begin
    # Inserting the 2parallel 4-node at a same-coloured edge pair in the 12-arm
    # cluster puts the α_s/2 into a BIGON between the new 4-node and a trivalent. That
    # region has no boundary gap; its EDGE neighbours all have the same distance (the
    # fusion is stuck), and the only sector of smaller distance is the opposite sector
    # of the same 4-node — also gapless.
    #
    # A gapless region like this is NOT a dead end for the fusion. It looks like one
    # under the V0 DISTANCE, which pushes the region down from 4 to 3 so that its two
    # edge neighbours lie at the SAME distance and there is no `d − 1` neighbour.
    # Under the edge distance, on which the fusion computes, the region has distance 4
    # and the neighbours 3, so the fusion fires.
    Leaf, NodePort, Edge = DiagrammaticHecke.Leaf, DiagrammaticHecke.NodePort, DiagrammaticHecke.Edge
    n, s, t = 3, 1, 2
    nodes = CircularNode[]
    for _ in 1:n; push!(nodes, circular_node([s, t, s, t, s, t])); end
    for _ in 1:n; push!(nodes, circular_node([t, t, t])); end
    push!(nodes, circular_node(fill(s, n)))
    edges = Edge[]
    for i in 1:n
        a = 4 * (i - 1) + 1
        push!(edges, Edge(s, Leaf(a + 2), NodePort(i, 1)))
        push!(edges, Edge(t, Leaf(a + 1), NodePort(i, 2)))
        push!(edges, Edge(s, Leaf(a),     NodePort(i, 3)))
        push!(edges, Edge(t, NodePort(i, 4), NodePort(n + mod1(i - 1, n), 2)))
        push!(edges, Edge(t, NodePort(i, 6), NodePort(n + i, 1)))
        push!(edges, Edge(t, Leaf(a + 3), NodePort(n + i, 3)))
    end
    for i in 1:n                     # the s-spokes LAST (e19..e21)
        push!(edges, Edge(s, NodePort(i, 5), NodePort(2n + 1, n + 1 - i)))
    end
    g12 = CircularGraph(CircularWord(repeat([s, t, s, t], n)), nodes, edges)

    f = circular_2parallel_apply(g12, 5, 10)      # the throwing pair (e5, e10)
    @test f !== nothing
    R = f.newregs[2]                          # the gapless region carrying α₂/2
    @test isempty(regions(f.graph)[R].gaps)
    m = CircularMorphismGraph(f.graph, 0, 0)
    labels = fill(one(SoergelPoly), region_count(f.graph))
    labels[R] = (1//2) * alpha(2)
    fdm = CircularDecoratedMorphism(m, labels)

    # The dead end under V0: all edge neighbours lie at the same distance.
    distV = circular_region_distances(m)
    adj, _ = circular_region_adjacency(f.graph)
    @test all(distV[R2] == distV[R] for (R2, _) in adj[R])
    # Under the EDGE distance it is not: the region is one step further out than its
    # neighbours — exactly the `d − 1` neighbour the fusion needs.
    distE = circular_region_distances_edges_only(m)
    @test all(distE[R2] == distE[R] - 1 for (R2, _) in adj[R])
    @test circular_fusion_step(fdm) !== nothing

    # and the whole reduction arrives back at the NORMAL FORM of 1·g12. The round trip
    # compares against reduce_to_circular_leave(g12), not against g12 literally, because
    # g12 itself reduces further.
    res = reduce_to_circular_leave(fdm) +
          reduce_to_circular_leave(CircularDecoratedMorphism(m,
              let l = fill(one(SoergelPoly), region_count(f.graph))
                  l[f.newregs[1]] = (1//2) * alpha(2); l
              end))
    expected = reduce_to_circular_leave(with_marks(g12, CircularMorphismGraph(g12, 0, 0)))
    @test Dict(circular_canonical_key(d) => c for (d, c) in pairs_of(res)) ==
          Dict(circular_canonical_key(d) => c for (d, c) in pairs_of(expected))
end

@testset "phase A per round — no label on an inner region" begin
    # The rule: as soon as a diagram contains polynomials, the
    # polynomials should first be moved to the marker — and one level above
    # `reduce_circular`. The block
    # `circular_extract_scalars` → `circular_fusion_step` therefore runs at the head of EVERY
    # round of `reduce_to_circular_leave`, not only at the head of the call and after a D4
    # hit. That makes it gapless: no graph step sees a non-trivial label on an inner
    # region.
    for r in double_leaves([1, 2, 1, 3], [1, 2, 1, 3])
        fm = circular_morphism(r.morphism)
        fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly),
                                            region_count(fm.graph)))
        res = reduce_to_circular_leave(fdm)
        for (d, _) in pairs_of(res)
            m2 = CircularMorphismGraph(d.graph, fm.cut1, fm.cut2)
            dist = circular_region_distances(m2)
            for R in 1:region_count(d.graph)
                # neither a scalar (which `circular_extract_scalars` pulls out) nor a
                # genuine polynomial (which the fusion pushes outwards)
                @test isone(d.region_labels[R]) || dist[R] < 0
            end
        end
    end

    # The guard before the block (`any(!isone, region_labels)`) is exactly equivalent:
    # both steps require a non-trivial label. Without it, an extra
    # `circular_region_distances` BFS would run per term and round.
    let r = first(double_leaves([1, 2, 1, 3], [1, 2, 1, 3]))
        fm = circular_morphism(r.morphism)
        fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly),
                                            region_count(fm.graph)))
        factor, rest = circular_extract_scalars(fdm)
        @test isone(factor)                       # nothing to pull out
        @test rest.region_labels == fdm.region_labels
        @test circular_fusion_step(fdm) === nothing
    end
end

@testset "the growth guard in reduce_to_circular_leave" begin
    # If the node count grows monotonically in the reduce step, the driver should
    # THROW rather than carry on — a rip cord for the loop "connect -> reduce ->
    # search again", not a substitute for a proof.
    nd = circular_node([1, 1, 1])
    g = CircularGraph(CircularWord([1, 1, 1]), [nd],
        Edge[Edge(1, Leaf(1), NodePort(1, 1)), Edge(1, Leaf(2), NodePort(1, 2)),
             Edge(1, Leaf(3), NodePort(1, 3))])
    m = CircularMorphismGraph(g, 0, 0)
    fdm = CircularDecoratedMorphism(m, fill(one(SoergelPoly), region_count(g)))

    # The guard itself: `limit` counts INCREASES, so it needs limit+1 entries. Strictly
    # increasing throws, a plateau does not.
    @test_throws ErrorException DiagrammaticHecke._circular_leave_growth_guard([1, 2, 3, 4], 3, fdm)
    @test DiagrammaticHecke._circular_leave_growth_guard([1, 2, 2, 3], 3, fdm) === nothing
    @test DiagrammaticHecke._circular_leave_growth_guard([1, 2, 3, 4], 4, fdm) === nothing # too short
    @test DiagrammaticHecke._circular_leave_growth_guard([1, 2, 3, 4], 0, fdm) === nothing # switched off

    # On an ordinary diagram the guard changes nothing, not even at the tightest
    # limit: a single trivalent never grows.
    @test reduce_to_circular_leave(fdm; node_growth_limit = 1) == reduce_to_circular_leave(fdm)
end

@testset "fusion on the EDGE distance — the star node" begin
    # The fusion pushes polynomials step by step through edges, and for that it needs
    # the distance that simply decreases by 1. `circular_fusion_step` therefore
    # computes on `circular_region_distances_edges_only` rather than on V0.
    #
    # The fixture: a single `:gen13` star with boundary word [1,1,3,3,1,1,3,3] and
    # f = α₂² on the fifth sector. Under V0 the distances are squeezed together
    # ([0,1,1,2,2,2,1,1]) and there is no −1 neighbour, so the fusion cannot fire.
    # Under the edge distance ([0,1,2,3,4,3,2,1]) the neighbour is there: two terms,
    # `act(3,f)`
    # and `demazure(3,f)`, both of the same degree as the input.
    Leaf, NodePort, Edge = DiagrammaticHecke.Leaf, DiagrammaticHecke.NodePort, DiagrammaticHecke.Edge
    w = [1, 1, 3, 3, 1, 1, 3, 3]
    n = length(w); slot(k) = mod1(2 - k, n)
    a = Vector{Int}(undef, n)
    for k in 1:n; a[slot(k)] = w[k]; end
    g13 = CircularGraph(CircularWord(w), [circular_node(a)],
                   Edge[Edge(w[k], Leaf(k), NodePort(1, slot(k))) for k in 1:n])
    m13 = CircularMorphismGraph(g13, 0, 1)
    Rf = circular_node_sector_regions(g13)[1][5]
    l = fill(one(SoergelPoly), region_count(g13)); l[Rf] = alpha(2)^2

    @test circular_region_distances(m13)            == [0, 1, 1, 2, 2, 2, 1, 1]
    @test circular_region_distances_edges_only(m13) == [0, 1, 2, 3, 4, 3, 2, 1]

    fdm = CircularDecoratedMorphism(m13, l)
    r = circular_fusion_step(fdm)
    @test r !== nothing                                   # under V0 there is no hit
    ts = collect(pairs_of(r))
    @test length(ts) == 2
    # term 1: the polynomial moves through the 3-edge; term 2: the broken-open strand
    # with two fresh 3-dots carries the Demazure derivative.
    labs = sort([string(only(x for x in d.region_labels if !isone(x))) for (d, _) in ts])
    @test labs == sort([string(act(3, alpha(2)^2)), string(demazure(3, alpha(2)^2))])
    @test count(d -> length(d.graph.nodes) == 3, first.(ts)) == 1   # +2 Dots
    # degree-preserving, like every polynomial rule
    expected_deg = circular_degree(CircularDecorated(g13, l))
    @test all(circular_term_degree(c, d) == expected_deg for (d, c) in ts)
end

@testset "apply_circular_d4: borderless BC region, endA in the LAST slot" begin
    # `apply_circular_d4` has two ways of naming the region that inherits `α_i` in the
    # trivalent term: over the boundary GAP of the BC region, and — when that region
    # has NO gap — geometrically, as the sector of the new trivalent between the two
    # `ei` spokes. The geometric branch must not take the CYCLIC SUCCESSOR of `slotA`
    # as `mod1(slotA, 3) + 1`: that is correct for slots 1 and 2 but gives 4 for slot
    # 3, a `BoundsError` on the 3-element `ends_sorted`. This is the minimal figure
    # that reaches the branch.
    #
    # THE FIGURE: two [1,1,1] trivalents joined by TWO edges (a bigon of colour 1
    # — its interior is the region without a boundary gap), each with a third arm
    # to a leaf, plus a dot hanging on leaf 3 in the surrounding region. `ei` is
    # one of the bigon edges: it separates the dot's region from the borderless
    # bigon interior, and `_circular_d4_cyclic_ends` puts its end `e.a` in slot 3.
    NP = DiagrammaticHecke.NodePort
    g_bigon = CircularGraph(CircularWord([1, 1, 1]),
                       CircularNode[circular_node([1, 1, 1]), circular_node([1, 1, 1]), circular_node([1])],
                       Edge[Edge(1, NP(1, 1), NP(2, 2)),      # ei — the bigon edge
                            Edge(1, NP(1, 2), NP(2, 1)),
                            Edge(1, Leaf(1), NP(1, 3)),
                            Edge(1, Leaf(2), NP(2, 3)),
                            Edge(1, Leaf(3), NP(3, 1))])      # the dot's capping edge
    @test euler(g_bigon) == 2 && isempty(check_wiring(g_bigon))
    @test isempty(DiagrammaticHecke.regions(g_bigon)[2].gaps)   # the BC region is borderless

    res = apply_circular_d4(circular_decorated(g_bigon), 3, 1)   # the geometric branch
    @test res !== nothing
    ts = collect(pairs_of(res))
    @test length(ts) == 3                                   # 2 dot terms − 1 trivalent
    # degree balance per term, and every result really is a planar diagram
    @test all(circular_term_degree(c, d) == circular_degree(g_bigon) for (d, c) in ts)
    @test all(is_wired(d.graph) for (d, _) in ts)
    # In the trivalent term the `α₁` lands on the sector that HAS NO boundary gap —
    # the BC region of the input.
    triv = only(d for (d, c) in ts if c == -one(SoergelPoly))
    lab = findall(!isone, triv.region_labels)
    @test length(lab) == 1 && triv.region_labels[only(lab)] == alpha(1)
    @test isempty(DiagrammaticHecke.regions(triv.graph)[only(lab)].gaps)
end

@testset "apply_circular_d4: second bigon edge is an invalid D4 move" begin
    # The same figure as above, but `ei` is the OTHER bigon
    # edge. Geometrically this choice is invalid — the other edge would simplify
    # first to a needle, so the whole diagram is 0 by another rule. `apply_circular_d4`
    # must return `nothing` instead of throwing.
    NP = DiagrammaticHecke.NodePort
    g_bigon = CircularGraph(CircularWord([1, 1, 1]),
                       CircularNode[circular_node([1, 1, 1]), circular_node([1, 1, 1]), circular_node([1])],
                       Edge[Edge(1, NP(1, 1), NP(2, 2)),
                            Edge(1, NP(1, 2), NP(2, 1)),
                            Edge(1, Leaf(1), NP(1, 3)),
                            Edge(1, Leaf(2), NP(2, 3)),
                            Edge(1, Leaf(3), NP(3, 1))])
    @test apply_circular_d4(circular_decorated(g_bigon), 3, 2) === nothing
end

