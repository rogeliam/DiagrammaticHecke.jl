# test/circulardotmerge.jl — the three rules C19/C20/C23
# (src/circular/rules/CircularDotMerge.jl), the dot policy
# (`circular_dot_may_pass`) and the 2k-saturated region words
# (`circular_region_words_2k`).

# 2k-armed node (colours 1/2) with dots on the slots `dots`.
function _cdm_node(k::Int, dots::Vector{Int})
    n = 2k
    a = [isodd(s) ? 1 : 2 for s in 1:n]
    free = [s for s in [mod1(2 - j, n) for j in 1:n] if !(s in dots)]
    nodes = CircularNode[CircularNode(:braid, a)]
    edges = Edge[Edge(a[s], Leaf(j), NodePort(1, s)) for (j, s) in enumerate(free)]
    for s in dots
        push!(nodes, circular_node([a[s]]))
        push!(edges, Edge(a[s], NodePort(1, s), NodePort(length(nodes), 1)))
    end
    return CircularGraph(CircularWord([a[s] for s in free]), nodes, edges)
end
_cdm_g(x) = x isa CircularGraph ? x : x.graph

# The two diagrams on the boundary word 12122121.
function _cdm_d1()
    W = CircularWord([1,2,1,2,2,1,2,1])
    nodes = [CircularNode(:mixed, [1,1,1]), CircularNode(:mono, [2,2,2]),
             CircularNode(:braid, [1,2,1,2,1,2])]
    edges = [Edge(1, Leaf(1), NodePort(1,1)), Edge(1, Leaf(8), NodePort(1,2)),
             Edge(1, NodePort(1,3), NodePort(3,1)), Edge(2, Leaf(4), NodePort(2,1)),
             Edge(2, NodePort(2,2), NodePort(3,4)), Edge(2, Leaf(5), NodePort(2,3)),
             Edge(2, Leaf(7), NodePort(3,2)), Edge(1, Leaf(6), NodePort(3,3)),
             Edge(1, Leaf(3), NodePort(3,5)), Edge(2, Leaf(2), NodePort(3,6))]
    return CircularGraph(W, nodes, edges)
end
function _cdm_d2()
    W = CircularWord([1,2,1,2,2,1,2,1])
    nodes = [CircularNode(:braid, [1,2,1,2,1,2,1,2,1,2]),
             CircularNode(:mixed, [1]), CircularNode(:mono, [2])]
    edges = [Edge(1, Leaf(1), NodePort(1,5)), Edge(2, Leaf(2), NodePort(1,4)),
             Edge(1, Leaf(3), NodePort(1,3)), Edge(2, Leaf(4), NodePort(1,2)),
             Edge(2, Leaf(5), NodePort(1,10)), Edge(1, Leaf(6), NodePort(1,9)),
             Edge(2, Leaf(7), NodePort(1,8)), Edge(1, Leaf(8), NodePort(1,7)),
             Edge(1, NodePort(1,1), NodePort(2,1)), Edge(2, NodePort(1,6), NodePort(3,1))]
    return CircularGraph(W, nodes, edges)
end

@testset "dot policy: circular_dots_reducible" begin
    # (a) two adjacent dots, (b) one colour keeps fewer than 3 arms.
    # Everything else stays put — dots on 12-braids are allowed.
    @test CIRCULAR_DOT_POLICY[] == true
    @test CIRCULAR_REGIONWORD_DOT_POLICY === CIRCULAR_DOT_POLICY
    for (k, D, expected, desc) in ((4, [1],      false, "one dot"),
                              (4, [1,2],    true,  "(a) adjacent"),
                              (4, [8,1],    true,  "(a) across the seam"),
                              (4, [1,5],    true,  "(b) colour 1 down to 2"),
                              (4, [1,4],    false, "3+3, no neighbours"),
                              (5, [1,6],    false, "D2: 4+4"),
                              (5, [1,3],    false, "3+5"),
                              (5, [1,3,5],  true,  "(b) colour 1 down to 2"),
                              (5, [1,2],    true,  "(a) at 10 arms"))
        g = _cdm_node(k, D)
        @test circular_dots_reducible(g, 1) == expected
        # and the driver acts accordingly: it computes something exactly then
        _, hist = reduce_circular(g)
        @test isempty(hist) == !expected
    end
    # dotless and minimal nodes are never locked
    @test !circular_dots_reducible(_cdm_node(4, Int[]), 1)      # no dot: nothing to pull
    @test circular_dots_reducible(_cdm_node(3, [1]), 1)         # 6-armed = minimal
end

@testset "dot policy circular_dot_may_pass (the minimal clause)" begin
    @test circular_dot_may_pass(CircularNode(:braid, [1,2,1,2,1,2]))
    @test circular_dot_may_pass(CircularNode(:braid, [2,3,2,3,2,3]))
    @test circular_dot_may_pass(CircularNode(:mixed, [1,3,1,3]))
    @test !circular_dot_may_pass(CircularNode(:braid, [1,2,1,2,1,2,1,2]))
    @test !circular_dot_may_pass(CircularNode(:braid, [1,2,1,2,1,2,1,2,1,2]))
    @test !circular_dot_may_pass(CircularNode(:mono, [2,2,2]))
end

@testset "C20 two_adjacent_dots" begin
    # bit-identical to C14 + C5, at all four sites including the seam
    for (k, ds) in ((4,[1,2]), (5,[1,2]), (4,[8,1]), (4,[3,4]))
        g = _cdm_node(k, ds)
        h20 = first(pairs_of(_fr_two_adjacent_dots(g)))[1]
        c14 = first(pairs_of(DiagrammaticHecke._fr_dot_on_gen12(g; policy = false)))[1]
        r, hist = reduce_circular(c14)
        @test h20 == _cdm_g(first(pairs_of(r))[1])
        @test hist == [:merge]
        @test euler(h20) == 2 && isempty(check_wiring(h20))
        @test circular_degree(h20) == circular_degree(g)
        @test arm_count(h20.nodes[1]) == 2k - 2
    end
    # opposite dots are NOT a hit
    @test _fr_two_adjacent_dots(_cdm_node(4, [1,5])) === nothing
    # at the 6-armed node C20 does not fire (that's C6)
    @test _fr_two_adjacent_dots(_cdm_node(3, [1,2])) === nothing
end

@testset "C19 trivalent_into_gen12 — the inverse of C14" begin
    g1, g2 = _cdm_d1(), _cdm_d2()
    step(g) = (c = _fr_trivalent_into_gen12(g); c === nothing ? nothing : first(pairs_of(c))[1])
    k1 = step(g1); @test k1 !== nothing
    @test euler(k1) == 2 && isempty(check_wiring(k1))
    @test circular_degree(k1) == circular_degree(g1)
    k2 = step(k1); @test k2 !== nothing
    @test euler(k2) == 2 && isempty(check_wiring(k2))
    @test k2 == g2                                  # D1 --C19,C19--> D2
    @test step(k2) === nothing
    # and back: C14 twice — WITH `policy = false`, since under the policy D2 is exactly
    # a normal form: its two dots are not adjacent and leave both colours at 4 arms.
    @test DiagrammaticHecke._fr_dot_on_gen12(g2) === nothing
    h1 = first(pairs_of(DiagrammaticHecke._fr_dot_on_gen12(g2; policy = false)))[1]
    h2 = first(pairs_of(DiagrammaticHecke._fr_dot_on_gen12(h1; policy = false)))[1]
    @test h2 == g1
end

@testset "circular_region_words_2k — saturation at the 2k-node" begin
    @test CIRCULAR_REGIONWORD_2K[] == true
    @test CIRCULAR_REGIONWORD_DOT_POLICY[] == true
    @test DiagrammaticHecke.CIRCULAR_REGIONWORD_TRIVALENT_START[] == false

    gap(g) = (rs = regions(g); [findfirst(R -> k in R.gaps, rs) for k in 1:length(g.word)])
    m1 = CircularMorphismGraph(_cdm_d1(), 0, 4)
    m2 = CircularMorphismGraph(_cdm_d2(), 0, 4)
    w1 = circular_region_words_2k(m1); w2 = circular_region_words_2k(m2)
    @test [w1[R] for R in gap(m1.graph)] ==
          [[1], [1,2], [1,2,1], [1,2,1,2], [1,2,1], [1,2], [1], Int[]]
    # D2: instead of 1212, only 121 remains — all words reduced
    @test [w2[R] for R in gap(m2.graph)] ==
          [[1], [1,2], [1,2,1], [1,2,1], [1,2,1], [1,2], [1], Int[]]
    @test all(is_reduced(x) for x in w2 if x !== nothing)
    @test circular_unreduced_transition(m2) === nothing      # D2 counts as reduced
    @test circular_unreduced_transition(m1) !== nothing

    # the numerical check on the bare 20-armed node:
    # ε, 1, 2, 12, 21 and 121 fifteen times
    g20 = _cdm_node(10, Int[])
    w20 = circular_region_words_2k(CircularMorphismGraph(g20, 0, 10))
    @test count(x -> x == Int[], w20) == 1
    @test count(x -> length(x) == 1, w20) == 2
    @test count(x -> length(x) == 2, w20) == 2
    @test count(x -> x == [1,2,1], w20) == 15
end


# 2k-armed braid of colours `colours`, carrying a two-coloured {1,3} node with
# arm word `armword` on slot 1. `armword[1]` is the connector and must have the
# braid colour at slot 1; the remaining arms go to leaves. Traversal as in
# `_cdm_node`, so the boundary word stays clockwise-consistent.
function _cdm_mixed_at_braid(k::Int, colours::Tuple{Int,Int}, armword::Vector{Int})
    n = 2k
    a = [isodd(s) ? colours[1] : colours[2] for s in 1:n]
    m = length(armword)
    armword[1] == a[1] || error("the connector must carry the braid colour at slot 1")
    nodes = CircularNode[CircularNode(:braid, a)]; edges = Edge[]; word = Int[]
    for j in 1:n
        s = mod1(2 - j, n)
        if s == 1
            push!(nodes, circular_node(armword)); v = length(nodes)
            push!(edges, Edge(a[1], NodePort(1, 1), NodePort(v, 1)))
            for t in m:-1:2
                push!(word, armword[t])
                push!(edges, Edge(armword[t], Leaf(length(word)), NodePort(v, t)))
            end
        else
            push!(word, a[s])
            push!(edges, Edge(a[s], Leaf(length(word)), NodePort(1, s)))
        end
    end
    return CircularGraph(CircularWord(word), nodes, edges)
end



@testset "C23 mixed_into_gen12 — break the mixed {1,3} node open and merge" begin
    step(g) = (c = _fr_mixed_into_gen12(g); c === nothing ? nothing : collect(pairs_of(c)))

    @test CIRCULAR_MIXED_MERGE[] == true
    @test :mixed_into_gen12 in CIRCULAR_WEIGHT_ASSERT_EXEMPT

    # 121212-braid + 13113 over a 1-arm, and the same figure with the colours
    # swapped: 323232-braid + 31331 over a 3-arm.
    for (colours, armword) in (((1, 2), [1, 3, 1, 1, 3]), ((3, 2), [3, 1, 3, 3, 1]))
        g = _cdm_mixed_at_braid(3, colours, armword)
        @test euler(g) == 2 && isempty(check_wiring(g))

        ts = step(g)
        @test ts !== nothing && length(ts) == 1            # one term, not a sum
        h, co = ts[1]
        @test co == 1
        @test euler(h) == 2 && isempty(check_wiring(h))
        @test circular_degree(h) == circular_degree(g)     # the degree balance closes
        # 2(k+1) arms, one new dot, and the two pieces A and B
        @test arm_count(h.nodes[findfirst(nd -> nd.kind === :braid, h.nodes)]) == 8
        @test count(nd -> arm_count(nd) == 1, h.nodes) == 1
        @test count(nd -> nd.kind === :mixed, h.nodes) == 2

        # the one step IS "split, then C19" — compared by canonical key
        gs = _circular_split_mixed_at(g, 2, 1)
        @test gs !== nothing && isempty(check_wiring(gs)) && euler(gs) == 2
        @test circular_degree(gs) == circular_degree(g)    # the split alone is degree-neutral
        h2 = first(pairs_of(_fr_trivalent_into_gen12(gs)))[1]
        @test circular_canonical_key(h2) == circular_canonical_key(h)

        # at five arms each piece keeps only two arms of colour c, so the rule
        # does not fire again — that is the termination argument in miniature
        @test _fr_mixed_into_gen12(h) === nothing

        # in the driver: C23 runs (it is exempt from the weight assert) and shows
        # up in the history
        @test :mixed_into_gen12 in reduce_circular(circular_decorated(g))[2]
    end

    # 13113 splits as A = (3,1), B = (1,3): two 1313 nodes on the trivalent's legs
    h = step(_cdm_mixed_at_braid(3, (1, 2), [1, 3, 1, 1, 3]))[1][1]
    @test [sort(collect(arms(nd))) for nd in h.nodes if nd.kind === :mixed] ==
          [[1, 1, 3, 3], [1, 1, 3, 3]]

    # seven arms: piece B keeps three 1-arms, so the rule fires again — and
    # stops, because every piece is strictly smaller than the node it came from
    x = step(_cdm_mixed_at_braid(3, (1, 2), [1, 3, 1, 3, 1, 3, 1]))[1][1]
    n = 0
    while (r = _fr_mixed_into_gen12(x)) !== nothing && n < 10
        x = first(pairs_of(r))[1]; n += 1
    end
    @test 1 <= n < 10
    @test all(count(==(1), arms(nd)) <= 2
              for nd in x.nodes if nd.kind === :mixed && arm_count(nd) > 3)

    # no false alarm: fewer than three arms of the braid colour, and the switch off
    @test step(_cdm_mixed_at_braid(3, (1, 2), [1, 3, 1, 3])) === nothing
    @test step(_cdm_mixed_at_braid(3, (1, 2), [1, 1, 3, 3])) === nothing
    let saved = CIRCULAR_MIXED_MERGE[]
        try
            CIRCULAR_MIXED_MERGE[] = false
            @test step(_cdm_mixed_at_braid(3, (1, 2), [1, 3, 1, 1, 3])) === nothing
        finally
            CIRCULAR_MIXED_MERGE[] = saved
        end
    end
end
