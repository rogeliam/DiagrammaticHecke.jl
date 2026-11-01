# test/circularmvalentwords.jl — the M-VALENT
# reading of the region words (`circular_region_mvalent_words`, switch
# `CIRCULAR_REGIONWORD_MVALENT`, default OFF). Node jumps at :gen2/:gen13 nodes
# plus the ε/fold collapse, the latter only at a TRIVALENT shared node.
#
# The fixtures are the target values: L4 `2,2,2,eps,2`, L5 keeps its transition,
# L7 loses it, `id_22` keeps it (a parallel bigon is no fold), and the −α₁ term of
# `id_1212` stays clean (the trivalent ε rule carries it).

@testset "circular_region_mvalent_words — the m-valent reading" begin
    @test CIRCULAR_REGIONWORD_MVALENT[] == false

    idfig(w) = let fm = circular_morphism(light_leaf_up(w))
        CircularMorphismGraph(fm.graph, fm.cut1, fm.cut2)
    end
    fixture(n) = let fp = joinpath(@__DIR__, "data", "plan36", "L$n.cd")
        isfile(fp) ? CircularMorphismGraph(load_circulardecorated(fp).graph, 0, 5) : nothing
    end

    # id_22: two PARALLEL 2-edges bound the middle region — they share no node,
    # so no fold, and there is no m-valent node either: both readings agree and
    # the trigger survives with the switch on.
    m22 = idfig([2, 2])
    wf, pf = circular_region_mvalent_words(m22)
    wl, pl = circular_region_words(m22)
    @test wf == wl && pf == pl

    old = CIRCULAR_REGIONWORD_MVALENT[]
    try
        CIRCULAR_REGIONWORD_MVALENT[] = true
        T = circular_unreduced_transition(m22)
        @test T !== nothing && T.case === :direct && T.colour == 2

        # L4 — the fixture for this reading. One 5-valent :gen2 node
        # whose sectors are all five regions; the target is `2,2,2,eps,2`
        # (the plain level-BFS reading gives `22,2,2,eps,22`).
        m4 = fixture(4)
        if m4 !== nothing
            w4, _ = circular_region_mvalent_words(m4)
            @test w4 == [[2], [2], [2], Int[], [2]]
            @test circular_unreduced_transition(m4) === nothing
        end

        # L5 KEEPS its transition: R1 hangs off the m-valent node by an own
        # 2-edge and reads `22`.
        m5 = fixture(5)
        if m5 !== nothing
            w5, _ = circular_region_mvalent_words(m5)
            @test w5 == [[2, 2], [2], [2], [2], Int[]]
            T5 = circular_unreduced_transition(m5)
            @test T5 !== nothing && T5.case === :direct && T5.colour == 2
        end

        # L6 unchanged (no m-valent node on the path), L7 loses its alarm.
        m6 = fixture(6)
        m6 === nothing || @test circular_unreduced_transition(m6) === nothing
        m7 = fixture(7)
        if m7 !== nothing
            w7, _ = circular_region_mvalent_words(m7)
            @test w7 == [[2], Int[], [2], [2], [2]]
            @test circular_unreduced_transition(m7) === nothing
        end

        # one rex-fusion step on id_1212, with the relation (121, 2). The terms
        # are picked by COEFFICIENT, not by position: the order of `pairs_of`
        # is not part of the relation.
        if isfile(circular_rex_cache_path([1, 2, 1], 2))
            CIRCULAR_REGIONWORD_MVALENT[] = false
            mid = idfig([1, 2, 1, 2])
            r = circular_rex_fusion_step(CircularDecoratedMorphism(
                    mid, fill(one(SoergelPoly), region_count(mid.graph))))
            CIRCULAR_REGIONWORD_MVALENT[] = true
            terms = collect(pairs_of(r.combo))
            @test length(terms) == 4
            term(c) = only(d for (d, k) in terms if k == c)

            # the −1 term, the false alarm of the level BFS, is silent
            m1 = CircularMorphismGraph(term(-one(SoergelPoly)).graph, mid.cut1, mid.cut2)
            @test circular_unreduced_transition(m1) === nothing
        end
    finally
        CIRCULAR_REGIONWORD_MVALENT[] = old
    end
end

