# test/circulargraphio.jl — circular/CircularGraphIO.jl. Round trip: save -> load ->
# hflip on the loaded object must agree with hflip on the original
# (no circular-composition analogue exists — see circular/CircularConvert.jl/Join.jl,
# `compose` exists only for `MorphismGraph`).

@testset "CircularGraphIO — save/load roundtrip (circular nodes)" begin
    # a merged CircularGraph (dots + a braid) built from a real conversion, so the
    # fixture is not hand-wired. e=[1,0,0,0,0,0] on both sides expresses [1] and
    # gives braid moves plus a genuine dot, so nodes/edges are non-trivial.
    w = double_leaf([1,2,1,3,2,1], [1,0,0,0,0,0], [1,2,1,3,2,1], [1,0,0,0,0,0])
    g0 = merge_all(circular(w.graph))
    @test !isempty(g0.nodes)

    # ---- CircularGraph roundtrip ----
    s = circulargraph_to_string(g0)
    g1 = circulargraph_from_string(s)
    @test g1 == g0                                   # circular_canonical_key equality
    @test g1.nodes == g0.nodes
    @test [nd.kind for nd in g1.nodes] == [nd.kind for nd in g0.nodes]
    @test is_wired(g1) == is_wired(g0)

    tmp = tempname() * ".fwg"
    save_circulargraph(tmp, g0)
    g2 = load_circulargraph(tmp)
    @test g2 == g0
    rm(tmp; force = true)

    # ---- CircularMorphismGraph roundtrip ----
    fm0 = circular_morphism(w)
    fm1 = circularmorphismgraph_from_string(circularmorphismgraph_to_string(fm0))
    @test fm1.graph == fm0.graph
    @test fm1.cut1 == fm0.cut1 && fm1.cut2 == fm0.cut2
    @test bottom(fm1) == bottom(fm0) && top(fm1) == top(fm0)

    tmpm = tempname() * ".fwg"
    save_circular_morphism(tmpm, fm0)
    fm2 = load_circular_morphism(tmpm)
    @test fm2.graph == fm0.graph && fm2.cut1 == fm0.cut1 && fm2.cut2 == fm0.cut2
    rm(tmpm; force = true)

    # a bare .fwg text (no `cuts` line) is an ERROR for the morphism reader
    @test_throws ErrorException circularmorphismgraph_from_string(circulargraph_to_string(g0))

    # hflip on the LOADED object must agree with hflip on the original — the
    # roundtrip test proper. No circular compose exists (checked: circular/CircularConvert.jl,
    # diagram/Join.jl only define hflip/flip/vflip for CircularMorphismGraph, no
    # compose analogue for CircularGraph/CircularMorphismGraph) — so the composition half
    # of the roundtrip is intentionally SKIPPED here.
    hf0 = hflip(fm0)
    hf1 = hflip(fm2)
    @test hf1.graph == hf0.graph
    @test hf1.cut1 == hf0.cut1 && hf1.cut2 == hf0.cut2
    @test bottom(hf1) == bottom(hf0) && top(hf1) == top(hf0)

    # ---- .fwgm batch container ----
    entries = [CircularMorphismEntry("case=1", circular_morphism(w)),
               CircularMorphismEntry("case=2", circular_morphism_graph(g0))]
    tmpb = tempname() * ".fwgm"
    save_circular_morphisms(tmpb, entries)
    loaded = load_circular_morphisms(tmpb)
    @test length(loaded) == 2
    @test [e.meta for e in loaded] == ["case=1", "case=2"]
    @test loaded[1].morphism.graph == entries[1].morphism.graph
    @test loaded[2].morphism.graph == entries[2].morphism.graph
    rm(tmpb; force = true)

    # ---- Zamo batch: one CircularMorphismGraph term per Zamo(i,j), i != j ----
    tmpz = tempname() * ".fwgm"
    written = save_zamo_circular_morphisms(tmpz; i_range = 1:2, j_range = 1:2)
    @test length(written) == 2                       # (1,2) and (2,1); (1,1)/(2,2) skipped
    @test Set(e.meta for e in written) == Set(["zamo=1,2", "zamo=2,1"])
    loadedz = load_circular_morphisms(tmpz)
    @test length(loadedz) == 2
    for (w1, w2) in zip(written, loadedz)
        @test w1.meta == w2.meta
        @test w1.morphism.graph == w2.morphism.graph
        @test w1.morphism.cut1 == w2.morphism.cut1 && w1.morphism.cut2 == w2.morphism.cut2
    end
    rm(tmpz; force = true)
end

