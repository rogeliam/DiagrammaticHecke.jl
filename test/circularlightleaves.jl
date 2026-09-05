# test/circularlightleaves.jl — circular LL/DL wrapper

@testset "Circular light/double leaves" begin
    x = [1, 2, 1]
    y = [1, 2, 1]
    z = [1]

    # circular_light_leaf gives a CircularMorphismGraph with the correct bottom/top
    e = [1, 1, 0]
    ll = circular_light_leaf(x, e)
    @test ll isa CircularMorphismGraph
    @test bottom(ll) == x
    @test top(ll) == expressed_word(x, e)

    # circular_double_leaf
    f = [1, 1, 0]
    dl = circular_double_leaf(x, e, y, f)
    @test dl isa CircularMorphismGraph
    @test bottom(dl) == x
    @test top(dl) == y

    # circular_light_leaves has the same count as light_leaves
    fll = circular_light_leaves(x, z)
    ll  = light_leaves(x, z)
    @test length(fll) == length(ll)
    @test all(fm isa CircularMorphismGraph for (_, fm) in fll)

    # circular_double_leaves has the same count as double_leaves
    fdl = circular_double_leaves(x, y)
    dl  = double_leaves(x, y)
    @test length(fdl) == length(dl)
    @test all(d.morphism isa CircularMorphismGraph for d in fdl)

    # every circular DL reduces without error (smoke test)
    for d in fdl
        combo = reduce_circular_full(d.morphism.graph)
        @test combo isa CircularComboR
    end

    # loading saved double leaves from a .wgm file. Written here and read back,
    # so the test carries its own fixture.
    path = joinpath(mktempdir(), "dl-1-1.wgm")
    dls  = double_leaves([1], [1])
    save_morphisms(path, [MorphismEntry("e=$(join(d.e)) f=$(join(d.f)) z=$(join(d.z)) " *
                                        "degree=$(d.degree)", d.morphism) for d in dls])

    loaded = circular_load_double_leaves(path)
    @test length(loaded) == length(dls)
    @test all(l.morphism isa CircularMorphismGraph for l in loaded)
    @test [l.degree for l in loaded] == [d.degree for d in dls]   # metadata survives
    for l in loaded                                               # every loaded DL reduces
        @test reduce_circular_full(l.morphism.graph) isa CircularComboR
    end
end

