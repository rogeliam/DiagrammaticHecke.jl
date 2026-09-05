# test/circularexplosion.jl — the explosion recorder
# (src/circular/CircularExplosion.jl).
#
# What is checked is the MECHANISM, not an actual explosion: that an abort
# mid-reduction writes the diagram plus marker to disk, that the file is
# readable again, and that the error is rethrown unchanged.
# All within a temporary directory — `data/explosions/` stays untouched.

@testset "recording explosions (CircularExplosion.jl)" begin
    tmp = mktempdir()
    old = CIRCULAR_EXPLOSION_DIR[]
    CIRCULAR_EXPLOSION_DIR[] = tmp
    try
        fm  = circular_morphism(light_leaf_up([1, 2, 1, 2]))
        m   = CircularMorphismGraph(fm.graph, fm.cut1, fm.cut2)
        fdm = CircularDecoratedMorphism(m, fill(one(SoergelPoly), region_count(m.graph)))

        @testset "nothing is written without an error" begin
            circular_capture_explosion(name = "ok") do
                reduce_to_circular_leave(fdm)
            end
            @test isempty(circular_explosions())
        end

        @testset "an abort writes diagram + profile, and rethrows" begin
            @test_throws ErrorException circular_capture_explosion(name = "test",
                                                                   verbose = false) do
                reduce_to_circular_leave(fdm)
                error("kuenstlich")
            end
            fs = circular_explosions()
            # The facility writes the DEEPEST diagram, and a second `-last` file
            # when the last diagram is a different one — which it is as soon as the
            # reduction recurses (rex fusion starts a fresh top-level reduction).
            @test 1 <= length(fs) <= 2
            @test any(occursin("deepest", f) for f in fs)
            txt = read(replace(fs[1], ".cd" => ".txt"), String)
            @test occursin("kuenstlich", txt)
            # the marker belongs to it — it is NOT part of the .cd format
            @test occursin("cut1", txt) && occursin("cut2", txt)
            @test occursin("switches", txt)
            d = load_circulardecorated(fs[1])
            @test euler(d.graph) == 2
            @test isempty(check_wiring(d.graph))
        end

        @testset "the switch turns the recording off" begin
            for f in readdir(tmp; join = true); rm(f); end
            CIRCULAR_EXPLOSION_ENABLED[] = false
            try
                @test_throws ErrorException circular_capture_explosion(name = "off",
                                                                       verbose = false) do
                    reduce_to_circular_leave(fdm)
                    error("kuenstlich")
                end
            finally
                CIRCULAR_EXPLOSION_ENABLED[] = true
            end
            @test isempty(circular_explosions())
        end
    finally
        CIRCULAR_EXPLOSION_DIR[] = old
        rm(tmp; recursive = true, force = true)
    end
end

