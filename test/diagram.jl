# diagram layer: graph / morphism text-format round-trips. (Faces, regions and
# the wiring check are covered by faces.jl.)

@testset "graph / morphism IO round-trips" begin
    g = r4_braid_dot()
    @test wordgraph_from_string(wordgraph_to_string(g)) isa WordGraph
    @test canonical_key(wordgraph_from_string(wordgraph_to_string(g))) == canonical_key(g)

    m = morphism_graph(braid(1, 3; m = 2), 0, 2)
    m2 = morphismgraph_from_string(morphismgraph_to_string(m))
    @test canonical_key(m2.graph) == canonical_key(m.graph)
    @test (m2.cut1, m2.cut2) == (m.cut1, m.cut2)
end
