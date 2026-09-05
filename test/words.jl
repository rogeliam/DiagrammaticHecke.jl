# words layer: CircularWord + Rules.

@testset "CircularWord (normal form, identity, order)" begin
    # rotation-equivalence: same element, same normal form
    @test CircularWord([1,3,1,2]) == CircularWord([3,1,2,1])
    @test compact(CircularWord([2,1])) == "12"          # canonical rotation
    @test compact(EMPTY) == "ε" && isempty(EMPTY)
    @test length(CircularWord([1,2,1])) == 3
    # stable total order: by length, then lexicographically
    @test EMPTY < CircularWord([1]) < CircularWord([2]) < CircularWord([1,1])
    # all rotations are listed (with repeats under cyclic symmetry)
    @test length(rotations(CircularWord([1,2,1,2]))) == 4
    @test_throws ArgumentError CircularWord([1,4])
end

@testset "rules (applied cyclically)" begin
    # worked example: 1312 --braid 121 on cyclic (3,4,1)--> 2321
    nbrs = Set(first(t) for t in neighbors(CircularWord([1,3,1,2])))
    @test CircularWord([2,3,2,1]) in nbrs

    # ii -> ε and ii -> i both exist
    n11 = Set(first(t) for t in neighbors(CircularWord([1,1])))
    @test EMPTY in n11
    @test CircularWord([1]) in n11

    # from ε, only increasing moves, and 11 is reachable
    nε = neighbors(EMPTY)
    @test all(sgn == +1 for (_, sgn) in nε)
    @test CircularWord([1,1]) in Set(first(t) for t in nε)

    # extended reducers (BASE_RULES do NOT contain these; EXTENDED_RULES does)
    ext = w -> moves(w, EXTENDED_RULES)
    @test any(m -> m.result == EMPTY && m.rule == :red_121212 && m.sign == -1,
              ext(CircularWord([1,2,1,2,1,2])))
    @test any(m -> m.result == CircularWord([2]) && m.rule == :red_12121,
              ext(CircularWord([1,2,1,2,1])))
    @test any(m -> m.result == CircularWord([2,1]) && m.rule == :red_1212,
              ext(CircularWord([1,2,1,2])))
    @test any(m -> m.result == EMPTY && m.rule == :red_232323,
              ext(CircularWord([2,3,2,3,2,3])))

    # down_neighbors uses the extended set but drops strictly-increasing moves
    @test isempty(down_neighbors(CircularWord([1])))      # 1 -> 11 is an increaser
    dns_121 = Set(first(t) for t in down_neighbors(CircularWord([1,2,1])))
    @test CircularWord([2,1,2]) in dns_121                # neutral braid is kept
    dns_121212 = Set(first(t) for t in down_neighbors(CircularWord([1,2,1,2,1,2])))
    @test EMPTY in dns_121212                             # reducing extra rule is kept
    @test all(sgn <= 0 for (_, sgn) in down_neighbors(CircularWord([1,2,1,2,1,2])))
end

@testset "braid-expander full expansions" begin
    # 1221 tiles fully but uses only {1,2}
    @test has_full_expansion(CircularWord([1,2,2,1]))
    @test !uses_all_letters(CircularWord([1,2,2,1]))
    @test uses_all_letters(CircularWord([1,2,3]))
    # a witness exists and its blocks are all braid-expanders
    ex = a_full_expansion(CircularWord([1,2,2,1]))
    @test ex !== nothing
    @test all(bt -> bt in BRAID_EXPANDERS, ex.blocks)

    # increasing_expansions returns every full increasing tiling (up to rotation)
    exps = increasing_expansions(CircularWord([1,2,2,1]))
    # 21211212 (12->2121, 21->1212), the docstring's witnessing expansion
    @test CircularWord([2,1,2,1,1,2,1,2]) in Set(exps)
    @test isempty(increasing_expansions(EMPTY))
end

@testset "up_neighbors (increases only, no ε-exit)" begin
    ups = Set(up_neighbors(CircularWord([1])))
    @test CircularWord([1,1]) in ups          # 1→11
    @test CircularWord([3,1,3]) in ups        # 1→313 (normal form 133)
    @test !(EMPTY in ups)                     # never decreases / hits ε
    @test all(length(x) > 1 for x in ups)     # every move increased length
    @test isempty(up_neighbors(EMPTY))        # ε-exits live in the seeds only
end
