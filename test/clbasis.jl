# test/clbasis.jl — circular/CircularDecoratedIO.jl (`.cd`/`.cdb`) + morphism/CLBasis.jl.
# The word pair is SMALL (`121`/`121`) so the basis computation costs seconds.

@testset "SoergelPoly ↔ text (.cd building blocks)" begin
    @test soergelpoly_to_string(zero(SoergelPoly)) == "0"
    @test soergelpoly_to_string(one(SoergelPoly)) == "1"
    for p in (zero(SoergelPoly), one(SoergelPoly), alpha(1),
              alpha(1) * alpha(2) + soergel_monomial((0, 0, 3), 5 // 2),
              -alpha(3) - alpha(3))
        @test soergelpoly_from_string(soergelpoly_to_string(p)) == p
    end
end

@testset "CircularDecorated ↔ text (.cd roundtrip)" begin
    w = double_leaf([1, 2, 1], [1, 0, 0], [1, 2, 1], [1, 0, 0])
    g = merge_all(circular(w.graph))
    labels = fill(one(SoergelPoly), region_count(g))
    isempty(labels) || (labels[1] = alpha(2))
    d0 = CircularDecorated(g, labels, alpha(1) + alpha(3))

    d1 = circulardecorated_from_string(circulardecorated_to_string(d0))
    @test d1 == d0
    @test circular_canonical_key(d1) == circular_canonical_key(d0)
    @test d1.outer_label == d0.outer_label
    @test d1.region_labels == d0.region_labels

    tmp = tempname() * ".cd"
    save_circulardecorated(tmp, d0)
    @test load_circulardecorated(tmp) == d0
    rm(tmp; force = true)
end

@testset "CLBasis — computation, cache, coordinates" begin
    x = [1, 2, 1]; y = [1, 2, 1]
    dir = mktempdir()

    b = cl_basis(x, y; zamo = false, dir = dir, recompute = true)
    @test b.x == x && b.y == y
    @test b.zamo == false
    @test b.ndl == length(double_leaves(x, y))
    @test !isempty(b.leaves)
    @test length(b.leaves) == length(b.degrees)
    @test b.degrees == [circular_degree(d) for d in b.leaves]
    # keys are pairwise distinct — otherwise they would not be coordinates.
    @test length(cl_basis_index(b)) == length(b.leaves)

    # the file is there, and the loader returns the same
    path = cl_cache_path(x, y; zamo = false, dir = dir)
    @test isfile(path)
    b2 = load_cl_basis(path; zamo = false)
    @test b2 !== nothing
    @test [circular_canonical_key(d) for d in b2.leaves] == [circular_canonical_key(d) for d in b.leaves]
    @test b2.degrees == b.degrees && b2.ndl == b.ndl

    # SWITCH SETTING: a zamo-free file is no good as a zamo basis
    @test load_cl_basis(path; zamo = true) === nothing
    @test load_cl_basis(joinpath(dir, "does-not-exist.cdb")) === nothing
    # and the cache name distinguishes the two settings
    @test cl_cache_path(x, y; zamo = true, dir = dir) != path

    # degree slice
    bf = cl_basis_filter(b; degree_max = 0)
    @test all(<=(0), bf.degrees)
    @test length(bf.leaves) <= length(b.leaves)

    # coordinates: the identity lies in the basis, with a scalar coefficient
    c = cl_reduce(light_leaf_up(copy(x)); zamo = false)
    co = cl_coordinates(c, b)
    @test isempty(co.rest)
    @test count(!iszero, co.coeffs) >= 1
    @test any(p -> !iszero(p) && iszero(degree(p)), co.coeffs)

    rm(dir; force = true, recursive = true)
end

