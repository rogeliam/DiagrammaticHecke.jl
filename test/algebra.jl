# algebra layer: the Soergel ring R = QQ[α₁,α₂,α₃], the A₃ reflection /
# Demazure operators, division-free linear algebra over R, and the pairing
# table reader.

@testset "Soergel ring R = QQ[α₁,α₂,α₃]" begin
    α1 = alpha(1)
    α2 = alpha(2)
    α3 = alpha(3)

    # basic normal form and arithmetic
    @test α1 == α1
    @test α1 != α2
    @test α1 + α2 == α2 + α1
    @test α1 - α1 == zero(R)
    @test α1 * α2 == α2 * α1
    @test (α1 + α2)^2 == α1^2 + 2 * α1 * α2 + α2^2
    @test hash(α1 + α2) == hash(α2 + α1)

    # one/zero helpers
    @test iszero(zero(R))
    @test iszero(α1 - α1)
    @test one(R) * α1 == α1
    @test degree(α1 * α2^2) == 6
    @test degree(zero(R)) == -1

    # reflection is an involution
    f = α1^2 + α2 * α3 - 3 // 2 * α1
    @test act(1, act(1, f)) == f
    @test act(2, act(2, f)) == f
    @test act(3, act(3, f)) == f

    # explicit A₃ Cartan checks
    @test act(1, α1) == -α1
    @test act(1, α2) == α1 + α2
    @test act(1, α3) == α3
    @test act(2, α1) == α1 + α2
    @test act(2, α2) == -α2
    @test act(2, α3) == α2 + α3
    @test act(3, α2) == α2 + α3
    @test act(3, α3) == -α3

    # Demazure standard identities
    @test demazure(1, α1) == 2
    @test demazure(2, α2^2) == zero(R)
    @test demazure(3, α3^2) == zero(R)
    @test δ(1, α1) == 2           # δ is an alias for demazure

    # Δ_s ∘ Δ_s = 0
    for s in 1:3
        @test iszero(demazure(s, demazure(s, f)))
    end

    # image of Δ_s is s-invariant
    for s in 1:3
        @test act(s, demazure(s, f)) == demazure(s, f)
    end

    # low-degree Demazure values
    @test demazure(1, α2) == -1
    @test demazure(1, α3) == 0
    @test demazure(2, α1) == -1
    @test demazure(2, α2) == 2
    @test demazure(2, α3) == -1

    # a non-trivial Demazure computation
    g = α1 * α2 + α2 * α3
    @test demazure(2, g) == 2 * (α1 + α2 + α3)

    # term access round-trips; compact printing
    p = 2 * α1^2 * α3 - α2
    @test poly_from_terms(poly_terms(p)) == p
    @test constant_term(p) === nothing
    @test constant_term(R(3 // 2)) == 3 // 2
    @test soergel_str(p) == "2·α₁^2α₃ - α₂"
end

@testset "soergel_det / is_invertible_over_frac_ring" begin
    α1 = alpha(1); α2 = alpha(2)

    # 0×0 and 1×1
    @test soergel_det(Matrix{SoergelPoly}(undef, 0, 0)) == one(R)
    @test soergel_det(reshape([α1], 1, 1)) == α1
    @test is_invertible_over_frac_ring(reshape([α1], 1, 1))
    @test !is_invertible_over_frac_ring(reshape([zero(R)], 1, 1))

    # 2×2 invertible: [[α1, 1], [1, α2]], det = α1*α2 - 1 ≠ 0
    M2 = [α1 one(R); one(R) α2]
    @test soergel_det(M2) == α1 * α2 - one(R)
    @test is_invertible_over_frac_ring(M2)

    # 2×2 singular: second row is α1 * first row
    M2s = [α1 α2; α1^2 α1*α2]
    @test iszero(soergel_det(M2s))
    @test !is_invertible_over_frac_ring(M2s)

    # 3×3 invertible: upper triangular with nonzero diagonal, det = α1*α2*α3
    α3 = alpha(3)
    M3 = [α1 one(R) one(R); zero(R) α2 one(R); zero(R) zero(R) α3]
    @test soergel_det(M3) == α1 * α2 * α3
    @test is_invertible_over_frac_ring(M3)

    # non-square is never invertible
    @test !is_invertible_over_frac_ring([α1 α2])

    # is_permutation_triangular: detects triangular form UP TO row/column
    # permutation — the fast path for large, sparse, almost-permutation matrices
    # where Laplace expansion has no chance.
    @test is_permutation_triangular(M3)                        # already triangular
    @test is_permutation_triangular(M3[[3, 1, 2], [2, 3, 1]])  # permuted
    @test is_permutation_triangular([zero(R) one(R); one(R) zero(R)])  # permutation matrix
    @test !is_permutation_triangular([one(R) one(R); one(R) one(R)])   # not peelable
    @test !is_permutation_triangular([one(R) one(R); zero(R) zero(R)]) # zero row
    # The fast path must not change the ANSWER, only the route to it:
    @test is_invertible_over_frac_ring(M3[[3, 1, 2], [2, 3, 1]])
    @test !is_invertible_over_frac_ring([one(R) one(R); one(R) one(R)])

    # soergel_det_peeled agrees with soergel_det and reports the dense core
    r = soergel_det_peeled(M2)
    @test r.computed && r.det == soergel_det(M2)
end

@testset "pairing table reader" begin
    # The ground-truth table (computed independently with Oscar) loads
    # and the known (121,121) entry reads back: 0:2 2:5 4:4 6:1.
    table = pairing_table()
    @test haskey(table, ([1,2,1], [1,2,1]))
    @test table[([1,2,1], [1,2,1])] == Dict(0 => 2, 2 => 5, 4 => 4, 6 => 1)
    @test expected_count([1,2,1], [1,2,1]) == 12
    @test expected_count([1,2,1], [1,2,1]; degree = 4) == 4
    @test expected_count([1,2,1], [1,2,1]; degree = 8) == 0
    @test expected_count(Int[], Int[]) == 1        # ⟨BS(ε), BS(ε)⟩ = 1 in degree 0
end
