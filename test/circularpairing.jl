# circularpairing.jl — the pairing ⟨DL_i, DL_j⟩_k as a circular-leaf
# linear combination (src/morphism/CircularPairing.jl).

@testset "CircularPairing / pairing(dl_i, dl_j)" begin
    x = [1, 2, 1]
    t = dl(x)

    # ---- Basic function: concatenate + decompose into circular leaves -----------------
    d0a, d0b = t[0][1], t[0][2]
    p = pairing(d0a, d0b)
    @test p isa CircularPairing
    @test p.x == x
    @test p.degree == d0a.degree + d0b.degree == 0
    @test all(!iszero(c) for c in values(p.coeffs))
    @test keys(p.coeffs) == keys(p.leaves)

    # self-composition of a degree-0 DL: a valid linear combination (possibly the empty
    # sum, e.g. when the idempotent part drops out — not an error, `CircularPairing` returns
    # it then too).
    ps = pairing(d0a, d0a)
    @test ps isa CircularPairing
    @test ps.degree == 0

    # not composable (different words) -> throws
    y = [1, 3, 2]
    ty = dl(y)
    @test_throws ErrorException pairing(d0a, ty[0][1])

    # ---- the "not yet perfect" notice ---------------------------------------
    # By test1_graded (test/pairing.jl), (121,121) is exact at EVERY degree (circular count ==
    # DL count), so NO notice may appear here.
    io = IOBuffer()
    show(io, MIME("text/plain"), p)
    out = String(take!(io))
    @test !occursin("⚠", out)
    @test occursin("CircularPairing", out)
end

