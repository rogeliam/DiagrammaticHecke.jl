# Focused checks for decorated diagrams: a WordGraph with a polynomial per inner cell,
# plus formal ℤ-combinations.

@testset "decorated diagrams" begin
    α1 = alpha(1)
    α2 = alpha(2)

    # a single strand (chord) has exactly two inner cells
    chord = WordGraph(CircularWord([1, 1]), Node[], Edge[Edge(1, Leaf(1), Leaf(2))])
    @test face_count(chord) == 2

    # default decoration: every cell labelled 1
    d0 = decorated(chord)
    @test d0.face_labels == [one(R), one(R)]
    @test face_label(d0, 1) == one(R)

    # explicit labels; wrong length is rejected
    d = decorated(chord, [α1, α2])
    @test face_label(d, 1) == α1
    @test face_label(d, 2) == α2
    @test_throws ArgumentError decorated(chord, [α1])

    # the single-polynomial shortcut needs exactly one inner cell
    onecell = dot(1)                       # dot has one inner cell
    @test face_count(onecell) == 1
    @test face_label(decorated(onecell, α1), 1) == α1
    @test_throws ArgumentError decorated(chord, α1)   # two cells → ambiguous

    # equality / hash key off BOTH structure and the label sequence
    @test decorated(chord, [α1, α2]) == decorated(chord, [α1, α2])
    @test decorated(chord, [α1, α2]) != decorated(chord, [α2, α1])
    @test hash(decorated(chord, [α1, α2])) == hash(decorated(chord, [α1, α2]))

    # ---- DiagramComboR arithmetic ----
    ca = DiagramComboR(d)
    @test length(ca) == 1
    @test coefficient(ca, d) == 1

    # same decorated term merges (coefficients add)
    @test coefficient(ca + ca, d) == 2

    # differing labels stay separate terms
    d2 = decorated(chord, [α2, α1])
    @test length(DiagramComboR(d) + DiagramComboR(d2)) == 2

    # zero coefficients are dropped: d - d = empty sum
    @test isempty(ca - ca)

    # scalar multiplication
    @test coefficient(3 * ca, d) == 3
    @test isempty(0 * ca)
    @test coefficient(-ca, d) == -1

    # decorated_svg writes each cell's label into the SVG
    svg = decorated_svg(decorated(onecell, α1))
    @test occursin("α₁", svg)
end

@testset "DecoratedDiagram image display" begin
    # DecoratedDiagram implements Base.show/showable for both image MIME types
    # (svg+xml and png, mirroring render/Tutte.jl), so cell labels (α_i) are visible
    # even where only image/png is read, unlike WordGraph/MorphismGraph.
    α1 = alpha(1)
    chord = WordGraph(CircularWord([1, 1]), Node[], Edge[Edge(1, Leaf(1), Leaf(2))])
    d = decorated(chord, [α1, one(R)])

    @test showable(MIME"image/svg+xml"(), d)
    @test showable(MIME"image/png"(), d)

    svg = sprint(show, MIME"image/svg+xml"(), d)
    @test svg == decorated_svg(d)
    @test occursin("α₁", svg)

    png = codeunits(sprint(show, MIME"image/png"(), d))
    @test !isempty(png)
    @test png[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]  # PNG magic bytes
end

# ---------------------------------------------------------------------------
# The coefficient ring as a TYPE PARAMETER.
#
# `DecoratedDiagram{P}` / `DiagramComboR{P}` must work over a ring OTHER than
# ℚ[α₁,α₂,α₃] without A₃ noticing anything. The touchstone is a second coefficient
# ring, here a quadratic extension of ℚ, K = ℚ(φ) with φ² = φ+1, to check that
# nothing in this layer assumes `SoergelPoly`. It is built INLINE here so the test
# does not depend on an external ring implementation.
#
# The two `Base.zero`/`Base.one` methods sit at TOP LEVEL, not inside the
# `@testset`: method definitions in a local scope only become visible in a newer
# world (world age) and could not be called from the same block. These two
# methods are exactly the contract a ring has to fulfil — `algebra/Ring.jl` does
# it for `SoergelPoly`, and this file does it inline for K.

const _AA5 = DiagrammaticHecke.AA
const _QQx5, _x5 = _AA5.polynomial_ring(_AA5.QQ, :φ)
const _K5, _kmap5 = _AA5.residue_field(_QQx5, _x5^2 - _x5 - 1)
const _PHI5 = _kmap5(_x5)
const _R5, (_B1, _B2, _B3) = _AA5.polynomial_ring(_K5, [:α₁, :α₂, :α₃])
const _P5 = _AA5.elem_type(_R5)

Base.zero(::Type{_P5}) = zero(_R5)
Base.one(::Type{_P5})  = one(_R5)

@testset "second coefficient ring (type parameter P)" begin
    chord = WordGraph(CircularWord([1, 1]), Node[], Edge[Edge(1, Leaf(1), Leaf(2))])
    @test face_count(chord) == 2

    # A₃ stays the default: with no ring given you get SoergelPoly
    dA = decorated(chord)
    @test dA isa DecoratedDiagram{SoergelPoly}
    @test poly_type(dA) === SoergelPoly

    # the same graph over K
    dK = DecoratedDiagram(chord, [_B2, one(_P5)])
    @test dK isa DecoratedDiagram{_P5}
    @test poly_type(dK) === _P5
    @test face_label(dK, 1) == _B2
    @test isone(outer_label(dK))
    @test decorated(chord, _P5) == DecoratedDiagram(chord, [one(_P5), one(_P5)])

    # labels over K separate terms just as they do over ℚ
    @test dK != DecoratedDiagram(chord, [_B1, one(_P5)])
    @test canonical_key(dK) != canonical_key(decorated(chord, _P5))

    # linear combination over K: φ·D + D = φ²·D — φ² = φ+1 is the touchstone
    cK = DiagramComboR(dK)
    @test cK isa DiagramComboR{_P5}
    @test poly_type(cK) === _P5
    cK2 = (_PHI5 * one(_P5)) * cK + cK
    @test length(cK2) == 1
    @test coefficient(cK2, dK) == (_PHI5^2) * one(_P5)
    @test coefficient(cK2, dK) == (_PHI5 + 1) * one(_P5)
    @test isempty(cK2 - cK2)
    @test iszero(coefficient(cK2 - cK2, dK))
    @test coefficient(zero(cK), dK) == zero(_P5)
    @test length(2 * cK) == 1 && coefficient(2 * cK, dK) == 2 * one(_P5)

    # and A₃ keeps computing unchanged
    cA = DiagramComboR(dA)
    @test cA isa DiagramComboR{SoergelPoly}
    @test coefficient(alpha(1) * cA, dA) == alpha(1)
end

