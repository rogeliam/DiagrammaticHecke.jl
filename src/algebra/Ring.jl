# algebra/Ring.jl — the coefficient ring R = QQ[α₁,α₂,α₃] on top of AbstractAlgebra
#
# The polynomials are AbstractAlgebra multivariate polynomials over QQ; only the
# Soergel-specific layer lives here: the simple roots `alpha(i)`, the Soergel
# grading `degree` (every α_i has degree 2), the A₃ reflection `act` and the
# Demazure operator `demazure`/`δ`.
#
# AbstractAlgebra is pure Julia with no binary dependencies — arithmetic,
# normal form, `==`/`hash` and exact division come from there.  Only
# `AbstractAlgebra` is imported (qualified as `AA`); no name of it enters the
# `DiagrammaticHecke` namespace, so `degree`, `act` etc. stay ours.

import AbstractAlgebra as AA

"""
    R

The polynomial ring `QQ[α₁,α₂,α₃]`.  `zero(R)` / `one(R)` / `R(3//2)` produce
its elements.
"""
const R, (_ALPHA1, _ALPHA2, _ALPHA3) = AA.polynomial_ring(AA.QQ, [:α₁, :α₂, :α₃])

"""
    SoergelPoly

The element type of [`R`](@ref) — a polynomial in `α₁,α₂,α₃` with coefficients
in `QQ` (stored as `Rational{BigInt}`, so the coefficients cannot overflow).
"""
const SoergelPoly = AA.elem_type(R)

const _ALPHAS = (_ALPHA1, _ALPHA2, _ALPHA3)

# AbstractAlgebra deliberately offers no parent-free `zero`/`one`, because in
# general the parent cannot be recovered from the type.  Here there is exactly
# one ring, so `zero(SoergelPoly)` is unambiguous — and the ~180 call sites in
# this package keep working unchanged.
Base.zero(::Type{SoergelPoly}) = zero(R)
Base.one(::Type{SoergelPoly})  = one(R)

# ---------------------------------------------------------------------------
# Constructors and term access
# ---------------------------------------------------------------------------

"""
    soergel_monomial(e::NTuple{3,Int}, c) -> SoergelPoly

The single monomial `c · α₁^e₁ α₂^e₂ α₃^e₃`.
"""
soergel_monomial(e::NTuple{3,Int}, c) =
    R([Rational{BigInt}(c)], [collect(e)])

"""
    poly_terms(p::SoergelPoly) -> Dict{NTuple{3,Int}, Rational{BigInt}}

The sparse term dictionary of `p`: exponent vector ⇒ coefficient, zero
coefficients removed.  Inverse of [`poly_from_terms`](@ref).  Used by the
persistence format and by the pivot search in the DL writer.
"""
function poly_terms(p::SoergelPoly)
    d = Dict{NTuple{3,Int}, Rational{BigInt}}()
    for i in 1:length(p)
        e = AA.exponent_vector(p, i)
        d[(e[1], e[2], e[3])] = AA.coeff(p, i)
    end
    return d
end

"""
    poly_from_terms(d) -> SoergelPoly

Build an element of `R` from a term dictionary as produced by
[`poly_terms`](@ref).
"""
function poly_from_terms(d::AbstractDict)
    res = zero(R)
    for (e, c) in d
        iszero(c) && continue
        res += soergel_monomial((e[1], e[2], e[3]), c)
    end
    return res
end

"""
    constant_term(p::SoergelPoly) -> Union{Rational{BigInt}, Nothing}

The constant coefficient of `p` if `p` is constant, otherwise `nothing`.
Note `constant_term(zero(R)) == 0`, not `nothing`.
"""
constant_term(p::SoergelPoly) =
    AA.is_constant(p) ? AA.constant_coefficient(p) : nothing

# ---------------------------------------------------------------------------
# Pretty printing
# ---------------------------------------------------------------------------

"""
    soergel_str(p::SoergelPoly) -> String

`p` in the compact notation used in the diagram renderings: `2·α₁^2α₃` instead
of AbstractAlgebra's `2*α₁^2*α₃`, terms sorted by descending total degree.
"""
# Fallback for a DIFFERENT coefficient ring: `DecoratedDiagram{P}` is generic in `P`,
# and the pretty printer below knows only the A₃ exponent/coefficient shape. Any
# other ring prints through its own `show`.
soergel_str(p) = string(p)

function soergel_str(p::SoergelPoly)
    iszero(p) && return "0"
    sorted = sort!(collect(poly_terms(p));
                   by = t -> let (a, b, c) = t[1]; (-(a + b + c), -a, -b, -c) end)
    io = IOBuffer()
    first_term = true
    for (e, c) in sorted
        abs_c = abs(c)
        neg = c < 0
        first_term ? (neg && print(io, '-')) : print(io, ' ', neg ? '-' : '+', ' ')
        coeff_str = denominator(abs_c) == 1 ? string(numerator(abs_c)) : string(abs_c)
        if e == (0, 0, 0)
            print(io, coeff_str)
        else
            isone(abs_c) || print(io, coeff_str, '·')
            for (idx, exp) in enumerate(e)
                exp == 0 && continue
                print(io, 'α', ('₁', '₂', '₃')[idx])
                exp > 1 && print(io, '^', exp)
            end
        end
        first_term = false
    end
    return String(take!(io))
end

# ---------------------------------------------------------------------------
# Generators, grading, and the A₃ action
# ---------------------------------------------------------------------------

"""
    alpha(i::Int) -> SoergelPoly

The simple root `α_i`.  Only `i ∈ {1,2,3}` is allowed.
"""
function alpha(i::Int)
    i in 1:3 || throw(ArgumentError("alpha index must be 1, 2 or 3, got $i"))
    return _ALPHAS[i]
end

"""
    degree(p::SoergelPoly) -> Int

The Soergel grading: every `α_i` has degree `2`.  The zero polynomial has
degree `-1`.
"""
degree(p::SoergelPoly) = iszero(p) ? -1 : 2 * AA.total_degree(p)

# Image of the generator α_k under the simple reflection s_s, A₃ Cartan matrix:
#   s_i(α_i) = -α_i,   s_i(α_j) = α_j + α_i for |i-j| = 1,   s_i(α_j) = α_j else.
function _reflect_generator(s::Int, k::Int)
    s == k && return -_ALPHAS[k]
    abs(s - k) == 1 && return _ALPHAS[k] + _ALPHAS[s]
    return _ALPHAS[k]
end

"""
    act(s::Int, f::SoergelPoly) -> SoergelPoly

Apply the simple reflection `s_i` (`i ∈ {1,2,3}`) to `f`, using the A₃ Cartan
matrix.  The action is extended linearly and multiplicatively from the simple
roots — i.e. it is the substitution `α_k ↦ s_i(α_k)`.
"""
function act(s::Int, f::SoergelPoly)
    s in 1:3 || throw(ArgumentError("reflection index must be 1, 2 or 3, got $s"))
    return AA.evaluate(f, [_reflect_generator(s, k) for k in 1:3])
end

"""
    demazure(s::Int, f::SoergelPoly) -> SoergelPoly
    δ(s::Int, f::SoergelPoly) -> SoergelPoly

The Demazure operator `Δ_s(f) = (f - s·f) / α_s`.  The division is exact on the
image of `f ↦ f - s·f`.  `δ` is a short alias.
"""
function demazure(s::Int, f::SoergelPoly)
    s in 1:3 || throw(ArgumentError("Demazure index must be 1, 2 or 3, got $s"))
    return AA.divexact(f - act(s, f), _ALPHAS[s])
end

"""Short alias for [`demazure`](@ref)."""
const δ = demazure
