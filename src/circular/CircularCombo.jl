# circular/CircularCombo.jl — formal R-linear combinations of CircularGraphs.
#
# `CircularCombo{T}` is keyed by `circular_canonical_key` and parametric in the
# coefficient ring `T`. `CircularComboR` fixes `T = SoergelPoly`: the circular
# rules produce α_i coefficients (barbell), which do not fit in `Int`.

"""
    CircularCombo{T}

A formal `T`-linear combination of circular diagrams, stored sparsely as
`circular_canonical_key => (representative, coefficient)` with zero coefficients
dropped. Terms with the same `circular_canonical_key` are merged — this is exactly
where isomorphic circular diagrams merge.

The representative is stored as `Any`, NOT as `CircularGraph`: once cell labels come
into play `CircularComboR` runs over `CircularDecorated`
(circular/CircularDecorated.jl), so the same combo type carries either
`CircularGraph` or `CircularDecorated` representatives, depending on how it was
constructed.
"""
struct CircularCombo{T}
    terms::Dict{Any,Tuple{Any,T}}
end

CircularCombo{T}() where {T} = CircularCombo{T}(Dict{Any,Tuple{Any,T}}())

"A single CircularGraph as a combination with coefficient `one(T)`."
function CircularCombo{T}(g::CircularGraph) where {T}
    c = CircularCombo{T}()
    c.terms[circular_canonical_key(g)] = (g, one(T))
    return c
end

Base.length(c::CircularCombo) = length(c.terms)
Base.isempty(c::CircularCombo) = isempty(c.terms)

"The coefficient of `g` in `c` (`zero(T)` if absent)."
coefficient(c::CircularCombo{T}, g::CircularGraph) where {T} =
    haskey(c.terms, circular_canonical_key(g)) ? c.terms[circular_canonical_key(g)][2] : zero(T)

"The (representative, coefficient) pairs, zero coefficients already excluded."
pairs_of(c::CircularCombo) = [(v[1], v[2]) for v in values(c.terms)]

"""
    _is_zero_diagram(g) -> Bool

If **any** region carries the label `0`, the whole diagram is `0` ("if a cell
has label 0, the whole diagram is 0").

Where this arises: `_circular_fuse_edge_regions` (circular/CircularDecorated.jl) sets
`glab * demazure(i, f)` in its second term. If `f` is constant, `demazure(i, f)
== 0` and the label becomes **exactly 0** — the term is the zero term. Without
this check, such a term keeps running and trips the guard in `apply_circular_d4`
("the dot's own region carries a non-trivial label"), even though the label
WAS zero.

A bare `CircularGraph` has no labels and is never zero for this reason.
"""
_is_zero_diagram(::CircularGraph) = false

# add `coeff * g` into `c` in place (isomorphic terms merge, 0 drops out).
function _add!(c::CircularCombo{T}, g, coeff::T) where {T}
    coeff == zero(T) && return c
    _is_zero_diagram(g) && return c
    k = circular_canonical_key(g)
    if haskey(c.terms, k)
        rep, old = c.terms[k]
        new = old + coeff
        new == zero(T) ? delete!(c.terms, k) : (c.terms[k] = (rep, new))
    else
        c.terms[k] = (g, coeff)
    end
    return c
end

function Base.:+(a::CircularCombo{T}, b::CircularCombo{T}) where {T}
    out = CircularCombo{T}(copy(a.terms))
    for (g, coeff) in pairs_of(b)
        _add!(out, g, coeff)
    end
    return out
end

Base.:-(a::CircularCombo{T}) where {T} = (zero(T) - one(T)) * a
Base.:-(a::CircularCombo{T}, b::CircularCombo{T}) where {T} = a + (-b)

function Base.:*(s::T, a::CircularCombo{T}) where {T}
    out = CircularCombo{T}()
    s == zero(T) && return out
    for (g, coeff) in pairs_of(a)
        _add!(out, g, s * coeff)
    end
    return out
end
Base.:*(a::CircularCombo{T}, s::T) where {T} = s * a

Base.:(==)(a::CircularCombo, b::CircularCombo) =
    Set(keys(a.terms)) == Set(keys(b.terms)) &&
    all(a.terms[k][2] == b.terms[k][2] for k in keys(a.terms))

_circular_rep_word(g::CircularGraph) = g.word
_circular_rep_word(d) = d.graph.word     # CircularDecorated (defined only in CircularDecorated.jl)

function Base.show(io::IO, c::CircularCombo{T}) where {T}
    if isempty(c)
        print(io, "0 (empty CircularCombo{$T})"); return
    end
    ps = sort(pairs_of(c); by = p -> letters(_circular_rep_word(p[1])))
    print(io, "CircularCombo{$T} with ", length(ps), " term(s):")
    for (g, coeff) in ps
        print(io, "\n  ", coeff, " · ", g)
    end
end

"""
    CircularComboR

Alias for `CircularCombo{SoergelPoly}` — the ring convention the circular rules (e.g.
barbell) compute in.
"""
const CircularComboR = CircularCombo{SoergelPoly}

# `SoergelPoly`-specific constructor/addition conveniences (an Int coefficient
# is auto-lifted to SoergelPoly), analogous to `DiagramComboR`.
CircularCombo{SoergelPoly}(g::CircularGraph) = CircularCombo{SoergelPoly}(Dict(circular_canonical_key(g) => (g, one(SoergelPoly))))
_add!(c::CircularComboR, g::CircularGraph, coeff::Integer) = _add!(c, g, coeff * one(SoergelPoly))
Base.:*(s::Integer, a::CircularComboR) = (s * one(SoergelPoly)) * a
Base.:*(a::CircularComboR, s::Integer) = s * a
