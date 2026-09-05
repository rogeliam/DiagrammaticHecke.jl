# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/Decorated.jl — decorated diagrams
#
# A decorated diagram is a `WordGraph` whose inner cells each carry a polynomial
# `f ∈ R = QQ[α₁,α₂,α₃]`, plus the OUTER region (field `outer_label`).
#
# THE OUTER LABEL. The cell machinery (diagram/Faces.jl) gives the outer face NO
# index (`cell_of_face[outer] = 0`, `face_count` counts inner cells only). Treating
# a polynomial that lands there as the term's COEFFICIENT is only correct in the
# distance-0 case (then it really is a global scalar, see `extract_scalars`); in
# general it must stay outside, in place. Hence the OWN field — not an extra slot in
# `face_labels`, which would have shifted every cell index. `==`/`hash`/
# `canonical_key` take it into account: `key(α·D) != key(D)` for an α sitting
# outside. Every MARKED morphism is a distance-0 case (`outer_cell_distance` is
# constantly 0, the boundary circle is a wall), but an UNMARKED `DecoratedDiagram`
# has no distance, and there the label stays in place.
#
# INDEXING CONVENTION (canonical). `face_labels[c]` is the label of the cell with id
# `c`, i.e. exactly the order of `inner_faces(g)`: cell ids are handed out during
# face tracing by the smallest dart id of the face cycle (free circles last), a pure
# deterministic function of the graph structure (diagram/Faces.jl). Structurally
# equal graphs get the same cell numbering; `==`/`hash`/`canonical_key` compare the
# labels in that order (unlike the `cells` attribute, labels are NOT neutral).
#
# THE COEFFICIENT RING IS A TYPE PARAMETER. `DecoratedDiagram{P}` /
# `DiagramComboR{P}` carry the polynomial element type `P`; **`P = SoergelPoly` is the
# default everywhere**, so every existing A₃ call site keeps working unchanged and
# still gets `DecoratedDiagram{SoergelPoly}`. The parameter exists so that a second
# coefficient ring can reuse this layer.
# What `P` must provide: `zero(P)`, `one(P)` (AbstractAlgebra offers no parent-free
# `zero`/`one`, so the ring's own file defines them — `algebra/Ring.jl` does it for
# `SoergelPoly`), plus `+`, `*`, `-`, `iszero`, `isone`,
# `==` and `hash`. Nothing here touches the ring itself; a label is only ever added,
# multiplied and compared.
#
# `DiagramComboR` is the formal R-linear combination of such diagrams — the analogue
# of `DiagramCombo{Int}` over `DecoratedDiagram` instead of `WordGraph`, with
# coefficients in R (polynomials may be pulled out as scalar factors). Same storage:
# sparse canonical_key ⇒ (representative, coefficient), zero coefficients dropped,
# isomorphic terms merged. Int coefficients stay compatible (converted to SoergelPoly).

# ---- the decorated graph -----------------------------------------------------

"""
    DecoratedDiagram

A `WordGraph` together with a vector `face_labels` of polynomials (`SoergelPoly`),
one per inner cell (canonical order of `inner_faces`, see the file header), **plus
the label of the outer region** (`outer_label`, default `1`). The label of cell `c`
is `face_labels[c]`; a `1` (`one(R)`) means unlabelled.

The third argument is optional: `DecoratedDiagram(g, labels)` is
`DecoratedDiagram(g, labels, one(R))`.
"""
struct DecoratedDiagram{P}
    graph::WordGraph
    face_labels::Vector{P}
    outer_label::P

    function DecoratedDiagram(g::WordGraph, labels::Vector{P},
                              outer::P = one(P)) where {P}
        length(labels) == face_count(g) || throw(ArgumentError(
            "expected $(face_count(g)) labels (one per inner cell), got $(length(labels))"))
        return new{P}(g, labels, outer)
    end
end

"""
    poly_type(d::DecoratedDiagram) -> Type
    poly_type(c::DiagramComboR)    -> Type

The coefficient ring's element type this diagram (combination) is decorated with —
`SoergelPoly` by default, or whatever ring a second coefficient type provides.
"""
poly_type(::DecoratedDiagram{P}) where {P} = P

"""
    decorated(g::WordGraph) -> DecoratedDiagram
    decorated(g::WordGraph, labels::Vector{SoergelPoly}[, outer]) -> DecoratedDiagram
    decorated(g::WordGraph, f::SoergelPoly) -> DecoratedDiagram

The diagram `g` with cell labels: `1` everywhere if none are given; a single
polynomial `f` is shorthand when `g` has exactly one inner cell. `outer` is the
label of the outer region (default `1`).
"""
decorated(g::WordGraph) = decorated(g, SoergelPoly)
decorated(g::WordGraph, ::Type{P}) where {P} =
    DecoratedDiagram(g, fill(one(P), face_count(g)))
decorated(g::WordGraph, labels::Vector{P}, outer::P = one(P)) where {P} =
    DecoratedDiagram(g, labels, outer)
function decorated(g::WordGraph, f::P) where {P}
    face_count(g) == 1 || throw(ArgumentError(
        "decorated(g, f) shorthand needs exactly one inner cell; found $(face_count(g)) — use the vector form"))
    return DecoratedDiagram(g, [f])
end

boundary(d::DecoratedDiagram) = boundary(d.graph)
Base.length(d::DecoratedDiagram) = length(d.graph)
face_count(d::DecoratedDiagram) = face_count(d.graph)
inner_faces(d::DecoratedDiagram) = inner_faces(d.graph)

"Label of the inner cell `c` (canonical id from `inner_faces`)."
face_label(d::DecoratedDiagram, c::Int) = d.face_labels[c]

"""
    outer_label(d::DecoratedDiagram) -> SoergelPoly

The label of the OUTER region (see the file header). `1` means unlabelled. The
outer region has no cell id, so it is NOT reachable via `face_label(d, c)`.
"""
outer_label(d::DecoratedDiagram) = d.outer_label

# Structural comparison without the `cells` attribute (a function of the structure
# anyway — excluded explicitly so it is clear what gets compared).
function Base.:(==)(a::DecoratedDiagram, b::DecoratedDiagram)
    a.face_labels == b.face_labels || return false
    a.outer_label == b.outer_label || return false
    return (a.graph.word, a.graph.nodes, a.graph.edges) ==
           (b.graph.word, b.graph.nodes, b.graph.edges)
end

Base.hash(d::DecoratedDiagram, h::UInt) =
    hash((d.graph.word, d.graph.nodes, d.graph.edges, d.face_labels, d.outer_label),
         hash(:DecoratedDiagram, h))

"""
    canonical_key(d::DecoratedDiagram)

The key under which decorated diagrams merge inside a `DiagramComboR`: the
structural `canonical_key` of the graph, the label sequence in canonical cell order
**and the outer label**. Same structure (same cell numbering) + same labels ⇒ same
term.

The outer label is part of it: if an `α` sits outside in place, `α·D` is a DIFFERENT
term from `D`.
"""
canonical_key(d::DecoratedDiagram) =
    (canonical_key(d.graph), Tuple(d.face_labels), d.outer_label)

function Base.show(io::IO, d::DecoratedDiagram)
    print(io, "DecoratedDiagram(", d.graph)
    for c in 1:face_count(d)
        f = d.face_labels[c]
        isone(f) && continue
        print(io, ", [", c, "] ↦ ", f)
    end
    isone(d.outer_label) || print(io, ", [outer] ↦ ", d.outer_label)
    print(io, ")")
end

# ---- formal R-linear combinations --------------------------------------------

"""
    DiagramComboR

A formal R-linear combination of decorated diagrams (`DecoratedDiagram`), stored
sparsely as `canonical_key ⇒ (representative, SoergelPoly coefficient)` with zero
coefficients dropped. Isomorphic terms (same structure + same cell labels in
canonical order) merge on addition. Int coefficients are converted to `SoergelPoly`
internally.
"""
struct DiagramComboR{P}
    terms::Dict{Any, Tuple{DecoratedDiagram{P}, P}}
end

DiagramComboR{P}() where {P} = DiagramComboR{P}(Dict{Any, Tuple{DecoratedDiagram{P}, P}}())
DiagramComboR() = DiagramComboR{SoergelPoly}()

"A single decorated diagram as a combination with coefficient 1."
DiagramComboR(d::DecoratedDiagram{P}) where {P} =
    DiagramComboR{P}(Dict{Any, Tuple{DecoratedDiagram{P}, P}}(
        canonical_key(d) => (d, one(P))))

poly_type(::DiagramComboR{P}) where {P} = P

"The neutral (empty) combination over the same ring as `c`."
Base.zero(c::DiagramComboR{P}) where {P} = DiagramComboR{P}()

Base.length(c::DiagramComboR) = length(c.terms)
Base.isempty(c::DiagramComboR) = isempty(c.terms)

"The coefficient of `d` in `c` (0 ∈ R if absent)."
coefficient(c::DiagramComboR{P}, d::DecoratedDiagram) where {P} =
    haskey(c.terms, canonical_key(d)) ? c.terms[canonical_key(d)][2] : zero(P)

"The (diagram, coefficient) pairs, zero coefficients already excluded."
pairs_of(c::DiagramComboR) = [(v[1], v[2]) for v in values(c.terms)]

# add `coeff * d` into `c` in place (isomorphic terms merge, zeros are dropped)
function _add!(c::DiagramComboR{P}, d::DecoratedDiagram{P}, coeff::P) where {P}
    iszero(coeff) && return c
    k = canonical_key(d)
    if haskey(c.terms, k)
        rep, old = c.terms[k]
        new = old + coeff
        iszero(new) ? delete!(c.terms, k) : (c.terms[k] = (rep, new))
    else
        c.terms[k] = (d, coeff)
    end
    return c
end
_add!(c::DiagramComboR{P}, d::DecoratedDiagram{P}, coeff::Integer) where {P} =
    _add!(c, d, coeff * one(P))

function Base.:+(a::DiagramComboR{P}, b::DiagramComboR{P}) where {P}
    out = DiagramComboR{P}(copy(a.terms))
    for (d, coeff) in pairs_of(b)
        _add!(out, d, coeff)
    end
    return out
end

Base.:-(a::DiagramComboR{P}) where {P} = (-one(P)) * a
Base.:-(a::DiagramComboR{P}, b::DiagramComboR{P}) where {P} = a + (-b)

function Base.:*(s::P, a::DiagramComboR{P}) where {P}
    out = DiagramComboR{P}()
    iszero(s) && return out
    for (d, coeff) in pairs_of(a)
        _add!(out, d, s * coeff)
    end
    return out
end
Base.:*(a::DiagramComboR{P}, s::P) where {P} = s * a
Base.:*(s::Integer, a::DiagramComboR{P}) where {P} = (s * one(P)) * a
Base.:*(a::DiagramComboR, s::Integer) = s * a

Base.:(==)(a::DiagramComboR, b::DiagramComboR) =
    Set(keys(a.terms)) == Set(keys(b.terms)) &&
    all(a.terms[k][2] == b.terms[k][2] for k in keys(a.terms))

function Base.show(io::IO, c::DiagramComboR)
    if isempty(c)
        print(io, "0 (empty DiagramComboR)"); return
    end
    ps = sort(pairs_of(c); by = p -> letters(p[1].graph.word))
    print(io, "DiagramComboR with ", length(ps), " term(s):")
    for (d, coeff) in ps
        print(io, "\n  ", coeff, " · ", d)
    end
end
