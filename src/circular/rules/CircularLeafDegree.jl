# circular/rules/CircularLeafDegree.jl — the DEGREE invariant of the Circular-leaf reduction
# and Test 1 PER DEGREE.
#
# Comparing only the TOTAL counts #DLs vs. #distinct circular leaves is not
# enough. The Soergel-Hom formula
# (grdim Hom(BS(x), BS(y)), data/pairing/a3-bs-pairing.txt) gives the
# dimension PER DEGREE, though — here are the two tools for that:
#
#  1. `circular_degree_violations` — checks that the reduction PRESERVES degree:
#     for every term `p · FL` of a reduced DL,
#
#         degree(DL) == degree(p) + circular_degree(FL)
#
#     This is the formula "f·D-degree = 2·deg(f) + deg(D)": the `2·` is
#     already built into `degree(::SoergelPoly)` (Ring.jl — every `α_i` has
#     Soergel degree 2), and `circular_degree` (CircularGraph.jl) is `2·#colours −
#     #arms` per node (braid fixed at 0). Verified on (121,121)
#     (a variant with the coefficient degree doubled again cannot be told
#     apart as long as all coefficients are constant — the barbell probe
#     `α ↦ degree 2 = two consumed dots` settles on this form).
#
#  2. `test1_graded` — Test 1 per degree: per degree d, the number of DLs,
#     the number of distinct circular leaves, and (if the pairing table knows the
#     pair) the Soergel coefficient. This shows in WHICH degree an excess
#     sits.

"""
    circular_degree(d::CircularDecorated) -> Int

Degree of a decorated Circular graph: `circular_degree` of the graph PLUS the Soergel
degrees of all region labels (label `1` contributes 0). A circular leaf in DL form
normally has only label `1` — the labels still count here, so the degree
invariant also holds on intermediate states.
"""
circular_degree(d::CircularDecorated) =
    circular_degree(d.graph) + sum(degree(f) for f in d.region_labels if !isone(f); init = 0)

"""
    circular_term_degree(coeff::SoergelPoly, fd) -> Int

Degree of the term `coeff · fd` of a `CircularComboR`: `degree(coeff) +
circular_degree(fd)`. `fd` may be `CircularGraph` or `CircularDecorated`.
"""
circular_term_degree(coeff::SoergelPoly, fd) = degree(coeff) + circular_degree(fd)

"""
    circular_degree_violations(expected::Int, combo::CircularComboR)
        -> Vector{@NamedTuple{term::Any, coeff::SoergelPoly, got::Int}}

The DEGREE INVARIANT of the reduction, made checkable: for every term
`coeff · fd` of `combo` (typically: the result of
[`reduce_to_circular_leave`](@ref) applied to a DL), `circular_term_degree(coeff, fd)
== expected` must hold (`expected` = `degree` of the source DL). Returns the
VIOLATIONS — empty means: the invariant holds on every term.
"""
function circular_degree_violations(expected::Int, combo::CircularComboR)
    out = @NamedTuple{term::Any, coeff::SoergelPoly, got::Int}[]
    for (fd, c) in pairs_of(combo)
        got = circular_term_degree(c, fd)
        got == expected || push!(out, (term = fd, coeff = c, got = got))
    end
    return out
end

"""
    circular_leaf_degrees(combos; scalar_only = true)
        -> (degrees::Vector{Int}, index::Dict{Any,Int}, dropped::Int)

The DISTINCT circular leaves of a collection of reduced double leaves, per the
counting rule: only the number of DISTINCT circular leaves counts, and a leaf that occurs
several times but only with polynomial coefficients is not a circular leaf.

A circular leaf is a diagram WITHOUT a prefactor. A leaf that occurs in the whole
collection ONLY with a non-scalar coefficient isn't its own leaf there, just
a summand in the decomposition of a higher-degree double leaf — it doesn't
count. That's exactly what `scalar_only = true` (the default) does: a key
survives if it occurs AT LEAST ONCE with `degree(coeff) == 0`.

`combos` is iterable over `CircularComboR`s (typically one `reduce_to_circular_leave`
per double leaf). Returns

* `degrees` — for each surviving leaf, its `circular_degree` (representative at
  the FIRST occurrence), in order of first occurrence: exactly the
  `leaf_degrees` argument of [`test1_graded`](@ref);
* `index`   — key ↦ position in `degrees` (for building a matrix/Test 2);
* `dropped` — how many distinct keys the rule counted out.

`scalar_only = false` gives the counting scheme where every distinct key
counts (thus `dropped == 0`) — so both numbers can be compared side by
side in one run.
"""
function circular_leaf_degrees(combos; scalar_only::Bool = true)
    order  = Any[]                       # keys in order of first occurrence
    gdeg   = Dict{Any,Int}()             # key -> circular_degree
    scalar = Dict{Any,Bool}()            # key -> ever seen with a scalar coefficient?
    for combo in combos, (d, c) in pairs_of(combo)
        k = circular_canonical_key(d)
        if !haskey(gdeg, k)
            push!(order, k)
            gdeg[k]   = circular_degree(d)
            scalar[k] = false
        end
        scalar[k] |= iszero(degree(c))
    end
    keep    = scalar_only ? [k for k in order if scalar[k]] : order
    degrees = Int[gdeg[k] for k in keep]
    index   = Dict{Any,Int}(k => i for (i, k) in enumerate(keep))
    return (degrees = degrees, index = index, dropped = length(order) - length(keep))
end

"""
    test1_graded(x, y, dl_degrees, leaf_degrees)
        -> (rows, consistent)

Test 1 PER DEGREE for the word pair `(x, y)`:

* `dl_degrees`   — the degrees of ALL double leaves of the pair
                   (`[d.degree for d in double_leaves(x, y)]`),
* `leaf_degrees` — ONE degree per DISTINCT circular leaf of the reduction
                   (`circular_degree` of the representative at first occurrence).

`rows` is the table sorted by degree,
`@NamedTuple{degree, n_dl, n_distinct, pairing, excess}`:

* `pairing` — the coefficient of `grdim Hom(BS(x), BS(y))` at this degree
  from the Hecke pairing ([`bs_pairing`](@ref)). Typed `Union{Int,Missing}`,
  though always an `Int` in practice: the pairing is computed, not looked up in
  a table, so EVERY word pair is known.
* `excess = n_distinct − n_dl` — the excess in this degree (sign rule):
  **positive = more distinct circular leaves than DLs**, negative = leaves
  merge.

`consistent` is `true` when `n_dl == n_distinct` holds in EVERY degree and both
match the pairing coefficient there. If
the ungraded Test 1 passes but `consistent` doesn't, the excess sits in
degrees that cancel each other out.
"""
function test1_graded(x::Vector{Int}, y::Vector{Int},
                      dl_degrees::AbstractVector{Int},
                      leaf_degrees::AbstractVector{Int})
    pairing = bs_pairing(x, y)

    degs = sort(unique(vcat(collect(dl_degrees), collect(leaf_degrees),
                            collect(keys(pairing)))))
    T = @NamedTuple{degree::Int, n_dl::Int, n_distinct::Int,
                    pairing::Union{Int,Missing}, excess::Int}
    rows = T[]
    consistent = true
    for d in degs
        n_dl = count(==(d), dl_degrees)
        n_di = count(==(d), leaf_degrees)
        p = get(pairing, d, 0)
        ok = n_dl == n_di && n_dl == p
        consistent &= ok
        push!(rows, (degree = d, n_dl = n_dl, n_distinct = n_di,
                     pairing = p, excess = n_di - n_dl))
    end
    return (rows = rows, consistent = consistent)
end
