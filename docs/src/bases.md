# The two bases

`Hom(BS(x), BS(y))` has two bases here, and moving between them is what the
package is for.

* the **double leaves** `DL`, glued from two light leaves — the classical basis;
* the **circular leaves** `CL`, the normal forms of the rewrite system.

## The double-leaves table

[`dl(x, y)`](@ref dl) builds every double leaf once and sorts them by degree.

```julia
x = [1, 2, 1]
t = dl(x, x)
dl_degrees(t), dl_ranks(t), expected_ranks(t)
ranks_match(t)          # the ranks agree with the Bott–Samelson pairing
t[2]                    # the basis elements of degree 2
dl_entries(t)           # all of them, flat
```

Each entry carries the two subexpressions `e`, `f` and the intermediate element
`z`, so it can be identified without looking at the picture.
[`circular_dl(x, y)`](@ref circular_dl) is the same basis read as circular
morphisms, index-parallel to `dl`.

## The circular-leaf basis

[`cl_reduce(m)`](@ref cl_reduce) reduces any morphism to a sum of circular
leaves; [`cl_basis(x, y)`](@ref cl_basis) collects the leaves that occur and
[`cl_coordinates`](@ref) reads a reduced sum off in those coordinates.

```julia
b = cl_basis(x, x)
b.leaves, b.degrees
coords = cl_coordinates(cl_reduce(h), b)
coords.coeffs
isempty(pairs_of(coords.rest))    # nothing left over
```

`rest` holds the terms that matched no basis element. An empty `rest` is the
statement that the answer is complete.

## The change of basis

[`dl_to_circular_matrix(x, y)`](@ref dl_to_circular_matrix) has one row per
double leaf and one column per circular leaf that occurs; row `k` is the
decomposition of the `k`-th double leaf.

```julia
M    = dl_to_circular_matrix(x, x)
Minv = circular_to_dl_matrix(x, x)      # `nothing` if it does not exist
circular_dl_triangular(circular_dl_matrix(x, x))
```

DL coordinates become CL coordinates by multiplying with `M`, and come back
through `Minv`. Two things go wrong in general, and both are the point rather
than a bug:

* **the matrix need not be square.** `M` has one row per double leaf and one
  column per circular leaf that occurs; for `21232` those numbers differ, and
  then `circular_to_dl_matrix` returns `nothing`.
* **the coefficients need not be scalars.** The first entries that are honest
  polynomials in `R` appear at `21232`; the small words give coefficient `1`
  throughout, which is an artefact of their size.

The smallest word where a double leaf is *not* a single circular leaf is
`1212`: one row of `M` has two nonzero entries, and the two leaves are
genuinely different diagrams, not two drawings of one.

## Coordinates of an arbitrary morphism

[`in_circular_dl_basis(h)`](@ref in_circular_dl_basis) writes any morphism in
the DL basis.

```julia
h = compose(f.morphism, g.morphism)     # `compose(f, g)` applies f FIRST
c = in_circular_dl_basis(h)
c.coeffs, c.exact, c.rest
```

`exact` and `rest` say how far to trust the answer, in the same sense as above.

## Composition and structure constants

[`compose_in_basis(f, g)`](@ref compose_in_basis) composes `f: x → y` with
`g: y → z` and returns the coordinates of the composite.
[`structure_constants(x, y, z)`](@ref structure_constants) runs that over every
pair.

```julia
sc = structure_constants([1, 2], [1, 2], [1, 2])
sc.exact, sc.table, sc.thrown
```

The table is keyed `(f, g)` with `f` applied first, so as a product the entry
is `g∘f`.

## The mirrors

[`hflip_dl_matrix(x, y)`](@ref hflip_dl_matrix) and
[`vflip_dl_matrix(x, y)`](@ref vflip_dl_matrix) decompose the two mirrors of
every basis element. `pure_swap` and `permutation` report the case where the
mirror only permutes the basis.

## The pairing

[`CircularPairing`](@ref) evaluates `⟨DLᵢ, DLⱼ⟩` on the circular side, which is
the same number [`bs_pairing`](@ref) computes in the Hecke algebra.

See [`examples/08-the-double-leaves-basis.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/08-the-double-leaves-basis.ipynb)
and [`examples/10-double-leaves-to-circular-leaves.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/10-double-leaves-to-circular-leaves.ipynb).
