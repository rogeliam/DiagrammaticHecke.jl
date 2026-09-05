# The Soergel ring

`src/algebra/`. The coefficient ring of every diagram label is

```math
R = \mathbb{Q}[\alpha_1, \alpha_2, \alpha_3],
```

the polynomial ring on the simple roots of A₃, graded so that every `αᵢ` has
**Soergel degree 2**.

## Polynomials

[`alpha(i)`](@ref alpha) is the simple root `αᵢ`; [`SoergelPoly`](@ref) is the
element type and `R` the ring.

```julia
α1, α2, α3 = alpha(1), alpha(2), alpha(3)
p = 2 * α1^2 * α3 - α2
degree(p)          # 6
soergel_str(p)     # the compact notation used inside diagram pictures
poly_terms(p)      # sparse term dictionary, inverse of poly_from_terms
```

## The A₃ action

[`act(s, f)`](@ref act) applies the simple reflection `s` through the A₃ Cartan
matrix. It is an involution.

```julia
act(1, α1) == -α1
act(1, act(1, f)) == f
```

## Demazure operators

[`demazure(s, f)`](@ref demazure), with the alias `δ`, is

```math
\Delta_s(f) = \frac{f - s\cdot f}{\alpha_s}.
```

The standard identities hold: `δ_s` lowers the degree by 2, `δ_s ∘ δ_s = 0`,
and the image of `δ_s` is `s`-invariant.

```julia
demazure(1, α1)          # 2
δ(2, δ(2, f)) == 0
act(2, δ(2, f)) == δ(2, f)
```

## Linear algebra over R

The circular-leaf study asks one question about a coefficient matrix over `R`:
is it invertible? `R` is not a field, so the answer is computed division-free.

```julia
M = [α1 one(R) one(R); zero(R) α2 one(R); zero(R) zero(R) α3]
soergel_det(M)
is_permutation_triangular(M[[3,1,2], [2,3,1]])   # a cheap sufficient test
is_invertible_over_frac_ring(M)
soergel_det_peeled(M)                            # (det, dense-core size, computed)
```

[`is_permutation_triangular`](@ref) catches the common case without expanding a
determinant at all: after a permutation of rows and columns the matrix is
triangular, so the determinant is the product of the diagonal.
[`soergel_det_peeled`](@ref) peels those rows and columns off first and reports
how big the dense core it actually had to expand was.

## The Bott–Samelson pairing

[`bs_pairing(x, y)`](@ref bs_pairing) is the graded dimension of
`Hom(BS(x), BS(y))`, computed in the Hecke algebra rather than shipped as a
table.

```julia
bs_pairing([1,2,1], [1,2,1])                    # 0:2 2:5 4:4 6:1
expected_count([1,2,1], [1,2,1])                # the value at v = 1
expected_count([1,2,1], [1,2,1]; degree = 4)
pairing_table()                                 # all 24 x 24 short-lex pairs
```

This is the number the double-leaves basis has to match, and
[`ranks_match`](@ref) is the check that it does.

## Where the polynomials sit

In a diagram every **region** carries one element of `R`. That is the whole
role of this layer: [`CircularDecorated`](@ref) is a diagram plus one
`SoergelPoly` per region, and the rewrite rules move those labels around.

See [`examples/02-the-soergel-ring.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/02-the-soergel-ring.ipynb).
