# algebra/SoergelMatrix.jl — determinants and invertibility over R = QQ[α₁,α₂,α₃]
#
# The one question the circular-leaf study asks about a coefficient matrix `M`
# with entries in `R`: is `M` invertible over the fraction field `Frac(R)`?
# `R` is a polynomial ring over a field, hence an integral domain, so this
# reduces to a question that stays inside `R`: `det(M) != 0`.
#
# The determinant itself comes from AbstractAlgebra (fraction-free elimination
# over the integral domain `R`).  What this file adds is the *sparsity* layer
# the AbstractAlgebra routines know nothing about: the matrices here are
# extremely thin — the 244×244 matrix of the `(123212, 123212)` sweep is 0.42 %
# dense, 237 of its 244 rows carry a single entry — and for such a matrix the
# determinant can be read off by peeling, without any elimination at all.

"""
    soergel_det(M::AbstractMatrix{SoergelPoly}) -> SoergelPoly

Determinant of a square matrix over `R = QQ[α₁,α₂,α₃]`, computed by
AbstractAlgebra over the integral domain `R`.  `1×1` returns the entry, `0×0`
returns `one(R)` (the empty-product convention).
"""
function soergel_det(M::AbstractMatrix{SoergelPoly})
    n, m = size(M)
    n == m || throw(ArgumentError("soergel_det: matrix must be square, got $(n)×$(m)"))
    n == 0 && return one(R)
    n == 1 && return M[1, 1]
    return AA.det(AA.matrix(R, Matrix{SoergelPoly}(M)))
end

"""
    is_permutation_triangular(M::AbstractMatrix{SoergelPoly}) -> Bool

`true` iff the rows and columns of the square matrix `M` can be permuted into
TRIANGULAR form with a nonzero diagonal.

**Method (peeling).** Repeatedly find a still-live row with EXACTLY ONE live
nonzero entry, pair it with that column, and strike both; if none exists, do
the same column-wise.  If this succeeds `n` times, `M` is triangular up to
permutation.

**Exact, not a heuristic.** Each peeling step is a Laplace expansion along a
row (or column) with a single entry: `det(M) = ±M[i,j]·det(Minor)`.  After `n`
steps, `det(M) = ±∏` of the paired entries, all of which are `≠ 0`.
"""
function is_permutation_triangular(M::AbstractMatrix{SoergelPoly})
    n, m = size(M)
    (n == m && n > 0) || return n == m
    row_alive = trues(n)
    col_alive = trues(n)
    for _ in 1:n
        found = false
        for i in 1:n                          # row with exactly one entry
            row_alive[i] || continue
            j0 = 0; cnt = 0
            for j in 1:n
                (col_alive[j] && !iszero(M[i, j])) || continue
                cnt += 1
                cnt > 1 && break
                j0 = j
            end
            if cnt == 1
                row_alive[i] = false; col_alive[j0] = false; found = true; break
            end
        end
        if !found
            for j in 1:n                      # else: column with exactly one entry
                col_alive[j] || continue
                i0 = 0; cnt = 0
                for i in 1:n
                    (row_alive[i] && !iszero(M[i, j])) || continue
                    cnt += 1
                    cnt > 1 && break
                    i0 = i
                end
                if cnt == 1
                    row_alive[i0] = false; col_alive[j] = false; found = true; break
                end
            end
        end
        found || return false
    end
    return true
end

"""
    soergel_det_peeled(M::AbstractMatrix{SoergelPoly}; maxcore = typemax(Int))
        -> (det, core, computed)

Determinant over `R = QQ[α₁,α₂,α₃]` that **exploits sparsity** before handing
the rest to [`soergel_det`](@ref).

**Method.** As long as there is a live row (or column) with EXACTLY ONE live
nonzero entry `M[i,j]`, that is a Laplace expansion along this row with a
single summand:

    det(M) = (−1)^(r+c) · M[i,j] · det(Minor)

`r`/`c` are the positions of `i`/`j` WITHIN the still-live rows/columns
respectively — the detail that is easy to get wrong.  The factor is
accumulated, the row and column struck, and the process repeats.  This is
**exact**, not an approximation.

**Termination cases.**
* A live row or column that is all zeros ⇒ `det = 0`, done immediately.
* No single-entry row/column left ⇒ a dense **core** remains, whose
  determinant is computed with [`soergel_det`](@ref).  `core` reports its
  size.  `maxcore` bounds that step: a core larger than the bound makes the
  function bail out with `computed = false` (`det` is then meaningless).  The
  default is no bound — the elimination copes with the cores that occur here;
  the parameter is kept for callers that would rather skip than wait.
"""
function soergel_det_peeled(M::AbstractMatrix{SoergelPoly};
                            maxcore::Int = typemax(Int))
    n, m = size(M)
    n == m || throw(ArgumentError(
        "soergel_det_peeled: matrix must be square, got $(n)×$(m)"))
    n == 0 && return (det = one(R), core = 0, computed = true)

    row_alive = trues(n)
    col_alive = trues(n)
    acc  = one(R)
    sgn  = 1
    left = n

    # position of `k` among the live indices (1-based)
    rank_of(k, alive) = count(x -> alive[x], 1:k)

    while left > 0
        hit = nothing                          # (i, j)
        # 1. row with exactly one live entry — also detects an all-zero row
        for i in 1:n
            row_alive[i] || continue
            j0 = 0; cnt = 0
            for j in 1:n
                (col_alive[j] && !iszero(M[i, j])) || continue
                cnt += 1
                cnt > 1 && break
                j0 = j
            end
            cnt == 0 && return (det = zero(R), core = 0, computed = true)
            if cnt == 1
                hit = (i, j0); break
            end
        end
        if hit === nothing
            # 2. else column-wise — also detects an all-zero column
            for j in 1:n
                col_alive[j] || continue
                i0 = 0; cnt = 0
                for i in 1:n
                    (row_alive[i] && !iszero(M[i, j])) || continue
                    cnt += 1
                    cnt > 1 && break
                    i0 = i
                end
                cnt == 0 && return (det = zero(R), core = 0, computed = true)
                if cnt == 1
                    hit = (i0, j); break
                end
            end
        end
        hit === nothing && break               # dense core reached

        i, j = hit
        isodd(rank_of(i, row_alive) + rank_of(j, col_alive)) && (sgn = -sgn)
        acc = acc * M[i, j]
        row_alive[i] = false; col_alive[j] = false
        left -= 1
    end

    left == 0 && return (det = sgn * acc, core = 0, computed = true)

    left > maxcore && return (det = zero(R), core = left, computed = false)

    rows = [i for i in 1:n if row_alive[i]]
    cols = [j for j in 1:n if col_alive[j]]
    return (det = sgn * acc * soergel_det(M[rows, cols]), core = left, computed = true)
end

"""
    is_invertible_over_frac_ring(M::AbstractMatrix{SoergelPoly}) -> Bool

`true` iff `M` is invertible over `Frac(R)`, `R = QQ[α₁,α₂,α₃]`.  `R` is an
integral domain, so this is equivalent to `det(M) != 0` in `R` itself — there
is no need to actually construct `Frac(R)`.  Non-square matrices are never
invertible (`false`).

**The triangular test first** ([`is_permutation_triangular`](@ref)): if it
succeeds, `det(M)` is the product of the paired entries and hence `≠ 0` —
without computing a determinant at all.  Otherwise
[`soergel_det_peeled`](@ref) strips the sparse rows and columns and only the
dense core reaches the elimination.
"""
function is_invertible_over_frac_ring(M::AbstractMatrix{SoergelPoly};
                                      maxcore::Int = typemax(Int))
    n, m = size(M)
    n == m || return false
    is_permutation_triangular(M) && return true
    r = soergel_det_peeled(M; maxcore = maxcore)
    r.computed && return !iszero(r.det)
    return !iszero(soergel_det(M))
end
