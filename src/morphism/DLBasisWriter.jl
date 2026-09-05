# morphism/DLBasisWriter.jl — writing an arbitrary morphism in the DL basis.
# It sits after DLBasis.jl in the include order because its signatures need
# DecoratedMorphism and CircularDecoratedMorphism/reduce_to_circular_leave.

# =============================================================================
# THE BASIS WRITER: `in_circular_dl_basis`
# =============================================================================
#
# WHAT THIS IS ABOUT. Write an arbitrary map `x → y`
# (DiagramComboR resp. CircularComboR) as a linear combination of the respective
# basis. The mechanics:
#
# * PLAIN: `reduce_to_dl` on every term, then match coefficients via the
#   canonical keys of the basis elements (the basis DLs are their own normal
#   forms — checked, not assumed).
# * CIRCULAR: `reduce_to_circular_leave` gives circular-LEAF coordinates; the basis elements
#   decompose as rows of the decomposition matrix `M`. So we solve for `c`
#   with `cᵀ·M = v` — exact when unitriangular (Gauss over `R`, pivots are
#   constants).
#
# What cannot be matched ends up in `rest`. The result is still returned, but
# `exact = false` plus an `@warn` say it is NOT a basis representation.
#
# A basis DL is expected to be its own normal form (1 term, coefficient 1, trivial
# labels); where that fails, the element is NOT included in the matching and a
# single warning is issued.
#
"""
    DLCombination

A map as a linear combination of the DL resp. circular-DL basis:

* `coeffs` — `(degree, j) => coefficient ∈ R`; `(degree, j)` names
  `dl(x,y)[degree][j]` (resp. its circular partner);
* `rest` — the terms that could NOT be matched to any basis element
  (`DiagramComboR`/`CircularComboR`; empty in the good case);
* `exact` — `true` if `rest` is empty and the matching was complete.

`basis` is `:dl` or `:circular_dl`.
"""
struct DLCombination
    x::Vector{Int}
    y::Vector{Int}
    basis::Symbol
    coeffs::Dict{Tuple{Int,Int},SoergelPoly}
    rest::Any
    exact::Bool
end

function Base.show(io::IO, ::MIME"text/plain", c::DLCombination)
    println(io, "DLCombination (", c.basis, ") $(_dlb_tok(c.x)) → $(_dlb_tok(c.y))")
    ks = sort!(collect(keys(c.coeffs)))
    isempty(ks) && println(io, "  0 (no basis contributions)")
    for k in ks
        println(io, "  [", k[1], "][", k[2], "]  ·  ", c.coeffs[k])
    end
    nrest = c.rest === nothing ? 0 : length(c.rest)
    nrest == 0 || println(io, "  ⚠ + ", nrest, " unmatched term(s) in `rest`")
    c.exact || println(io, "  ⚠ NOT a basis representation (exact = false)")
end

Base.show(io::IO, c::DLCombination) =
    print(io, "DLCombination(", c.basis, ", ", _dlb_tok(c.x), " → ", _dlb_tok(c.y),
          ", ", length(c.coeffs), " contributions", c.exact ? "" : ", INEXACT", ")")

# The writer is `in_circular_dl_basis` below; the four functions further down
# (`hflip_dl_matrix`, `vflip_dl_matrix`, `compose_in_basis`,
# `structure_constants`) all go through it. Matching runs on the circular side,
# where the bigon-zero rule is available as C5 and C1.

# ---- leaf rotation: normalize cut1 to 0 ---------------------------------------
#
# The basis elements (`double_leaf`) have `cut1 = 0` (bottom = leaves
# 1..|x|). A general morphism — e.g. `flip(DL)` with swapped cuts — numbers
# its leaves differently; the canonical keys then NEVER match even though it
# is the same morphism. So leaves are rotated by `cut1` before reduction
# (the same principle applies on the circular side). The boundary word
# rotates along so leaf and edge colours stay consistent (like `_mirror_graph`,
# diagram/Join.jl).

function _dlb_rot_leaves(g::WordGraph, r::Int)
    n = length(g.word)
    r = mod(r, n)
    r == 0 && return g
    cols = [leaf_colour(g, k) for k in 1:n]
    neww = CircularWord([cols[mod1(k + r, n)] for k in 1:n])
    rp(p) = p isa Leaf ? Leaf(mod1(p.k - r, n)) : p
    return WordGraph(neww, g.nodes, Edge[Edge(e.colour, rp(e.a), rp(e.b)) for e in g.edges])
end

function _dlb_rot_leaves(g::CircularGraph, r::Int)
    n = length(letters(g.word))
    r = mod(r, n)
    r == 0 && return g
    cols = collect(letters(g.word))
    neww = CircularWord([cols[mod1(k + r, n)] for k in 1:n])
    rp(p) = p isa Leaf ? Leaf(mod1(p.k - r, n)) : p
    return CircularGraph(neww, g.nodes, Edge[Edge(e.colour, rp(e.a), rp(e.b)) for e in g.edges])
end

# Rotation amount of a morphism (sentinel "empty bottom" ⇒ 0).
_dlb_cut_rot(cut1::Int, n::Int) = (n == 0 || cut1 < 0) ? 0 : mod(cut1, n)

# ---- circular: solving against the decomposition matrix ----------------------------

# Constant term of a polynomial (`nothing` if not constant) — for the pivots.
_dlb_const(p::SoergelPoly) = constant_term(p)

"""
    _dlb_solve_left(M, v) -> Union{Nothing, Vector{SoergelPoly}}

Solves `cᵀ·M = v` exactly over `R` (Gauss on `Mᵀ`); pivots must be constants
(the diagonal blocks of the decomposition matrices carry constant
coefficients). `nothing` if no constant pivot is found or the system has no
unique solution.
"""
function _dlb_solve_left(M::Matrix{SoergelPoly}, v::Vector{SoergelPoly})
    n, m = size(M)
    (m == length(v) && n == m) || return nothing
    # augmented matrix of Mᵀ·c = v
    A = Matrix{SoergelPoly}(undef, m, n + 1)
    for i in 1:m, j in 1:n
        A[i, j] = M[j, i]
    end
    for i in 1:m
        A[i, n + 1] = v[i]
    end
    piv = zeros(Int, n)
    row = 1
    for col in 1:n
        # find a constant, nonzero pivot
        k = findfirst(i -> begin c = _dlb_const(A[i, col]); c !== nothing && !iszero(c) end,
                      row:m)
        k === nothing && return nothing
        k = row + k - 1
        A[row, :], A[k, :] = A[k, :], A[row, :]
        pc = _dlb_const(A[row, col])
        inv_p = soergel_monomial((0, 0, 0), 1 // pc)
        for j in col:(n + 1)
            A[row, j] = inv_p * A[row, j]
        end
        for i in 1:m
            i == row && continue
            f = A[i, col]
            iszero(f) && continue
            for j in col:(n + 1)
                A[i, j] = A[i, j] - f * A[row, j]
            end
        end
        piv[col] = row
        row += 1
    end
    # consistency: remaining rows must be 0 = 0
    for i in row:m
        iszero(A[i, n + 1]) || return nothing
    end
    return [A[piv[j], n + 1] for j in 1:n]
end

"""
    in_circular_dl_basis(c::CircularComboR, x, y)        -> DLCombination
    in_circular_dl_basis(fdm::CircularDecoratedMorphism) -> DLCombination
    in_circular_dl_basis(m::MorphismGraph)          -> DLCombination

Writes a map `x → y` in the CIRCULAR-DL basis: bring it to circular-leaf coordinates
with `reduce_to_circular_leave` and solve against the decomposition matrix `M` from
the decomposition matrix (`cᵀ·M = v`, exact when unitriangular). **Expensive** — needs
[`circular_dl_matrix`](@ref).

Terms whose circular leaf doesn't occur in `M` at all end up in `rest` (the "stray"
finding); then — and when solving fails (matrix not
square/unique) — `exact = false`, with a warning.
"""
function in_circular_dl_basis(c::CircularComboR, x::Vector{Int}, y::Vector{Int})
    dec = circular_dl_matrix(x, y)
    index = Dict{Any,Int}(k => i for (i, k) in enumerate(dec.col_key))
    v = fill(zero(SoergelPoly), length(dec.col_key))
    rest = CircularComboR()
    for (d0, coeff0) in pairs_of(c)
        # Same as in `in_dl_basis`: the outer label is a scalar FACTOR of the
        # term, not a cell of the diagram, so it belongs in the coefficient
        # before matching against the columns of `M` (which carry no labels).
        # Without this the barbell case lost its `α₁` silently — see the box in
        # `in_dl_basis`.
        d, coeff = isone(d0.outer_label) ? (d0, coeff0) :
                   (CircularDecorated(d0.graph, d0.region_labels), coeff0 * d0.outer_label)
        k = circular_canonical_key(d)
        if haskey(index, k)
            v[index[k]] = v[index[k]] + coeff
        else
            rest = rest + coeff * CircularComboR(d)
        end
    end
    sol = _dlb_solve_left(dec.M, v)
    coeffs = Dict{Tuple{Int,Int},SoergelPoly}()
    if sol !== nothing
        # row i of the matrix -> (degree, j) of the table (same order as dl_entries)
        seen = Dict{Int,Int}()
        for (i, en) in enumerate(dec.rows)
            j = seen[en.degree] = get(seen, en.degree, 0) + 1
            iszero(sol[i]) || (coeffs[(en.degree, j)] = sol[i])
        end
    end
    exact = sol !== nothing && isempty(rest)
    exact || (_DLB_WARN[] && @warn """
          in_circular_dl_basis($(_dlb_tok(x)), $(_dlb_tok(y))): \
          $(sol === nothing ? "the system cᵀ·M = v has no unique solution" :
            "$(length(rest)) circular-leaf term(s) do not occur in M (`rest`)") \
          — not a basis representation. Warning can be disabled with \
          `dl_warnings!(false)`.""")
    return DLCombination(x, y, :circular_dl, coeffs, rest, exact)
end

function in_circular_dl_basis(fdm::CircularDecoratedMorphism)
    x, y = bottom(fdm.m), top(fdm.m)
    r = _dlb_cut_rot(fdm.m.cut1, length(letters(fdm.m.graph.word)))
    if r != 0
        # same normalization as in `in_dl_basis`: rotate cut1 to 0 so the keys
        # are comparable to the `double_leaf` basis elements.
        all(isone, fdm.region_labels) || error(
            "in_circular_dl_basis: leaf rotation (cut1 = $(fdm.m.cut1)) with " *
            "nontrivial region labels is not supported")
        g2 = _dlb_rot_leaves(fdm.m.graph, r)
        fm2 = CircularMorphismGraph(g2, 0, length(x))
        fdm = CircularDecoratedMorphism(fm2, fill(one(SoergelPoly), region_count(g2)))
    end
    return in_circular_dl_basis(reduce_to_circular_leave(fdm), x, y)
end

function in_circular_dl_basis(m::MorphismGraph)
    fm = circular_morphism(m)
    return in_circular_dl_basis(
        CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph))))
end

# ---- basis change DL ↔ Circular -----------------------------------------------------

"""
    dl_to_circular_matrix(x, y) -> Matrix{SoergelPoly}

The basis-change matrix DL ↦ Circular-Leaves: `circular(DL_k) = Σ_l M_kl · CircularDL_l` —
exactly `circular_dl_matrix(x, y).M` (rows in `dl_entries` order, columns in the
decomposition's `col_key` order).
"""
dl_to_circular_matrix(x::Vector{Int}, y::Vector{Int}) = circular_dl_matrix(x, y).M

"""
    circular_to_dl_matrix(x, y) -> Union{Nothing, Matrix{SoergelPoly}}

The way back: the exact inverse `M⁻¹` of the decomposition matrix (column-l
circular leaf as a combination of the `circular(DL_k)`). Exists when unitriangular
(`circular_dl_matrix`); otherwise — matrix not square or not invertible over `R` —
`nothing`, with a warning.
"""
function circular_to_dl_matrix(x::Vector{Int}, y::Vector{Int})
    dec = circular_dl_matrix(x, y)
    n, m = size(dec.M)
    inv_ok = n == m
    cols = Vector{Vector{SoergelPoly}}()
    if inv_ok
        for l in 1:m
            e = [i == l ? one(SoergelPoly) : zero(SoergelPoly) for i in 1:m]
            sol = _dlb_solve_left(dec.M, e)     # cᵀ·M = e_l  ⇒  column l of M⁻¹ᵀ
            sol === nothing && (inv_ok = false; break)
            push!(cols, sol)
        end
    end
    if !inv_ok
        _DLB_WARN[] && @warn """
              circular_to_dl_matrix($(_dlb_tok(x)), $(_dlb_tok(y))): the \
              decomposition matrix is not invertible (over R) — no way back \
              Circular ↦ DL. Warning can be disabled with `dl_warnings!(false)`."""
        return nothing
    end
    Minv = Matrix{SoergelPoly}(undef, n, n)
    for l in 1:n, k in 1:n
        Minv[k, l] = cols[l][k]
    end
    return Minv
end

# =============================================================================
# THE HFLIP DECOMPOSITION: the cellular anti-involution
# =============================================================================
#
# WHAT THIS IS ABOUT. `hflip(DL_{(e,f)}: x→y)` — axis =
# HORIZONTAL, swaps bottom↔top; in the plain stack this exact map is called
# `flip(::MorphismGraph)` (naming note in diagram/Join.jl: the old
# `MorphismGraph` names are swapped relative to the circular convention) — is a map
# `y → x`; its decomposition in the DL basis of `(y, x)` as a matrix.
# EXPECTATION (cell structure): a pure swap `(e, f) ↦ (f, e)` with coefficient
# 1 — hflip should be the cellular ANTI-INVOLUTION (Elias–Williamson).
#
# `hflip_dl_matrix` computes with `flip` = bottom↔top; `vflip_dl_matrix` below is
# the true left-right mirror.
#
# The decomposition needs the cut1 normalization above; without it almost every
# case falsely looks like a "stray".

"""
    hflip_dl_matrix(x, y) -> NamedTuple

The decomposition of `flip(DL_k)` (bottom↔top, map `y → x`) in the DL basis of
`(y, x)`:

* `M` — the matrix over `R`; row `k` = `dl_entries(dl(x,y))[k]`, column `l` =
  `dl_entries(dl(y,x))[l]`;
* `pure_swap` — `true` if EVERY row carries exactly the entry 1 on the
  `(e,f) ↦ (f,e)` partner (the cellular anti-involution);
* `partner` — the column index of this partner per row (`0` where there is
  none);
* `exact` — `false` as soon as some `flip(DL)` did not lie exactly in the
  basis (then `M` is incomplete and `pure_swap == false`).

The decomposition runs via [`in_circular_dl_basis`](@ref) (expensive: it needs
`circular_dl_matrix(y, x)`).
"""
function hflip_dl_matrix(x::Vector{Int}, y::Vector{Int})
    rows = dl_entries(dl(x, y))
    cols = dl_entries(dl(y, x))
    # (degree, j) -> flat column index, and the (f,e) partner per row
    flat = Dict{Tuple{Int,Int},Int}()
    seen = Dict{Int,Int}()
    for (l, en) in enumerate(cols)
        j = seen[en.degree] = get(seen, en.degree, 0) + 1
        flat[(en.degree, j)] = l
    end
    partner = [something(findfirst(c -> c.e == en.f && c.f == en.e && c.z == en.z, cols), 0)
               for en in rows]

    n, m = length(rows), length(cols)
    M = fill(zero(SoergelPoly), n, m)
    exact = true
    for (k, en) in enumerate(rows)
        c = in_circular_dl_basis(flip(en.morphism))
        exact &= c.exact
        for (key, coeff) in c.coeffs
            haskey(flat, key) || (exact = false; continue)
            M[k, flat[key]] = coeff
        end
    end
    pure = exact && all(partner[k] != 0 &&
                        all(M[k, l] == (l == partner[k] ? one(SoergelPoly) : zero(SoergelPoly))
                            for l in 1:m)
                        for k in 1:n)
    return (; M, pure_swap = pure, partner, exact)
end

hflip_dl_matrix(x::Vector{Int}; kwargs...) = hflip_dl_matrix(x, x; kwargs...)

# =============================================================================
# THE TRUE VFLIP DECOMPOSITION: the left-right mirror
# =============================================================================
#
# `vflip` = axis VERTICAL: bottom stays bottom, top stays top, only the letter
# order reverses. On a `MorphismGraph` that map is called `hflip` (diagram/Join.jl,
# where the axis names are the other way round). `vflip(DL: x→y)` is therefore a map
# `reverse(x) → reverse(y)`, and its decomposition in the DL basis of
# `(reverse(x), reverse(y))` is what is computed here. It is NOT a bijection on the
# DLs — a DL generally maps to a LINEAR COMBINATION — so there is no `pure_swap`;
# `permutation` says whether the decomposition happens to be a permutation matrix
# anyway.

"""
    vflip_dl_matrix(x, y) -> NamedTuple

The decomposition of the left-right mirror `vflip(DL_k)` (map
`reverse(x) → reverse(y)`; `hflip(::MorphismGraph)` in the plain stack) in the
DL basis of `(reverse(x), reverse(y))`:

* `M` — the matrix over `R`; row `k` = `dl_entries(dl(x,y))[k]`, column `l` =
  `dl_entries(dl(reverse(x), reverse(y)))[l]`;
* `permutation` — `true` if every row carries exactly ONE entry 1 (a bijection on
  the DLs, which is generally NOT the case);
* `exact` — `false` as soon as some `vflip(DL)` did not lie exactly in the
  basis (then `M` is incomplete).

The decomposition runs via [`in_circular_dl_basis`](@ref) (expensive).
"""
function vflip_dl_matrix(x::Vector{Int}, y::Vector{Int})
    rows = dl_entries(dl(x, y))
    cols = dl_entries(dl(reverse(x), reverse(y)))
    flat = Dict{Tuple{Int,Int},Int}()
    seen = Dict{Int,Int}()
    for (l, en) in enumerate(cols)
        j = seen[en.degree] = get(seen, en.degree, 0) + 1
        flat[(en.degree, j)] = l
    end
    n, m = length(rows), length(cols)
    M = fill(zero(SoergelPoly), n, m)
    exact = true
    for (k, en) in enumerate(rows)
        v = hflip(en.morphism)          # plain name: hflip = left↔right
        c = in_circular_dl_basis(v)
        exact &= c.exact
        for (key, coeff) in c.coeffs
            haskey(flat, key) || (exact = false; continue)
            M[k, flat[key]] = coeff
        end
    end
    perm = exact && all(count(l -> !iszero(M[k, l]), 1:m) == 1 &&
                        all(iszero(M[k, l]) || isone(M[k, l]) for l in 1:m)
                        for k in 1:n)
    return (; M, permutation = perm, exact)
end

vflip_dl_matrix(x::Vector{Int}; kwargs...) = vflip_dl_matrix(x, x; kwargs...)

# =============================================================================
# THE STRUCTURE CONSTANTS OF COMPOSITION
# =============================================================================
#
# WHAT THIS IS ABOUT. For basis elements `f: x→y`,
# `g: y→z` the decomposition `g∘f = Σ c^{(m,n)} · DL^{x→z}_{(m,n)}` — the
# STRUCTURE CONSTANTS of the category in the DL basis (basis + anti-involution
# (0.4) + triangular form (0.2) together = "cellular basis", Elias–Williamson).
# Chapter 3 (the pairing) is the special case `x = z` with scalar extraction.
#
# ⚠ A composition can reduce to a circular leaf outside the basis. Such terms
# are not hidden: they end up in `rest`, `exact = false`, with a warning.

"""
    compose_in_basis(f::MorphismGraph, g::MorphismGraph) -> DLCombination

The decomposition of `g∘f` in the DL basis: `f: x → y` and `g: y → z` (order
as in [`compose`](@ref): **f is applied first**), result over `(x, z)`, via
[`in_circular_dl_basis`](@ref) (expensive).

Throws if `top(f) != bottom(g)` (not composable).
"""
function compose_in_basis(f::MorphismGraph, g::MorphismGraph)
    top(f) == bottom(g) || error(
        "compose_in_basis: top(f) = $(_dlb_tok(top(f))) doesn't match " *
        "bottom(g) = $(_dlb_tok(bottom(g)))")
    h = compose(f, g)                       # g∘f — f is applied first
    return in_circular_dl_basis(h)
end

compose_in_basis(f::DLEntry, g::DLEntry; kwargs...) =
    compose_in_basis(f.morphism, g.morphism; kwargs...)

"""
    structure_constants(x, y, z) -> NamedTuple

All structure constants `g∘f = Σ c · DL^{x→z}` for `f ∈ dl(x,y)`,
`g ∈ dl(y,z)`:

* `table` — `Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}}, DLCombination}`, key =
  (`(degree_f, j_f)`, `(degree_g, j_g)`);
* `exact` — `true` if ALL compositions lay exactly in the basis. Where not, the
  respective `DLCombination` (rest/exact) says what is missing, visibly;
* `thrown` — the pairs whose reduction THROWS (`apply_d4` handles only
  boundary-leaf strand ends); missing from `table`, and `exact` is then `false`.

**Expensive** (|dl(x,y)|·|dl(y,z)| reductions) — meant for small words.
"""
function structure_constants(x::Vector{Int}, y::Vector{Int}, z::Vector{Int})
    tf, tg = dl(x, y), dl(y, z)
    table = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},DLCombination}()
    thrown = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
    exact = true
    for df in dl_degrees(tf), (jf, enf) in enumerate(tf[df])
        for dg in dl_degrees(tg), (jg, eng) in enumerate(tg[dg])
            c = try
                compose_in_basis(enf.morphism, eng.morphism)
            catch
                push!(thrown, ((df, jf), (dg, jg)))
                exact = false
                continue
            end
            exact &= c.exact
            table[((df, jf), (dg, jg))] = c
        end
    end
    isempty(thrown) || (_DLB_WARN[] && @warn """
          structure_constants($(_dlb_tok(x)), $(_dlb_tok(y)), $(_dlb_tok(z))): \
          $(length(thrown)) composition(s) throw during reduction \
          (apply_d4 limit) and are missing from the table.""" pairs = thrown)
    return (; table, exact, thrown)
end
