# src/morphism/DLBasis.jl — the Double Leaves basis API.
#
# WHAT FOR. Light/Double Leaves are normally addressed via their 01-sequences:
# `light_leaves(x, z)` returns pairs `(e, morphism)`, `double_leaves(x, y)` a
# flat list of NamedTuples. Everything built on top of them — the circular
# counterpart, the basis writer, the flip decompositions, the structure
# constants — instead wants **the basis per degree**:
#
#     t = dl(x, y)      # DLTable
#     t[0]              # all Double Leaves of degree 0
#     t[0][1]           # the first one
#     t[0][1].e         # its 01-sequence on x
#
# `show(t)` prints the GRADED RANK TABLE — exactly what the Hecke ground truth
# `expected_count(x, y; degree = d)` predicts; the comparison is available as
# `expected_ranks`/`ranks_match`.
#
# ORDER WITHIN A DEGREE. `double_leaves` iterates over a `Set` of expressed
# elements; its iteration order is a hash order and thus not something a basis
# INDEX can rely on. So entries here are always sorted lexicographically by
# `(z, e, f)` — `t[i][j]` is a reproducible name for a basis element across
# sessions.
#
# CACHE. A per-session cache per table, keyed by the words; only what is requested
# is computed — no precomputation, no disk. `clear_dl_cache!()` clears it.

# ---- entries -----------------------------------------------------------------

"""
    LLEntry

A Light Leaf as a basis element: the 01-sequence `e` on the word, the expressed
element `z`, the degree (`defect`), and the morphism itself.
"""
struct LLEntry
    e::Vector{Int}
    z::Vector{Int}
    degree::Int
    morphism::MorphismGraph
end

"""
    DLEntry

A Double Leaf as a basis element: the 01-sequences `e` (on `x`) and `f` (on
`y`), the common expressed element `z`, the degree
`defect(x,e) + defect(y,f)`, and the morphism `x → y`.
"""
struct DLEntry
    e::Vector{Int}
    f::Vector{Int}
    z::Vector{Int}
    degree::Int
    morphism::MorphismGraph
end

Base.show(io::IO, en::LLEntry) =
    print(io, "LL(", _dlb_tok(en.e), " → ", _dlb_tok(en.z), ", degree ", en.degree, ")")
Base.show(io::IO, en::DLEntry) =
    print(io, "DL(", _dlb_tok(en.e), "|", _dlb_tok(en.f),
          " over ", _dlb_tok(en.z), ", degree ", en.degree, ")")
_dlb_tok(w::Vector{Int}) = isempty(w) ? "ε" : join(w)

# ---- tables -------------------------------------------------------------------

"""
    LLTable

The Light Leaves of a word, sorted by degree. `t[d]` is the vector of all
entries of degree `d` (empty if there are none), `t[d][j]` the `j`-th one.
`zword === nothing` means "over all reachable `z`".
"""
struct LLTable
    word::Vector{Int}
    zword::Union{Nothing,Vector{Int}}
    by_degree::Dict{Int,Vector{LLEntry}}
end

"""
    DLTable

The Double Leaves of a pair of words, sorted by degree — same access pattern
as `LLTable`.
"""
struct DLTable
    x::Vector{Int}
    y::Vector{Int}
    by_degree::Dict{Int,Vector{DLEntry}}
end

"""
    CircularDLEntry

A Double Leaf as a circular diagram — same data as [`DLEntry`](@ref), but with
`CircularMorphismGraph` instead of `MorphismGraph`. Unreduced; the reduction to circular
leaves is done by [`circular_dl_matrix`](@ref), which sits here
because `_DLB_TABLE` needs to know the table type.
"""
struct CircularDLEntry
    e::Vector{Int}
    f::Vector{Int}
    z::Vector{Int}
    degree::Int
    morphism::CircularMorphismGraph
end

"""
    CircularDLTable

The circular Double Leaves of a pair of words, sorted by degree — same access
pattern as [`DLTable`](@ref), and **index-parallel** to it: `circular_dl` and `dl`
sort by the same key `(z, e, f)`, so `circular_dl(x,y)[d][j]` corresponds to
`dl(x,y)[d][j]`. Whether that correspondence survives REDUCTION is answered by
[`is_circular_dl_basis`](@ref).
"""
struct CircularDLTable
    x::Vector{Int}
    y::Vector{Int}
    by_degree::Dict{Int,Vector{CircularDLEntry}}
end

Base.show(io::IO, en::CircularDLEntry) =
    print(io, "CircularDL(", _dlb_tok(en.e), "|", _dlb_tok(en.f),
          " over ", _dlb_tok(en.z), ", degree ", en.degree, ")")

const _DLB_TABLE = Union{LLTable,DLTable,CircularDLTable}

Base.getindex(t::_DLB_TABLE, d::Integer) = get(t.by_degree, Int(d), _dlb_empty(t))
_dlb_empty(::LLTable) = LLEntry[]
_dlb_empty(::DLTable) = DLEntry[]

"""
    dl_degrees(t) -> Vector{Int}

The occupied degrees, ascending.
"""
dl_degrees(t::_DLB_TABLE) = sort!(collect(keys(t.by_degree)))

"""
    dl_ranks(t) -> Vector{Pair{Int,Int}}

`degree => count`, ascending by degree — the graded rank table that `show`
prints.
"""
dl_ranks(t::_DLB_TABLE) = [d => length(t[d]) for d in dl_degrees(t)]

Base.length(t::_DLB_TABLE) = sum(length(v) for v in values(t.by_degree); init = 0)

"""
    dl_entries(t) -> Vector

All entries, by degree and then in table order — flat, for loops.
"""
dl_entries(t::_DLB_TABLE) = reduce(vcat, (t[d] for d in dl_degrees(t)); init = _dlb_empty(t))

# ---- construction --------------------------------------------------------------

# Deterministic order within a degree: (z, e, f) lexicographic. `double_leaves`
# iterates over a Set — without this sort, `t[i][j]` would not be a stable name.
_dlb_by_degree(entries_::Vector{LLEntry}) = _dlb_group(entries_, en -> (en.z, en.e))
_dlb_by_degree(entries_::Vector{DLEntry}) = _dlb_group(entries_, en -> (en.z, en.e, en.f))

function _dlb_group(entries_::Vector{T}, key) where {T}
    out = Dict{Int,Vector{T}}()
    for en in entries_
        push!(get!(out, en.degree, T[]), en)
    end
    for v in values(out)
        sort!(v; by = key)
    end
    return out
end

const _LL_CACHE = Dict{Tuple{Vector{Int},Union{Nothing,Vector{Int}}},LLTable}()
const _DL_CACHE = Dict{Tuple{Vector{Int},Vector{Int}},DLTable}()

"""
    clear_dl_cache!()

Clears the session cache of `ll`/`dl`.
"""
function clear_dl_cache!()
    empty!(_LL_CACHE)
    empty!(_DL_CACHE)
    empty!(_FDL_CACHE)          # defined below
    empty!(_FDL_MATRIX_CACHE)
    return nothing
end

"""
    ll(x)        -> LLTable
    ll(x, z)     -> LLTable

The Light Leaves on word `x`, sorted by degree: `ll(x)` over all expressed
elements, `ll(x, z)` only over `z`. Built via `light_leaves`; the result is
kept in the session cache.

```julia
t = ll([1,2,1])
t[0]            # the Light Leaves of degree 0
t[0][1].e       # their 01-sequence
```
"""
function ll(x::Vector{Int})
    get!(_LL_CACHE, (x, nothing)) do
        es = LLEntry[]
        for e in subexpressions(length(x))
            z = expressed_word(x, e)
            push!(es, LLEntry(e, z, defect(x, e), light_leaf(x, e)))
        end
        LLTable(x, nothing, _dlb_by_degree(es))
    end
end

function ll(x::Vector{Int}, z::Vector{Int})
    get!(_LL_CACHE, (x, z)) do
        es = [LLEntry(e, z, defect(x, e), m) for (e, m) in light_leaves(x, z)]
        LLTable(x, z, _dlb_by_degree(es))
    end
end

"""
    dl(x, y)  -> DLTable
    dl(x)     -> DLTable        # = dl(x, x)

The Double Leaves `x → y`, sorted by degree. Built via `double_leaves`; the
result is kept in the session cache.

```julia
t = dl([1,2,1])
t                       # graded rank table
t[0][1].e, t[0][1].f    # 01-sequences of the first degree-0 basis element
```
"""
function dl(x::Vector{Int}, y::Vector{Int})
    haskey(_DL_CACHE, (x, y)) && return _DL_CACHE[(x, y)]
    es = [DLEntry(d.e, d.f, d.z, d.degree, d.morphism) for d in double_leaves(x, y)]
    t = DLTable(x, y, _dlb_by_degree(es))
    _DL_CACHE[(x, y)] = t
    # Warn only on the FIRST build — the cache returns the same object afterward.
    return _dlb_warn_if_not_basis(t)
end

dl(x::Vector{Int}) = dl(x, x)

# ---- comparison against the Hecke ground truth ---------------------------------

"""
    expected_ranks(t::DLTable) -> Vector{Pair{Int,Int}}

The rank table predicted by the Hecke pairing table (`expected_count`), in the
same shape as `dl_ranks`. Throws if `(x, y)` is not in the table.
"""
function expected_ranks(t::DLTable)
    degs = bs_pairing(t.x, t.y)
    return [d => degs[d] for d in sort!(collect(keys(degs)))]
end

"""
    ranks_match(t::DLTable) -> Bool

`true` if the measured rank table matches the Hecke prediction degree by
degree. Defined for every word pair, since the pairing is computed rather than
looked up ([`bs_pairing`](@ref)).
"""
ranks_match(t::DLTable) = dl_ranks(t) == expected_ranks(t)

"""
    is_dl_basis(t::DLTable) -> Bool

Is the table actually a **basis**? The comparison of its graded ranks against
the Hecke ground truth.

This is the precondition for **any indexing**: if the graded ranks don't
match, `t[d][j]` is not a name for a basis element, just a slot in a list.

The pairing is COMPUTED for every word pair ([`bs_pairing`](@ref)), so the answer
is always `true` or `false`.
"""
is_dl_basis(t::DLTable) = dl_ranks(t) == expected_ranks(t)

# ---- warnings -------------------------------------------------------------------
#
# The principle for every table in this layer: one whose ranks don't match is
# NOT a basis, and any assignment built on it (index, basis change, structure
# constant) is meaningless. It is still built and returned — one has to look at
# it closely to find out why — but it says LOUDLY that it isn't one.
# `circular_dl` hangs off the same mechanism: there the question is not
# "do the ranks match" but "is the decomposition matrix DL ↦ Circular-DL
# (uni)triangular", and if it isn't there is no 1:1 correspondence
# `circular_dl(x)[i][j]` ↔ `dl(x)[i][j]`.
#
# Warned once per word pair (the cache only builds once anyway);
# `dl_warnings!(false)` silences it — for callers that log the finding themselves.

const _DLB_WARN = Ref(true)

"""
    dl_warnings!(on::Bool) -> Bool

Turns the basis warnings of `dl`/`circular_dl` on or off (returns the previous
state). Off only for callers that record the finding themselves.
"""
function dl_warnings!(on::Bool)
    old = _DLB_WARN[]
    _DLB_WARN[] = on
    return old
end

function _dlb_warn_if_not_basis(t::DLTable)
    _DLB_WARN[] || return t
    is_dl_basis(t) && return t
    @warn """
          dl($(_dlb_tok(t.x)), $(_dlb_tok(t.y))): the graded ranks do NOT \
          match the Hecke ground truth — this table is not a basis, \
          `t[d][j]` names no basis element. Warning can be disabled with \
          `dl_warnings!(false)`.""" measured = dl_ranks(t) expected = expected_ranks(t)
    return t
end

# ---- show -------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", t::LLTable)
    head = t.zword === nothing ? "LLTable $(_dlb_tok(t.word))" :
                                 "LLTable $(_dlb_tok(t.word)) → $(_dlb_tok(t.zword))"
    println(io, head, "  (", length(t), " Light Leaves)")
    _dlb_print_ranks(io, t)
end

function Base.show(io::IO, ::MIME"text/plain", t::DLTable)
    println(io, "DLTable $(_dlb_tok(t.x)) → $(_dlb_tok(t.y))  (", length(t), " Double Leaves)")
    _dlb_print_ranks(io, t)
    # The comparison against the Hecke ground truth is the actual point of the
    # table, so it is always printed — a silent nothing must not look like a pass.
    if is_dl_basis(t)
        println(io, "  Hecke: matches (basis)")
    else
        println(io, "  Hecke: MISMATCH — expected ",
                join(("$d:$n" for (d, n) in expected_ranks(t)), " "))
        println(io, "  ⚠ NOT a basis — t[d][j] names no basis element")
    end
end

function _dlb_print_ranks(io::IO, t::_DLB_TABLE)
    ds = dl_degrees(t)
    if isempty(ds)
        println(io, "  (empty)")
        return
    end
    println(io, "  degree ", join((lpad(string(d), 4) for d in ds)))
    println(io, "  count  ", join((lpad(string(length(t[d])), 4) for d in ds)))
end

Base.show(io::IO, t::LLTable) = print(io, "LLTable(", _dlb_tok(t.word), ", ", length(t), ")")
Base.show(io::IO, t::DLTable) =
    print(io, "DLTable(", _dlb_tok(t.x), " → ", _dlb_tok(t.y), ", ", length(t), ")")

# =============================================================================
# THE CIRCULAR COUNTERPART: `circular_dl` and the decomposition matrix
# =============================================================================
#
# WHAT THIS IS ABOUT. `dl(x, y)` is the Double Leaves
# basis. On the circular side sits the same list, just read as `CircularMorphismGraph`
# (`circular_double_leaves`) — that's `circular_dl(x, y)`, a pure reinterpretation: same
# entries, same degrees, same order.
#
# THE ACTUAL QUESTION is different: reducing a Double Leaf in the circular picture
# (`reduce_to_circular_leave`), it decomposes into CIRCULAR LEAVES,
#
#     circular(DL_k)  =  Σ_l  M_kl · CircularDL_l ,
#
# and the CONJECTURE is: `M` is (uni)triangular, so there is a 1:1
# correspondence DL ↔ Circular-DL, and `circular_dl(x)[d][j]` names the partner of
# `dl(x)[d][j]`.
#
# The two properties to check on `M` are that it is square (column count == row
# count) and invertible over Frac(R); `circular_dl_triangular` checks the stronger
# TRIANGULAR FORM.
#
# WARNING BEHAVIOR. For the other tables the basis
# question is cheap (rank comparison against the Hecke table). Here it costs a
# FULL reduction of every Double Leaf — milliseconds at length 2, minutes to
# hours at length 6. An `@warn` on build would force that cost on every
# `circular_dl` call. So: `circular_dl` builds silently, `show` says "triangular form
# NOT checked", and `is_circular_dl_basis(t)` computes it — THAT one warns, with a
# reason.

_dlb_empty(::CircularDLTable) = CircularDLEntry[]
_dlb_by_degree(entries_::Vector{CircularDLEntry}) = _dlb_group(entries_, en -> (en.z, en.e, en.f))

const _FDL_CACHE = Dict{Tuple{Vector{Int},Vector{Int}},CircularDLTable}()

"""
    circular_dl(x, y) -> CircularDLTable
    circular_dl(x)    -> CircularDLTable        # = circular_dl(x, x)

The circular counterpart to [`dl`](@ref): the same Double Leaves, read as
`CircularMorphismGraph` (via `circular_double_leaves`), sorted by degree and
index-parallel to `dl(x, y)`.

Construction is cheap and **checks nothing** — the question "is this a basis?"
costs a full reduction per Double Leaf here and lives in
[`is_circular_dl_basis`](@ref) instead.
"""
function circular_dl(x::Vector{Int}, y::Vector{Int})
    get!(_FDL_CACHE, (x, y)) do
        es = [CircularDLEntry(d.e, d.f, d.z, d.degree, d.morphism)
              for d in circular_double_leaves(x, y)]
        CircularDLTable(x, y, _dlb_by_degree(es))
    end
end

circular_dl(x::Vector{Int}) = circular_dl(x, x)

"""
    CircularDLDecomposition

The decomposition of the DL basis into circular leaves, `circular(DL_k) = Σ_l M_kl · CircularDL_l`:

* `rows` — the Double Leaves in table order (`dl_entries(dl(x,y))`),
  `row_degree` their degrees;
* `col_key` — the canonical keys of the circular leaves that occur, in order of
  first occurrence, `col_leaf` a representative diagram for each;
* `col_degree` — the smallest row degree at which the column occurs (a circular
  leaf has no 01-sequence of its own, this is the substitute);
* `M` — the coefficient matrix over `R`.

Same matrix the rank and triangularity checks build.
"""
struct CircularDLDecomposition
    x::Vector{Int}
    y::Vector{Int}
    rows::Vector{DLEntry}
    row_degree::Vector{Int}
    col_key::Vector{Any}
    col_leaf::Vector{Any}
    col_degree::Vector{Int}
    M::Matrix{SoergelPoly}
end

const _FDL_MATRIX_CACHE = Dict{Tuple{Vector{Int},Vector{Int}},CircularDLDecomposition}()

"""
    circular_dl_matrix(x, y) -> CircularDLDecomposition

Computes the decomposition matrix: reduce each Double Leaf `x → y` as a circular
diagram with `reduce_to_circular_leave` and merge terms via `circular_canonical_key`.
**Expensive** (a full reduction per Double Leaf); the result is kept in the
session cache.
"""
function circular_dl_matrix(x::Vector{Int}, y::Vector{Int})
    get!(_FDL_MATRIX_CACHE, (x, y)) do
        rows = dl_entries(dl(x, y))
        combos = Any[]
        for en in rows
            fm = circular_morphism(en.morphism)
            fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
            push!(combos, reduce_to_circular_leave(fdm))
        end

        col_key = Any[]; col_leaf = Any[]; index = Dict{Any,Int}()
        for c in combos, (d, _) in pairs_of(c)
            k = circular_canonical_key(d)
            haskey(index, k) && continue
            index[k] = length(col_key) + 1
            push!(col_key, k); push!(col_leaf, d)
        end

        n, m = length(rows), length(col_key)
        M = fill(zero(SoergelPoly), n, m)
        for (i, c) in enumerate(combos), (d, coeff) in pairs_of(c)
            M[i, index[circular_canonical_key(d)]] = coeff
        end

        rdeg = [en.degree for en in rows]
        cdeg = [minimum(rdeg[i] for i in 1:n if !iszero(M[i, j])) for j in 1:m]
        CircularDLDecomposition(x, y, rows, rdeg, col_key, col_leaf, cdeg, M)
    end
end

circular_dl_matrix(x::Vector{Int}) = circular_dl_matrix(x, x)
circular_dl_matrix(t::CircularDLTable) = circular_dl_matrix(t.x, t.y)

"""
    circular_dl_triangular(dec::CircularDLDecomposition) -> NamedTuple

Checks the triangularity conjecture on a finished decomposition matrix and
returns the REASON, not just a yes/no:

* `square` — row count == column count? Otherwise no 1:1 correspondence can
  exist, and `excess_by_degree` says how many circular leaves are surplus in which
  degree;
* `triangular` — sorting columns by `(col_degree, first occurrence)`, is
  `M_kl = 0` whenever the column has a SMALLER degree than the row?
* `block_diagonal` — nothing above either, i.e. the matrix decomposes into
  blocks per degree (the sharper, measured statement);
* `unitriangular` — `triangular`, and there is a column permutation with all
  `1`s on the diagonal;
* `partner` — that very correspondence: `partner[k]` is the column index of
  the circular-leaf partner of `rows[k]` (`0` where there is none);
* `perm` — the corresponding column permutation (empty if none exists).

⚠ THE CORRESPONDENCE IS MATCHED, NOT GUESSED. First-occurrence order does not
always give the correct diagonal — for `2132 → 2132`, two degree-2 rows are
merely swapped there, even though the correspondence exists. A perfect
matching on the entries `== 1` is therefore searched per degree block
(Kuhn's algorithm; the blocks are small).
"""
function circular_dl_triangular(dec::CircularDLDecomposition)
    n, m = length(dec.rows), length(dec.col_key)
    excess = Dict{Int,Int}()
    for d in dec.col_degree; excess[d] = get(excess, d, 0) + 1; end
    for d in dec.row_degree; excess[d] = get(excess, d, 0) - 1; end
    excess_by_degree = [d => excess[d] for d in sort!(collect(keys(excess))) if excess[d] != 0]

    if n != m
        return (; square = false, triangular = false, block_diagonal = false,
                unitriangular = false, partner = zeros(Int, n), perm = Int[],
                excess_by_degree)
    end

    lower_free = all(iszero(dec.M[i, j])
                     for i in 1:n, j in 1:m if dec.col_degree[j] < dec.row_degree[i])
    upper_free = all(iszero(dec.M[i, j])
                     for i in 1:n, j in 1:m if dec.col_degree[j] > dec.row_degree[i])

    # A perfect matching on the one-entries, per degree block.
    partner = zeros(Int, n)
    colof = zeros(Int, m)        # column -> assigned row
    complete = true
    for d in sort!(unique(dec.row_degree))
        rows_d = [i for i in 1:n if dec.row_degree[i] == d]
        cols_d = [j for j in 1:m if dec.col_degree[j] == d]
        length(rows_d) == length(cols_d) || (complete = false)
        for i in rows_d
            seen = Set{Int}()
            _dlb_augment!(i, cols_d, seen, partner, colof, dec) || (complete = false)
        end
    end

    perm = complete ? [partner[i] for i in 1:n] : Int[]
    uni = complete && (n == 0 || all(isone(dec.M[i, partner[i]]) for i in 1:n))
    return (; square = true, triangular = lower_free,
            block_diagonal = lower_free && upper_free,
            unitriangular = lower_free && uni, partner, perm, excess_by_degree)
end

# Kuhn step: try to assign row `i` to a still-free (or freeable) column from
# `cols` with entry `1`.
function _dlb_augment!(i::Int, cols::Vector{Int}, seen::Set{Int},
                       partner::Vector{Int}, colof::Vector{Int},
                       dec::CircularDLDecomposition)
    for j in cols
        (isone(dec.M[i, j]) && !(j in seen)) || continue
        push!(seen, j)
        if colof[j] == 0 || _dlb_augment!(colof[j], cols, seen, partner, colof, dec)
            colof[j] = i
            partner[i] = j
            return true
        end
    end
    return false
end

"""
    is_circular_dl_basis(t::CircularDLTable) -> Bool or missing
    is_circular_dl_basis(x, y)          -> Bool or missing

Does the correspondence `circular_dl(x,y)[d][j]` ↔ `dl(x,y)[d][j]` hold? Checked via
the triangularity conjecture on the decomposition matrix
([`circular_dl_matrix`](@ref)): `true` if it is square and unitriangular, `false`
otherwise — with an `@warn` naming the REASON (how many circular leaves are surplus
in which degree, or why triangular form fails). `missing` if the reduction
throws, i.e. the question cannot be answered at all.

**Expensive** — the matrix is a full reduction per Double Leaf. To silence the
warning, use `dl_warnings!(false)`.
"""
function is_circular_dl_basis(t::CircularDLTable)
    dec = try
        circular_dl_matrix(t.x, t.y)
    catch err
        _DLB_WARN[] && @warn """
              circular_dl($(_dlb_tok(t.x)), $(_dlb_tok(t.y))): the reduction to \
              circular leaves THROWS — whether the correspondence holds is NOT \
              decided.""" error = err
        return missing
    end
    r = circular_dl_triangular(dec)
    r.unitriangular && return true
    _DLB_WARN[] && @warn """
          circular_dl($(_dlb_tok(t.x)), $(_dlb_tok(t.y))): the decomposition matrix \
          circular(DL) = Σ M · CircularDL is NOT unitriangular — there is no 1:1 \
          correspondence, `circular_dl(…)[d][j]` is NOT the partner of \
          `dl(…)[d][j]`. Warning can be disabled with `dl_warnings!(false)`.""" reason =
              (r.square ? "square, but triangular form or diagonal violated" :
                          "too many circular leaves") excess_by_degree =
              r.excess_by_degree triangular_form = r.triangular
    return false
end

is_circular_dl_basis(x::Vector{Int}, y::Vector{Int}) = is_circular_dl_basis(circular_dl(x, y))

function Base.show(io::IO, ::MIME"text/plain", t::CircularDLTable)
    println(io, "CircularDLTable $(_dlb_tok(t.x)) → $(_dlb_tok(t.y))  (", length(t),
                " Circular Double Leaves)")
    _dlb_print_ranks(io, t)
    dec = get(_FDL_MATRIX_CACHE, (t.x, t.y), nothing)
    if dec === nothing
        println(io, "  triangular form: NOT checked (costs a full reduction per ",
                    "Double Leaf) — `is_circular_dl_basis(t)` computes it")
    else
        r = circular_dl_triangular(dec)
        println(io, "  decomposition: ", length(dec.rows), " DLs → ", length(dec.col_key),
                    " circular leaves; unitriangular: ", r.unitriangular,
                    r.square ? "" : "  ⚠ excess $(r.excess_by_degree)")
    end
end

Base.show(io::IO, t::CircularDLTable) =
    print(io, "CircularDLTable(", _dlb_tok(t.x), " → ", _dlb_tok(t.y), ", ", length(t), ")")

Base.show(io::IO, dec::CircularDLDecomposition) =
    print(io, "CircularDLDecomposition(", _dlb_tok(dec.x), " → ", _dlb_tok(dec.y), ", ",
          length(dec.rows), "×", length(dec.col_key), ")")
