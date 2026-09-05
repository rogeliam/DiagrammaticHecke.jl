# dlbasis.jl — the DL BASIS API, stage 0.1 (src/morphism/DLBasis.jl)
#
# Three things are to be shown:
#   1. the table contains EXACTLY the same leaves as `double_leaves`,
#   2. `t[d][j]` is a STABLE name (deterministic order),
#   3. the graded rank table agrees with the Hecke ground truth.

@testset "DL basis API stage 0.1 (ll / dl, graded rank tables)" begin
    x = [1, 2, 1]

    # ---- 1. same content as `double_leaves` ------------------------------------
    t = dl(x)
    @test t isa DLTable
    @test t.x == x && t.y == x
    @test length(t) == length(double_leaves(x, x))          # 12, see the morphism test
    @test dl(x) === dl(x, x)                                # dl(x) = dl(x,x), from the cache

    # degree by degree the same count as `double_leaves`'s degree filter
    for d in dl_degrees(t)
        @test length(t[d]) == length(double_leaves(x, x; degree = d))
    end
    @test dl_ranks(t) == [0 => 2, 2 => 5, 4 => 4, 6 => 1]   # by hand, cf. test/pairing.jl
    @test length(t[99]) == 0                                # unpopulated degree: empty, not an error

    # every entry carries its 01-sequences along (the way back to `double_leaves`)
    en = t[0][1]
    @test en isa DLEntry
    @test en.degree == 0
    @test en.z == expressed_word(x, en.e) == expressed_word(x, en.f)
    @test en.degree == defect(x, en.e) + defect(x, en.f)
    @test bottom(en.morphism) == x && top(en.morphism) == x

    # dl_entries: flat, all of them, by degree
    @test length(dl_entries(t)) == length(t)
    @test [e.degree for e in dl_entries(t)] == sort([e.degree for e in dl_entries(t)])

    # ---- 2. the order is deterministic -------------------------------------
    # `double_leaves` runs over a Set; without the sorting in `_dlb_by_degree`,
    # `t[d][j]` would not be a reproducible name. Freshly built must give the same as
    # from the cache.
    clear_dl_cache!()
    t2 = dl(x)
    @test [(e.z, e.e, e.f) for e in dl_entries(t2)] == [(e.z, e.e, e.f) for e in dl_entries(t)]
    # and within a degree it is ascending
    @test issorted([(e.z, e.e, e.f) for e in t[2]])

    # ---- 3. comparison with the Hecke ground truth -------------------------
    @test expected_ranks(t) == [0 => 2, 2 => 5, 4 => 4, 6 => 1]
    @test ranks_match(t)
    # also for an unequal pair from the table
    t3 = dl([1, 2, 1], [1, 3, 2])
    @test ranks_match(t3)
    @test sum(n for (_, n) in dl_ranks(t3)) == expected_count([1,2,1], [1,3,2])
    # EVERY word pair is covered, not just the 24×24 shortlex representatives:
    # the pairing is computed, not looked up. A non-shortlex reduced expression of
    # the same element (212 instead of 121) is therefore computed too, and gives the
    # same ranks as 121 — the pairing depends on the ELEMENT, not on the word.
    @test expected_ranks(dl([2,1,2], [2,1,2])) == expected_ranks(dl([1,2,1], [1,2,1]))
    @test ranks_match(dl([2,1,2], [2,1,2]))

    # ---- ll: the same access form -------------------------------------------
    lt = ll(x, [1])
    @test lt isa LLTable
    @test lt.zword == [1]
    @test length(lt) == length(light_leaves(x, [1]))         # 2
    @test dl_ranks(lt) == [0 => 1, 2 => 1]                   # cf. the morphism test
    @test Set(e.e for e in dl_entries(lt)) == Set([[1,0,0], [0,0,1]])
    @test all(e.z == [1] for e in dl_entries(lt))

    lta = ll(x)                                              # over ALL z
    @test lta.zword === nothing
    @test length(lta) == length(subexpressions(length(x)))   # 8, one per 01-sequence
    @test sum(length(lta[d]) for d in dl_degrees(lta)) == 8

    # `show` prints the rank table (not the leaves) — only the shape is checked
    s = sprint(show, MIME"text/plain"(), t)
    @test occursin("DLTable", s) && occursin("degree", s) && occursin("count", s)
    @test occursin("Hecke: matches", s)

    # ---- 4. "not assignable" is SAID, not concealed ------------------------
    # If the assignment does not hold, the API must warn instead of quietly giving an
    # index. The pairing is computed for every pair, so only true/false verdicts
    # remain.
    @test is_dl_basis(t) === true                       # checked and fine
    @test is_dl_basis(dl([2,1,2], [2,1,2])) === true     # resolved, no `missing`

    # and `show` states the verdict either way
    s2 = sprint(show, MIME"text/plain"(), dl([2,1,2], [2,1,2]))
    @test occursin("matches", s2)
    @test !occursin("no table entry", s2)

    # The FALSE case, produced artificially, so that the warning path is exercised.
    bad = DLTable(x, x, Dict(0 => t[0][1:1]))         # one degree, one leaf too few
    @test is_dl_basis(bad) === false
    s3 = sprint(show, MIME"text/plain"(), bad)
    @test occursin("MISMATCH", s3)
    @test occursin("NOT a basis", s3)
    # and it really does warn, not merely print
    @test_logs (:warn,) DiagrammaticHecke._dlb_warn_if_not_basis(bad)

    # the switch for sweeps returns the previous state and restores it
    old = dl_warnings!(false)
    @test old === true
    @test_logs DiagrammaticHecke._dlb_warn_if_not_basis(bad)         # silent: NO message
    @test dl_warnings!(old) === false
    @test dl_warnings!(true) === true                 # on again, the initial state
end

# ---------------------------------------------------------------------------
# Stage 0.2 — the circular counterpart and the decay matrix
# ---------------------------------------------------------------------------
#
# Over ALL 400 pairs from `a3_words()` of length <= 4 the decay matrix is square,
# block-diagonal and UNITRIANGULAR, so the triangularity conjecture holds there in the
# sharper form (blocks per degree). The assignment is MATCHED: at `2132 -> 2132` two
# rows of degree 2 are swapped.
#
# Only a few small cases stand here — the matrix costs a full reduction per double
# leaf.

@testset "DL basis API stage 0.2 (circular_dl, decay matrix, triangular form)" begin
    x = [1, 2, 1]

    # ---- 1. circular_dl is index-parallel to dl ---------------------------------
    ft = circular_dl(x)
    @test ft isa CircularDLTable
    @test circular_dl(x) === circular_dl(x, x)                    # from the cache
    @test dl_ranks(ft) == dl_ranks(dl(x))               # same graded ranks
    for d in dl_degrees(ft), j in eachindex(ft[d])
        @test ft[d][j].e == dl(x)[d][j].e               # same order
        @test ft[d][j].f == dl(x)[d][j].f
        @test ft[d][j].degree == dl(x)[d][j].degree
    end
    @test ft[0][1].morphism isa CircularMorphismGraph        # the circular image, unreduced
    @test length(ft[99]) == 0

    # `show` says expressly that the expensive check did NOT run
    clear_dl_cache!()
    s = sprint(show, MIME"text/plain"(), circular_dl(x))
    @test occursin("NOT checked", s)

    # ---- 2. the decay matrix -----------------------------------------------
    dec = circular_dl_matrix(x)
    @test dec isa CircularDLDecomposition
    @test length(dec.rows) == length(dl(x))                       # 12 DLs
    @test size(dec.M) == (length(dec.rows), length(dec.col_key))
    @test dec.row_degree == [en.degree for en in dl_entries(dl(x))]
    # every row really is a decay: at least one term
    @test all(any(!iszero, dec.M[i, :]) for i in 1:length(dec.rows))
    # the column degree is the smallest row degree in which the column occurs
    @test all(dec.col_degree[j] in dec.row_degree for j in 1:length(dec.col_key))

    # ---- 3. the triangularity conjecture -----------------------------------
    r = circular_dl_triangular(dec)
    @test r.square && r.triangular && r.block_diagonal && r.unitriangular
    @test sort(r.perm) == 1:length(dec.rows)                      # a permutation
    @test all(isone(dec.M[i, r.partner[i]]) for i in 1:length(dec.rows))
    @test isempty(r.excess_by_degree)
    @test is_circular_dl_basis(circular_dl(x)) === true
    @test is_circular_dl_basis(x, x) === true

    # the matched case: two rows of degree 2 are swapped
    @test is_circular_dl_basis([2, 1, 3, 2], [2, 1, 3, 2]) === true

    # ---- 4. the WARNING PATH with a reason ---------------------------------
    # Artificial: one column too many (this is what the known `distinct > dl` case
    # looks like). The warning must name the reason, not just say "does not work".
    n = length(dec.rows)
    Mbad = hcat(dec.M, fill(zero(SoergelPoly), n, 1))
    Mbad[1, end] = one(SoergelPoly)
    bad = CircularDLDecomposition(x, x, dec.rows, dec.row_degree,
                             vcat(dec.col_key, [:kuenstlich]),
                             vcat(dec.col_leaf, [nothing]),
                             vcat(dec.col_degree, [dec.row_degree[1]]), Mbad)
    rb = circular_dl_triangular(bad)
    @test rb.square === false
    @test rb.unitriangular === false
    @test rb.excess_by_degree == [dec.row_degree[1] => 1]         # +1 in this degree
end

# ---------------------------------------------------------------------------
# Stage 0.3 — the basis writer and the change of basis DL <-> circular
# ---------------------------------------------------------------------------
#
# `in_circular_dl_basis`: reduce_to_circular_leave gives circular-leaf coordinates `v`, and cT*M = v
# is solved against the decay matrix from stage 0.2.
# Anything unassignable lands in `rest`, `exact = false`, with @warn.

@testset "DL basis API stage 0.3 (in_circular_dl_basis / change of basis)" begin
    x = [1, 2, 1]
    t = dl(x)

    en = t[0][1]
    e2 = t[2][3]

    # ---- circular: basis element -> partner coefficient 1 -----------------------
    cf = in_circular_dl_basis(en.morphism)
    @test cf.basis === :circular_dl && cf.exact
    @test cf.coeffs == Dict((0, 1) => one(SoergelPoly))
    @test isempty(cf.rest)

    # and an element of higher degree (does the solution really use the matrix?)
    c4 = in_circular_dl_basis(e2.morphism)
    @test c4.exact
    @test c4.coeffs[(2, 3)] == one(SoergelPoly)

    # ---- 4. change of basis: M and M^-1 are exactly inverse ----------------
    M = dl_to_circular_matrix(x, x)
    Minv = circular_to_dl_matrix(x, x)
    @test Minv !== nothing
    n = size(M, 1)
    P = M * Minv                       # Matrix{SoergelPoly} product via Base
    @test all(P[i, j] == (i == j ? one(SoergelPoly) : zero(SoergelPoly))
              for i in 1:n, j in 1:n)

    # ---- 5. the solver itself (small, direct) ------------------------------
    a1 = alpha(1)
    Mk = [one(SoergelPoly) a1; zero(SoergelPoly) one(SoergelPoly)]
    v = [alpha(2), alpha(2) * a1 + one(SoergelPoly)]
    sol = DiagrammaticHecke._dlb_solve_left(Mk, v)
    @test sol == [alpha(2), one(SoergelPoly)]              # c1*1=α2; c1*α1+c2=...
    # non-constant pivot: a clean `nothing`, no throw
    Mbad = [a1 zero(SoergelPoly); zero(SoergelPoly) one(SoergelPoly)]
    @test DiagrammaticHecke._dlb_solve_left(Mbad, v) === nothing
end

# ---------------------------------------------------------------------------
# Stage 0.4 — the hflip decay (the cellular anti-involution)
# ---------------------------------------------------------------------------
#
# It computes with flip = bottom<->top, i.e. the HFLIP (axis = the horizon). Over all
# pairs from a3_words() of length <= 3, flip(DL_{(e,f)}: x->y) in the DL basis of
# (y,x) is the PURE SWAP (e,f) -> (f,e). This requires the cut1
# normalisation (_dlb_rot_leaves) in the writer — flip swaps the cuts, and without the
# rotation no canonical key matches.

@testset "DL basis API stage 0.4 (hflip decay = pure swap)" begin
    x = [1, 2, 1]

    # symmetric pair: permutation matrix, partner = the (f,e) swap
    r = hflip_dl_matrix(x, x)
    @test r.exact && r.pure_swap
    n = size(r.M, 1)
    @test sort(r.partner) == collect(1:n)              # a permutation
    rows = dl_entries(dl(x, x)); cols = dl_entries(dl(x, x))
    @test all(cols[r.partner[k]].e == rows[k].f &&
              cols[r.partner[k]].f == rows[k].e for k in 1:n)
    # involution: flipped twice is the identity matrix
    @test all((r.M * r.M)[i, j] == (i == j ? one(SoergelPoly) : zero(SoergelPoly))
              for i in 1:n, j in 1:n)

    # unequal pair: rows dl(x,y), columns dl(y,x)
    r2 = hflip_dl_matrix([1], [1, 2, 1])
    @test r2.exact && r2.pure_swap
    @test size(r2.M) == (length(dl([1], x)), length(dl(x, [1])))

    # the cut1 normalisation itself: the flip of a basis element can be written
    # exactly
    en = dl(x)[2][1]
    c = in_circular_dl_basis(flip(en.morphism))
    @test c.exact && length(c.coeffs) == 1

    # circular route (expensive, smallest pair): the same claim via in_circular_dl_basis
    r3 = hflip_dl_matrix([1], [1])
    @test r3.exact && r3.pure_swap
end

# ---------------------------------------------------------------------------
# The GENUINE vflip decay (left<->right)
# ---------------------------------------------------------------------------
#
# vflip(DL: x->y) is a map reverse(x) -> reverse(y) (plain name:
# hflip(::MorphismGraph)); it decays in the DL basis of (reverse(x), reverse(y)). There
# is NO bijection — no pure swap here.
#
# The decomposition runs through the circular writer. The two pairs 121/121 and
# 232/232 run through with `permutation = false`.

@testset "DL basis API (the genuine vflip decay)" begin
    # one-letter word: reverse is trivial, vflip = the identity of the basis
    r = vflip_dl_matrix([1], [1])
    @test r.exact
    @test size(r.M) == (length(dl_entries(dl([1], [1]))),
                        length(dl_entries(dl([1], [1]))))
    r2 = vflip_dl_matrix([1, 2], [1, 2])
    @test r2.exact
    @test size(r2.M) == (4, 4)
    r3 = vflip_dl_matrix([1, 2], [2, 1])
    @test r3.exact

    # --- exact in the circular path AND not a permutation --------------------
    r121 = vflip_dl_matrix([1, 2, 1], [1, 2, 1])
    @test r121.exact
    @test !r121.permutation
    @test size(r121.M) == (12, 12)
    r232 = vflip_dl_matrix([2, 3, 2], [2, 3, 2])
    @test r232.exact
    @test !r232.permutation

    # the default IS circular: setting it explicitly gives the same
    @test vflip_dl_matrix([1, 2, 1], [1, 2, 1]).M == r121.M
end

# ---------------------------------------------------------------------------
# vflip is NOT a bijection of the DL basis — counterexample x=[1,1]
# ---------------------------------------------------------------------------
#
# The smallest case where this becomes visible: x=[1,1] (dot + identity). Unlike
# [1,2,1]/[2,3,2] this case is cheap (8x8, no apply_d4 boundary) and exact.
@testset "vflip is not a bijection — x=[1,1]" begin
    r11 = vflip_dl_matrix([1, 1], [1, 1])
    @test r11.exact
    @test !r11.permutation
    @test size(r11.M) == (8, 8)
    # at least one row with more than one non-zero entry
    @test any(count(l -> !iszero(r11.M[k, l]), 1:8) > 1 for k in 1:8)
end

# ---------------------------------------------------------------------------
# Stage 0.5 — structure constants of the composition
# ---------------------------------------------------------------------------
#
# At x=y=z=121 (144 compositions) 142 land in the table, 2 throw (the apply_d4
# boundary) and 6 are not exact. The tests stay with the small, exact cases.

@testset "DL basis API stage 0.5 (compose_in_basis / structure_constants)" begin
    t = dl([1])
    # id∘id = id, dot∘dot = α₁·dot (barbell), id then dots = dots
    @test compose_in_basis(t[0][1], t[0][1]).coeffs == Dict((0, 1) => one(SoergelPoly))
    @test compose_in_basis(t[2][1], t[2][1]).coeffs == Dict((2, 1) => alpha(1))
    @test compose_in_basis(t[0][1], t[2][1]).coeffs == Dict((2, 1) => one(SoergelPoly))
    # non-composable throws with a clear message
    @test_throws ErrorException compose_in_basis(dl([1])[0][1], dl([2])[0][1])

    r = structure_constants([1], [1], [1])
    @test r.exact && length(r.table) == 4 && isempty(r.thrown)
    @test r.table[((0, 1), (0, 1))].coeffs == Dict((0, 1) => one(SoergelPoly))

    r2 = structure_constants([1, 2], [1, 2], [1, 2])
    @test r2.exact && length(r2.table) == 16 && isempty(r2.thrown)

    # the circular route confirms the barbell case
    cf = compose_in_basis(t[2][1], t[2][1])
    @test cf.exact && cf.coeffs == Dict((2, 1) => alpha(1))
end

