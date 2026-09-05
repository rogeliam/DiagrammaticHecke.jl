# morphism/CircularPairing.jl — the pairing ⟨DL_i, DL_j⟩_k as a function in the
# package.
#
# WHAT THIS IS ABOUT. Two Double Leaves over the same word `x` (with `z = ε`)
# are endomorphisms of `x`, hence composable. `pairing(dl_i, dl_j)` composes
# them and decomposes the result into CIRCULAR LEAVES (not back into the DL basis
# like `compose_in_basis`/`in_circular_dl_basis` in DLBasisWriter.jl): as long
# as the circular basis can have more circular leaves than
# DLs (`circular_dl_triangular(...).square == false`), the assignment
# circular-leaf → DL is not always unique; the pairing should therefore live
# directly in circular-leaf coordinates, with a hint WHEN that ambiguity occurs in
# the affected degree — not a blocker, the linear combination is always
# returned.
#
# DEGREE BOUND. `reduce_to_circular_leave`
# (circular/rules/CircularDecoratedRules.jl) is a deeply nested recursion over several
# term-generating steps (fusion, P1, 2parallel, Zamo, D4), each
# with its own termination argument and its own `node_growth_limit`. Some
# steps (fusion, D4) change the degree temporarily before it drops again, so
# a degree filter placed inside that recursion could cut off an intermediate
# term silently, before its degree drops back down. `reduce_to_circular_leave`
# therefore runs untouched to fixpoint, and `pairing` calls it normally. The
# degree bound is applied AFTER reduction, as a diagnostic (`degree` on
# `CircularPairing`) and for the visibility of the mismatch hint.

"""
    CircularPairing

The decomposition of the composition `dl_i ∘ dl_j` into CIRCULAR LEAVES (not into
the DL basis — see [`compose_in_basis`](@ref) for that):

* `x` — the common word (`dl_i`/`dl_j` are endomorphisms of `x`);
* `coeffs` — `circular_canonical_key(::CircularDecorated) => coefficient ∈ R`;
* `leaves` — same key => a representative `CircularDecorated`, for drawing;
* `degree` — `dl_i.degree + dl_j.degree`, the degree of the composition.
"""
struct CircularPairing
    x::Vector{Int}
    coeffs::Dict{Any,SoergelPoly}
    leaves::Dict{Any,Any}
    degree::Int
end

"""
    pairing(dl_i::DLEntry, dl_j::DLEntry) -> CircularPairing

The structure constant ⟨DL_i, DL_j⟩_k: composes two Double Leaves
over the same word (`compose(dl_i.morphism, dl_j.morphism)` — `dl_i` is
applied first, as everywhere in the package), reduces with
[`reduce_to_circular_leave`](@ref) and reads off the circular-leaf coefficients
directly — WITHOUT decomposing back into the DL basis, since that assignment
is exactly what's not unique under a circular/DL mismatch.

Throws if `top(dl_i.morphism) != bottom(dl_j.morphism)` (not composable — in
particular `dl_i`/`dl_j` must be endomorphisms of the same word, `x = z`).
"""
function pairing(dl_i::DLEntry, dl_j::DLEntry)
    top(dl_i.morphism) == bottom(dl_j.morphism) || error(
        "pairing: top(dl_i) = $(_dlb_tok(top(dl_i.morphism))) doesn't match " *
        "bottom(dl_j) = $(_dlb_tok(bottom(dl_j.morphism))) — dl_i/dl_j must " *
        "be endomorphisms of the same word.")
    h = compose(dl_i.morphism, dl_j.morphism)
    fm = circular_morphism(h)
    fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
    combo = reduce_to_circular_leave(fdm)

    coeffs = Dict{Any,SoergelPoly}()
    leaves = Dict{Any,Any}()
    for (d0, c0) in pairs_of(combo)
        # THE OUTER LABEL BELONGS IN THE COEFFICIENT: a polynomial with no
        # cell to sit in stays in `outer_label` instead of becoming a factor
        # — the barbell from `dot∘dot` carries its `α₁` there.
        # `circular_canonical_key` sees the label, so without pulling it in,
        # such terms would count as SEPARATE leaves with coefficient 1, and
        # the value of the pairing coefficient would be lost (`pairing` on
        # `121212` would give `1` instead of `α₁²α₂ + α₁α₂²`).
        # Pulling it in here is exactly right: the outer region is not a
        # cell, its polynomial is a SCALAR factor of the term — and the
        # pairing value is a scalar. The rule governs the DIAGRAM, not
        # the coordinates (literally the same rationale as in `in_dl_basis`,
        # morphism/DLBasisWriter.jl).
        d, c = isone(d0.outer_label) ? (d0, c0) :
               (CircularDecorated(d0.graph, d0.region_labels), c0 * d0.outer_label)
        k = circular_canonical_key(d)
        coeffs[k] = get(coeffs, k, zero(SoergelPoly)) + c
        leaves[k] = d
    end
    filter!(p -> !iszero(p.second), coeffs)

    return CircularPairing(bottom(dl_i.morphism), coeffs, leaves,
                       dl_i.degree + dl_j.degree)
end

# ---- the mismatch hint -----------------------------------------------------
#
# As long as #circular-leaves != #DLs (in the affected degree), a HINT
# appears — the decomposition is printed regardless. Once the circular/DL
# mismatch is fixed, `circular_dl_triangular` returns an empty
# `excess_by_degree` and the hint disappears automatically.
#
# Costs a `circular_dl_matrix` computation (expensive, see DLBasis.jl); it runs
# cached via `_FDL_MATRIX_CACHE`, so it's only actually paid for on the first
# `show` of a word.

function _pairing_excess_hint(p::CircularPairing)
    dec = try
        circular_dl_matrix(p.x, p.x)
    catch
        return nothing        # reduction throws — no statement possible, no hint forced
    end
    r = circular_dl_triangular(dec)
    isempty(r.excess_by_degree) && return nothing
    hits = [d => n for (d, n) in r.excess_by_degree if d <= p.degree]
    isempty(hits) && return nothing
    return hits
end

function Base.show(io::IO, ::MIME"text/plain", p::CircularPairing)
    println(io, "CircularPairing $(_dlb_tok(p.x)) (degree ", p.degree, ")")
    ks = sort!(collect(keys(p.coeffs)); by = k -> string(k))
    if isempty(ks)
        println(io, "  0 (empty sum)")
    else
        for k in ks
            println(io, "  ", p.coeffs[k], "  ·  ", p.leaves[k])
        end
    end
    hint = _pairing_excess_hint(p)
    hint === nothing || println(io,
        "  ⚠ circular basis still has more circular leaves than DLs in degree(s) ",
        join(("$d (+$n)" for (d, n) in hint), ", "),
        " — the assignment circular-leaf → DL is not yet unique there. ",
        "Once fixed, this hint disappears on its own.")
end

Base.show(io::IO, p::CircularPairing) =
    print(io, "CircularPairing(", _dlb_tok(p.x), ", degree ", p.degree, ", ",
          length(p.coeffs), " circular leaf(s))")
