# algebra/Pairing.jl
#
# The BOTT–SAMELSON SOERGEL PAIRING, computed in the Hecke algebra over the
# Coxeter group.
#
# WHAT IT IS. For words x = i1…in, y = j1…jm over the generators,
#
#     pairing(x, y) = ⟨bs(x), bs(y)⟩ = grdim Hom(BS(x), BS(y)),
#
# where `bs(w) = b_{w1} · … · b_{wk}` in the Hecke algebra H of the Coxeter
# group. This is the ground truth that `double_leaves(x, y)` (our purely
# diagrammatic construction) is checked against — see test/pairing.jl.
#
# THE ONE IDENTITY THAT MAKES IT CHEAP. The Soergel pairing is
# `⟨h, h'⟩ = ε(h · ω(h'))` with `ε` the coefficient of `δ_e` and `ω` the
# antiautomorphism fixing each `δ_s`. Every `b_s` is ω-invariant, so
# `ω(bs(y)) = bs(reverse(y))`, and a product of `bs`s is the `bs` of the
# concatenated word. Hence
#
#     pairing(x, y) = [δ_e] bs(x ++ reverse(y))
#
# — ONE left-to-right pass multiplying by generators, then read off one
# coefficient. No general Hecke multiplication is needed.
#
# THE ALGEBRA. H is free over Z[v, v^-1] on the standard basis {δ_w}, with
# `δ_s² = (v^-1 - v) δ_s + δ_e` (EMTW convention) and `b_s = δ_s + v·δ_e`.
# Multiplying a basis element by `b_s` therefore needs only a length comparison:
#
#     δ_w · b_s = δ_ws + v·δ_w        if ℓ(ws) > ℓ(w)     (s not a right descent)
#               = δ_ws + v^-1·δ_w     if ℓ(ws) < ℓ(w)     (s a right descent)
#
# (down case: δ_w·δ_s = (v^-1 - v)δ_w + δ_ws, and the +v·δ_w from b_s's second
# term cancels the -v·δ_w, leaving v^-1.)
#
# GRADING. Degrees are v-exponents and MAY BE NEGATIVE — e.g.
# `pairing(121321, 121321)` has a term in degree -2. A negative degree is a real
# leaf with more trivalents than dots. Beware Julia's `%` on negatives: use
# `mod(d, 2)`.

# ---- the group, built once --------------------------------------------------

const _PAIRING_GROUP = Ref{Any}(nothing)
function _pairing_group()
    _PAIRING_GROUP[] === nothing && (_PAIRING_GROUP[] = coxeter_group_a3())
    return _PAIRING_GROUP[][1]                       # (W, gens) -> W
end

# ---- Hecke elements: Dict(w => Dict(v-degree => coefficient)) ----------------

_hecke_addcoeff!(poly::Dict{Int,Int}, d::Int, c::Int) = begin
    c == 0 && return poly
    n = get(poly, d, 0) + c
    n == 0 ? delete!(poly, d) : (poly[d] = n)
    return poly
end

"Multiply the Hecke element `h` on the right by `b_s`."
function _hecke_times_b(h::Dict{K,Dict{Int,Int}}, s::Int) where {K}
    out = Dict{K,Dict{Int,Int}}()
    for (w, poly) in h
        ws = right_multiply(w, s)
        shift = coxeter_length(ws) < coxeter_length(w) ? -1 : 1
        for (d, c) in poly
            _hecke_addcoeff!(get!(out, ws, Dict{Int,Int}()), d, c)          # δ_ws
            _hecke_addcoeff!(get!(out, w,  Dict{Int,Int}()), d + shift, c)  # v^±1 δ_w
        end
    end
    for (w, poly) in collect(out)
        isempty(poly) && delete!(out, w)
    end
    return out
end

"`bs(word) = b_{word[1]} · … · b_{word[end]}` as a Hecke element."
function _hecke_bs(word::AbstractVector{<:Integer})
    W = _pairing_group()
    e = one(W)
    h = Dict{typeof(e),Dict{Int,Int}}(e => Dict(0 => 1))
    for s in word
        h = _hecke_times_b(h, Int(s))
    end
    return h
end

const _PAIRING_CACHE = Dict{Tuple{Vector{Int},Vector{Int}},Dict{Int,Int}}()

"""
    bs_pairing(x::Vector{Int}, y::Vector{Int}) -> Dict{Int,Int}

`grdim Hom(BS(x), BS(y))` as `Dict(v-degree => coefficient)`, zero coefficients
omitted. Computed as `[δ_e] bs(x ++ reverse(y))` in the Hecke algebra (see the
file header) and memoised.

Unlike the table this replaced, it is defined for **any** pair of words, not
only for the reduced ones that happened to be tabulated.
"""
function bs_pairing(x::Vector{Int}, y::Vector{Int})
    get!(_PAIRING_CACHE, (x, y)) do
        W = _pairing_group()
        get(_hecke_bs(vcat(x, reverse(y))), one(W), Dict{Int,Int}())
    end
end

"""
    pairing_table() -> Dict{Tuple{Vector{Int},Vector{Int}},Dict{Int,Int}}

The pairing for every pair of short-lex reduced words of the group (576 pairs in
A₃), in the shape the tests iterate over. Computed on first call and cached.

For a single pair prefer [`bs_pairing`](@ref) — it does not build the table and
accepts non-reduced words too.
"""
const _PAIRING_TABLE = Ref{Union{Nothing,Dict{Tuple{Vector{Int},Vector{Int}},Dict{Int,Int}}}}(nothing)
function pairing_table()
    _PAIRING_TABLE[] !== nothing && return _PAIRING_TABLE[]
    W = _pairing_group()
    words = [Int[Int(c) for c in short_lex(w)] for w in enumerate_whole_group(W)]
    table = Dict{Tuple{Vector{Int},Vector{Int}},Dict{Int,Int}}()
    for x in words, y in words
        table[(x, y)] = bs_pairing(x, y)
    end
    _PAIRING_TABLE[] = table
    return table
end

"""
    expected_count(x::Vector{Int}, y::Vector{Int}; degree = nothing) -> Int

The expected number of double leaves, i.e. `grdim Hom(BS(x), BS(y))` evaluated at
`v = 1`, or its degree-`d` slice.

- `degree = nothing` (default): total over all degrees;
- `degree = d`: just the coefficient at degree `d` (`0` if absent).
"""
function expected_count(x::Vector{Int}, y::Vector{Int}; degree::Union{Nothing,Int} = nothing)
    degs = bs_pairing(x, y)
    degree === nothing && return sum(values(degs); init = 0)
    return get(degs, degree, 0)
end
