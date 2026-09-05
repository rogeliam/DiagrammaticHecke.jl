# circular/rules/CircularRexRelation.jl — the rex-cycle relation
# `id_{w·i} = [route B] − [remainder terms ⊗ id_i]` as a DERIVED, CACHED function.
#
# DERIVED, NOT COPIED (the `_zamo_build_z2` principle, ZamoTermRules.jl): both
# routes are recomputed from the rule set on first access, and the following
# situation is asserted — route A
# (reduce the cycle `w → w' → w` tensored with the i-strand directly) contains
# the identity summand on `w·i`, route B (FORCE the fusion at an i-pair first)
# does not, and A ≠ B. If the rule set drifts, the assertions fail loudly
# instead of a frozen fixture drifting silently.
#
# NOT every i-pair carries the relation: on `212321` the pair (1,9)
# is INERT — its route B still contains the identity and equals A. The pairs
# are tried in order and the FIRST one with the asserted situation wins.
#
# THE RELATION. With `ID = cl_reduce(id_{w·i})`, `A` and `B` as above (all three
# normalised the same way, outer labels pulled into the coefficients):
#
#     rhs := B − A + ID     and the claim is   ID == rhs   (as morphisms),
#
# where `rhs` contains NO term of `ID` (asserted): the identity terms of `A`
# cancel exactly against `ID`. So `rhs` expresses the identity on `w·i` in
# other diagrams — the right-hand side splices in.
#
# DISK CACHE: a derivation costs 10 s (one braid move) to
# 84 s (three moves) plus the inert pairs tried before the right one — cached
# as versioned text under `data/cache/rex/`, template `.cd`/`.cdb`
# (CircularDecoratedIO.jl). The cache stores BOTH sides; the assertions run at
# DERIVATION time, not at load time (same policy as the `.cdb` basis files).

"Default directory for the cached rex relations."
circular_rex_cache_dir() = joinpath(dirname(dirname(dirname(@__DIR__))), "data", "cache", "rex")

"File name of the cached relation for `(w, i)`."
circular_rex_cache_path(w::Vector{Int}, i::Int;
                        dir::AbstractString = circular_rex_cache_dir()) =
    joinpath(dir, string("rex-", _fdb_word_token(w), "-i", i, ".rexrel"))

const _REXREL_VERSION = 1

"""
    CIRCULAR_REXREL_DERIVE_WITHOUT_REX

`true` (default): while a rex relation is DERIVED (cache miss), the rex fusion is
switched off for that derivation — see the note in [`circular_rex_relation`](@ref)
(Q33b: the derivation would otherwise ask for the relation it is deriving).
`false` leaves the rex-fusion criterion on during derivation too, for measurements only.
"""
const CIRCULAR_REXREL_DERIVE_WITHOUT_REX = Ref(true)

# ---- (de)serialisation -------------------------------------------------------

function _rexrel_combo_to_io(io::IO, c::CircularComboR)
    for (d0, coeff) in sort(pairs_of(c); by = p -> string(circular_canonical_key(p[1])))
        d = d0 isa CircularDecorated ? d0 :
            CircularDecorated(d0, fill(one(SoergelPoly), region_count(d0)))
        println(io, "coeff ", soergelpoly_to_string(coeff))
        print(io, circulardecorated_to_string(d))
        println(io, "---")
    end
end

function _rexrel_combo_from_lines(lines::Vector{<:AbstractString})
    out = CircularComboR()
    coeff = nothing
    body = String[]
    for line in lines
        if startswith(strip(line), "coeff ")
            coeff = soergelpoly_from_string(strip(line)[7:end])
            empty!(body)
        elseif strip(line) == "---"
            coeff === nothing && error(".rexrel: term without a coeff line")
            d = circulardecorated_from_string(join(body, '\n'))
            out = out + coeff * CircularComboR(d)
            coeff = nothing
        else
            push!(body, line)
        end
    end
    coeff === nothing || error(".rexrel: trailing term without closing ---")
    return out
end

function _rexrel_save(path::AbstractString, w::Vector{Int}, i::Int,
                      pair::Tuple{Int,Int}, id::CircularComboR, rhs::CircularComboR)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "#rexrel ", _REXREL_VERSION)
        println(io, "w ", join(w, ' '))
        println(io, "i ", i)
        println(io, "pair ", pair[1], ' ', pair[2])
        println(io, "[id]")
        _rexrel_combo_to_io(io, id)
        println(io, "[rhs]")
        _rexrel_combo_to_io(io, rhs)
    end
    return path
end

function _rexrel_load(path::AbstractString, w::Vector{Int}, i::Int)
    isfile(path) || return nothing
    lines = split(read(path, String), '\n')
    length(lines) >= 5 || return nothing
    v = parse(Int, split(strip(lines[1]))[2])
    v == _REXREL_VERSION || error("unsupported .rexrel version $v (expected $_REXREL_VERSION)")
    strip(lines[2]) == "w " * join(w, ' ') && strip(lines[3]) == "i $i" || error(
        ".rexrel header does not match (w, i) — wrong file at $path")
    pr = parse.(Int, split(strip(lines[4]))[2:3])
    iid = findfirst(==("[id]"), strip.(lines))
    irh = findfirst(==("[rhs]"), strip.(lines))
    (iid === nothing || irh === nothing) && error(".rexrel: missing [id]/[rhs] section")
    id  = _rexrel_combo_from_lines(collect(lines[iid+1:irh-1]))
    rhs = _rexrel_combo_from_lines(collect(lines[irh+1:end]))
    return (id = id, rhs = rhs, pair = (pr[1], pr[2]))
end

# ---- the derivation ----------------------------------------------------------

"The pure 2parallel surgery at a GIVEN pair, as the two-term sum with α_i/2 —
the route-B move."
function _rexrel_forced_2parallel(fdm::CircularDecoratedMorphism, ei::Int, ej::Int)
    g = fdm.m.graph
    f = circular_2parallel_apply(g, ei, ej)
    f === nothing && return nothing
    half = (1//2) * alpha(g.edges[ei].colour)
    out = CircularComboR()
    for r in f.newregs
        labels = _circular_transfer_labels(g, f.graph, fdm.region_labels;
                                           overrides = Dict(r => half))
        out = out + CircularComboR(CircularDecorated(f.graph, labels))
    end
    return out
end

"All pairs of parallel `col`-coloured edges (common region, not touching)."
function _rexrel_colour_pairs(g::CircularGraph, col::Int)
    adj, _ = circular_region_adjacency(g)
    out = Tuple{Int,Int}[]
    for R in 1:region_count(g)
        lst = adj[R]
        for a in 1:length(lst), b in (a + 1):length(lst)
            ei = lst[a][2]; ej = lst[b][2]
            (g.edges[ei].colour == col && g.edges[ej].colour == col) || continue
            _circular_edges_touch(g.edges[ei], g.edges[ej]) && continue
            push!(out, (min(ei, ej), max(ei, ej)))
        end
    end
    return unique!(out)
end

"Reduce every term of a (labelled) combination fully, normalised as cl_reduce."
function _rexrel_reduce_combo(c::CircularComboR, cut1::Int, cut2::Int)
    acc = CircularComboR()
    for (dd, coeff) in pairs_of(c)
        m3 = CircularMorphismGraph(dd.graph, cut1, cut2)
        acc = acc + coeff * reduce_to_circular_leave(
                  CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label))
    end
    return _cl_normalize(acc)
end

_rexrel_keys(c::CircularComboR) = Set(circular_canonical_key(d) for (d, _) in pairs_of(c))

function _rexrel_derive(w::Vector{Int}, i::Int; maxpairs::Int = 8)
    is_reduced(w) || throw(ArgumentError("circular_rex_relation: w = $w is not reduced"))
    is_reduced(vcat(w, i)) && throw(ArgumentError(
        "circular_rex_relation: w·i is reduced — i = $i is no right descent of $w"))

    wp, _ = braid_to_end_with(w, i)
    f = braid_top_to(light_leaf_up(w), wp)
    cyc_i = tensor(compose(f, flip(f)), identity_strand(i))
    fm = circular_morphism(cyc_i)
    g = fm.graph
    fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(g)))

    ID = cl_reduce(light_leaf_up(vcat(w, i)))
    idk = _rexrel_keys(ID)

    # route A: reduce directly — must contain the identity summand.
    A = _cl_normalize(reduce_to_circular_leave(fdm))
    kA = _rexrel_keys(A)
    issubset(idk, kA) || error(
        "circular_rex_relation($(join(w)), $i): route A does not contain the " *
        "identity summand; the relation cannot be derived with the current rule set")

    # route B: force the fusion, first pair with the asserted situation wins.
    tried = 0
    for (ei, ej) in _rexrel_colour_pairs(g, i)
        tried >= maxpairs && break
        c = _rexrel_forced_2parallel(fdm, ei, ej)
        c === nothing && continue
        tried += 1
        B = _rexrel_reduce_combo(c, fm.cut1, fm.cut2)
        (isempty(intersect(_rexrel_keys(B), idk)) && A != B) || continue   # inert pair
        rhs = B - A + ID
        isempty(intersect(_rexrel_keys(rhs), idk)) || error(
            "circular_rex_relation($(join(w)), $i): the identity terms of route A " *
            "do not cancel against ID at pair ($ei,$ej)")
        return (id = ID, rhs = rhs, pair = (ei, ej))
    end
    error("circular_rex_relation($(join(w)), $i): no i-pair with route B free of " *
          "the identity and A ≠ B among the first $tried constructible pairs; " *
          "the relation is not derivable here")
end

"""
    circular_rex_relation(w::Vector{Int}, i::Int; dir = circular_rex_cache_dir(),
                          recompute = false, maxpairs = 8)
        -> (id::CircularComboR, rhs::CircularComboR, pair::Tuple{Int,Int})

The rex-cycle relation `id == rhs` for the pair `(w, i)` —
`id` is the reduced identity on `w·i`, `rhs` expresses it in OTHER diagrams
(file head: `rhs = route B − remainder terms`, no `id` term left, asserted).
`pair` is the i-pair of route B that carried the fusion.

`w` must be reduced with right descent `i` (checked). Derived on first call
(10–90 s, file head) and cached on disk under `dir`; `recompute = true`
ignores and overwrites the cache.
"""
function circular_rex_relation(w::Vector{Int}, i::Int;
                               dir::AbstractString = circular_rex_cache_dir(),
                               recompute::Bool = false, maxpairs::Int = 8)
    path = circular_rex_cache_path(w, i; dir = dir)
    if !recompute
        got = _rexrel_load(path, w, i)
        got === nothing || return got
    end
    # ⚠ DERIVE WITH THE CRITERION OFF. The derivation reduces `id_{w·i}` and
    # the two routes with the FULL driver. With `CIRCULAR_REX_FUSION_ENABLED`
    # on, the driver's region-word criterion fires inside that very
    # reduction — the identity diagram of the UNREDUCED word `w·i` is
    # exactly the figure it fires on — and asks for the relation `(w, i)`
    # that is being derived right now: cache miss, `_rexrel_derive` again,
    # StackOverflow (observed at `id_1213213`/`id_1213212`, every relation
    # that is not yet on disk). Deriving a relation and applying it are two
    # different jobs; the cut is at this boundary. Switching
    # `CIRCULAR_REXREL_DERIVE_WITHOUT_REX[]` off reaches that reduction path
    # too, for measurement.
    rex_was = CIRCULAR_REX_FUSION_ENABLED[]
    CIRCULAR_REXREL_DERIVE_WITHOUT_REX[] && (CIRCULAR_REX_FUSION_ENABLED[] = false)
    r = try
        _rexrel_derive(w, i; maxpairs = maxpairs)
    finally
        CIRCULAR_REX_FUSION_ENABLED[] = rex_was
    end
    _rexrel_save(path, w, i, r.pair, r.id, r.rhs)
    return r
end
