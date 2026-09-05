# morphism/CLBasis.jl — the Circular-L basis of a word pair, COMPUTED and CACHED.
#
# WHAT'S HERE. Three things every script needs:
#
#  1. `cl_reduce(m)` — reduces a morphism `x → y` to its circular-leaf sum:
#     `circular_morphism` → `CircularDecoratedMorphism` → (cut1 rotated to 0, as in
#     `in_circular_dl_basis`) → `reduce_to_circular_leave`, with the zamo switch as an
#     argument and the OUTER LABEL PULLED INTO THE COEFFICIENT (same rule and
#     rationale as in `in_circular_dl_basis`/`pairing`).
#  2. `cl_basis(x, y)` — every double leaf of `(x, y)` reduced once this way,
#     yielding the list of DISTINCT circular leaves per the counting rule of
#     `circular_leaf_degrees` (a leaf counts if it occurs at least once with a
#     SCALAR coefficient). Result is cached as `.cdb` (circular/CircularDecoratedIO.jl) —
#     the expensive part is the reductions, and the same basis is needed
#     repeatedly across many computations.
#  3. `cl_coordinates(combo, b)` — read a reduced sum off in this basis:
#     coefficient vector plus the terms whose leaf is NOT in the basis
#     (`rest`). A nonempty `rest` is a FINDING, not an error — it means the
#     pair's basis is incomplete, or the reduction leaves it.
#
# DEGREE BOUND. `cl_basis(x, y; degree_max = 0)` crops the RESULT to leaves
# with `circular_degree <= degree_max` after the computation (a degree-0
# morphism can only carry leaves of degree ≤ 0, since `deg(coeff) =
# −circular_degree ≥ 0`). The FULL basis is always what gets cached — cropping is
# free, recomputing takes hours.

"""
Default directory for the `.cdb` basis files.

A cache is generated output, so it is NOT written inside the package tree:
`DH_CACHE_DIR` if that environment variable is set, otherwise a
`diagrammatichecke-cl` folder under `tempdir()`. Every entry point takes a
`dir =` keyword, so a caller can put it anywhere.
"""
cl_cache_dir() = get(ENV, "DH_CACHE_DIR",
                     joinpath(tempdir(), "diagrammatichecke-cl"))

"""
    cl_cache_path(x, y; zamo = false, dl_degree_max = nothing, dir = cl_cache_dir())

File name of a basis: word pair + switch setting (zamo-free is a DIFFERENT
basis) + the DL degree bound it was computed under (`-dl<g>`). The bound
belongs in the NAME, not just the header: a basis computed from the degree-≤-0
DLs is not a cropped version of the full one but a smaller computation — and
it must not overwrite the full one.
"""
cl_cache_path(x::Vector{Int}, y::Vector{Int}; zamo::Bool = false,
                dl_degree_max::Union{Nothing,Int} = nothing,
                dir::AbstractString = cl_cache_dir()) =
    joinpath(dir, string(_fdb_word_token(x), "-", _fdb_word_token(y),
                         zamo ? "" : "-zamofrei",
                         dl_degree_max === nothing ? "" : "-dl$(dl_degree_max)",
                         ".cdb"))

# The outer label belongs in the COEFFICIENT, not the key (see the rationale
# box in `pairing`, morphism/CircularPairing.jl).
function _cl_normalize(c::CircularComboR)
    out = CircularComboR()
    for (d0, coeff0) in pairs_of(c)
        d, coeff = isone(d0.outer_label) ? (d0, coeff0) :
                   (CircularDecorated(d0.graph, d0.region_labels), coeff0 * d0.outer_label)
        out = out + coeff * CircularComboR(d)
    end
    return out
end

"""
    cl_reduce(m::MorphismGraph; zamo = false) -> CircularComboR
    cl_reduce(fdm::CircularDecoratedMorphism; zamo = false) -> CircularComboR

The circular-leaf sum of `m`. `zamo` sets [`CIRCULAR_ZAMO_ENABLED`](@ref) for the
duration of the reduction (default **off** — measures zamo-free) and
restores the previous setting afterward. Leaf rotation (`cut1 != 0`) is
normalized to 0 as in `in_circular_dl_basis`, so the keys are comparable to those
of the double leaves.
"""
function cl_reduce(fdm::CircularDecoratedMorphism; zamo::Bool = false)
    r = _dlb_cut_rot(fdm.m.cut1, length(letters(fdm.m.graph.word)))
    if r != 0
        all(isone, fdm.region_labels) || error(
            "cl_reduce: leaf rotation (cut1 = $(fdm.m.cut1)) with nontrivial " *
            "region labels is not supported")
        g2 = _dlb_rot_leaves(fdm.m.graph, r)
        fm2 = CircularMorphismGraph(g2, 0, length(bottom(fdm.m)))
        fdm = CircularDecoratedMorphism(fm2, fill(one(SoergelPoly), region_count(g2)))
    end
    old = CIRCULAR_ZAMO_ENABLED[]
    CIRCULAR_ZAMO_ENABLED[] = zamo
    try
        return _cl_normalize(reduce_to_circular_leave(fdm))
    finally
        CIRCULAR_ZAMO_ENABLED[] = old
    end
end

function cl_reduce(m::MorphismGraph; zamo::Bool = false)
    fm = circular_morphism(m)
    return cl_reduce(
        CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)));
        zamo = zamo)
end

"""
    cl_basis(x, y; zamo = false, dl_degree_max = nothing, degree_max = nothing,
               recompute = false, dir = cl_cache_dir(), progress = false) -> CLBasis

The Circular-L basis of the word pair `(x, y)`: every double leaf of
`double_leaves(x, y)` run once through [`cl_reduce`](@ref), then the
DISTINCT circular leaves per the counting rule of [`circular_leaf_degrees`](@ref) (a
leaf counts if it occurs at least once with a scalar coefficient), in order of
first occurrence.

**`dl_degree_max` is the cheap route**: restricting to degree ≤ 0 is enough,
obtained by reducing only the DLs of degree ≤ 0, i.e. far fewer. Only double
leaves with `degree <= dl_degree_max` are reduced at
all. For `(121321, 321323)` that's **2 instead of 204** — seconds instead of
hours. The cost is named honestly: whatever a higher-degree double leaf would
contribute to lower-degree leaves is invisible to this basis. `dl_degree_max =
nothing` (default) computes all of them.

Cached as `.cdb` under `dir` — the loader only accepts the file if its SWITCH
SETTING matches what was requested (`load_cl_basis`); the DL bound is
encoded in the file name (`-dl0`). `recompute = true` forces recomputation.
`degree_max` crops the RESULT. `progress = true` writes a progress line per
double leaf to `stderr` — for long runs (rule: visible progress).

**Note on the name**: whether this set actually IS a basis is not established;
it is more precisely a `leaf_set`. The function is named `cl_basis` because
that name is public and appears in the header of the `.cdb` files.
"""
function cl_basis(x::Vector{Int}, y::Vector{Int}; zamo::Bool = false,
                    dl_degree_max::Union{Nothing,Int} = nothing,
                    degree_max::Union{Nothing,Int} = nothing,
                    recompute::Bool = false,
                    dir::AbstractString = cl_cache_dir(),
                    progress::Bool = false)
    path = cl_cache_path(x, y; zamo = zamo, dl_degree_max = dl_degree_max, dir = dir)
    b = recompute ? nothing : load_cl_basis(path; zamo = zamo)
    if b === nothing
        b = _cl_basis_compute(x, y; zamo = zamo, dl_degree_max = dl_degree_max,
                                progress = progress)
        save_cl_basis(path, b)
    end
    return degree_max === nothing ? b : cl_basis_filter(b; degree_max = degree_max)
end

"""
    cl_basis_filter(b::CLBasis; degree_max) -> CLBasis

`b` cropped to the leaves with `circular_degree <= degree_max`. `ndl`/`zamo`/
`pkgversion` are kept as-is (provenance of the full run).
"""
function cl_basis_filter(b::CLBasis; degree_max::Int)
    keep = [i for i in eachindex(b.degrees) if b.degrees[i] <= degree_max]
    return CLBasis(b.x, b.y, b.leaves[keep], b.degrees[keep], b.zamo,
                     b.pkgversion, b.ndl)
end

function _cl_basis_compute(x::Vector{Int}, y::Vector{Int}; zamo::Bool,
                             dl_degree_max::Union{Nothing,Int} = nothing,
                             progress::Bool)
    dls = double_leaves(x, y)
    dl_degree_max === nothing || filter!(dl -> dl.degree <= dl_degree_max, dls)
    order  = Any[]                      # keys in order of first occurrence
    reps   = Dict{Any,CircularDecorated}()
    gdeg   = Dict{Any,Int}()
    scalar = Dict{Any,Bool}()
    for (i, dl) in enumerate(dls)
        progress && print(stderr, "\r  cl_basis($(_dlb_tok(x))→$(_dlb_tok(y))): ",
                          "DL $i/$(length(dls)) (degree $(dl.degree)), ",
                          length(order), " leaves    ")
        combo = cl_reduce(dl.morphism; zamo = zamo)
        for (d, c) in pairs_of(combo)
            k = circular_canonical_key(d)
            if !haskey(gdeg, k)
                push!(order, k); reps[k] = d; gdeg[k] = circular_degree(d); scalar[k] = false
            end
            scalar[k] |= iszero(degree(c))
        end
    end
    progress && println(stderr)
    keep = [k for k in order if scalar[k]]
    return CLBasis(x, y, CircularDecorated[reps[k] for k in keep], Int[gdeg[k] for k in keep],
                     zamo, _cl_pkgversion(), length(dls))
end

# Package version from Project.toml — informational only, for the `.cdb` header.
function _cl_pkgversion()
    p = joinpath(dirname(dirname(@__DIR__)), "Project.toml")
    isfile(p) || return "?"
    for line in eachline(p)
        startswith(strip(line), "version") && return strip(split(line, '=')[2], [' ', '"'])
    end
    return "?"
end

"""
    cl_coordinates(c::CircularComboR, b::CLBasis)
        -> (coeffs::Vector{SoergelPoly}, rest::CircularComboR)

`c` in the coordinates of basis `b`: `coeffs[i]` is the coefficient of the
`i`-th leaf, `rest` collects the terms whose leaf does NOT occur in `b`. A
nonempty `rest` is a finding (incomplete basis, or a leaf violating the degree
bound), not an error — the caller decides.
"""
function cl_coordinates(c::CircularComboR, b::CLBasis)
    index = cl_basis_index(b)
    coeffs = fill(zero(SoergelPoly), length(b.leaves))
    rest = CircularComboR()
    for (d, coeff) in pairs_of(c)
        k = circular_canonical_key(d)
        if haskey(index, k)
            coeffs[index[k]] = coeffs[index[k]] + coeff
        else
            rest = rest + coeff * CircularComboR(d)
        end
    end
    return (coeffs = coeffs, rest = rest)
end
