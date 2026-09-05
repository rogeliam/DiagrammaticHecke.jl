# circular/CircularDecoratedIO.jl — save/load a CircularDecorated (graph + region labels) as
# versioned, git-diffable text; on top of that the BASIS FILE `.cdb`.
#
# FILE EXTENSIONS vs. MAGIC LINE. The file extensions are `.cd`/`.cdb`/`.cdc`, but
# the magic lines inside are `#fd` / `#fdb` / `#fwg` and the constants are
# `FD_VERSION` / `FDB_VERSION` — CONTENT and extension are decoupled, so that
# already-written cache and fixture files stay readable even if the extension
# changes.
#
# WHY. Computing the circular-L basis of a word pair costs minutes to hours (every
# double leaf once through `reduce_to_circular_leave`). Across the seven
# composition stages, and across sessions, the same basis is needed again and
# again — so it belongs on disk. What gets saved is the REPRESENTATIVE of each
# leaf (graph plus labels), not the canonical key: the key is a pure function
# of the representative (`circular_canonical_key`) and is recomputed on load, so a
# later change to the key construction never silently carries old keys
# forward.
#
# THE SWITCH SETTING LIVES IN THE HEADER. A basis computed with zamo disabled
# (`CIRCULAR_ZAMO_ENABLED[] = false`) is a DIFFERENT basis from the normal one; the
# loader compares the header line `zamo on|off` against the requested setting
# and discards the file on mismatch (the caller then recomputes). The same
# holds for the package version: it is recorded in the header but only
# REPORTED, never enforced — whether an older version still fits is the
# caller's call.
#
# FORMAT — `.cd` (one decorated graph):
#
#     #fd 1
#     labels 1 | a1 | 2,0,0:3/2
#     outer 1
#     #fwg 1
#     word 121321
#     …                       (verbatim `.fwg` body, see CircularGraphIO.jl)
#
# `labels` is the label sequence in region order, separated by ` | ` (an empty
# sequence => the line is omitted). A polynomial is `0`, `1`, or a list
# `a,b,c:num/den` (monomial α₁^a α₂^b α₃^c), separated by `;` — sorted so the
# file diffs stably.
#
# FORMAT — `.cdb` (one basis): a header, then `.cd` entries separated by `---`.
#
#     #fdb 1
#     pkgversion 0.10.33
#     zamo off
#     x 121321
#     y 321323
#     ndl 204
#     ---
#     deg 0
#     #fd 1
#     …

# Word <-> token as in GraphIO.jl, but here over `Vector{Int}` (the word pairs
# `x`/`y` are letter lists, not `CircularWord`s).
_fdb_word_token(w::Vector{Int}) = isempty(w) ? "eps" : join(w)
_fdb_token_word(tok::AbstractString) =
    tok == "eps" ? Int[] : [Int(c - '0') for c in tok]

const FD_VERSION  = 1
const FDB_VERSION = 1

# ---- SoergelPoly <-> Text ----------------------------------------------------

"""
    soergelpoly_to_string(p::SoergelPoly) -> String

`p` as stable text: `0` for the zero polynomial, `1` for the unit, otherwise
the exponent-sorted list `a,b,c:num/den` (monomial α₁^a α₂^b α₃^c), separated
by `;`. Inverse: [`soergelpoly_from_string`](@ref).
"""
function soergelpoly_to_string(p::SoergelPoly)
    iszero(p) && return "0"
    isone(p) && return "1"
    t = poly_terms(p)
    es = sort!(collect(keys(t)))
    return join(("$(e[1]),$(e[2]),$(e[3]):$(numerator(t[e]))/$(denominator(t[e]))"
                 for e in es), ";")
end

"""
    soergelpoly_from_string(s::AbstractString) -> SoergelPoly

Inverse of [`soergelpoly_to_string`](@ref).
"""
function soergelpoly_from_string(s::AbstractString)
    t = strip(s)
    t == "0" && return zero(SoergelPoly)
    t == "1" && return one(SoergelPoly)
    d = Dict{NTuple{3,Int},Rational{BigInt}}()
    for part in split(t, ';')
        isempty(strip(part)) && continue
        mono, coeff = split(part, ':')
        a, b, c = (parse(Int, x) for x in split(mono, ','))
        num, den = split(coeff, '/')
        d[(a, b, c)] = Rational{BigInt}(parse(BigInt, num), parse(BigInt, den))
    end
    return poly_from_terms(d)
end

# ---- CircularDecorated <-> Text ---------------------------------------------------

"""
    circulardecorated_to_string(d::CircularDecorated) -> String

`d` in the `.cd` format (see file header): a labels line, an `outer` line, and
below that the `.fwg` body of the graph, verbatim. Inverse:
[`circulardecorated_from_string`](@ref); round-trips exactly.
"""
function circulardecorated_to_string(d::CircularDecorated)
    io = IOBuffer()
    println(io, "#fd $FD_VERSION")
    isempty(d.region_labels) ||
        println(io, "labels ", join((soergelpoly_to_string(f) for f in d.region_labels), " | "))
    println(io, "outer ", soergelpoly_to_string(d.outer_label))
    print(io, circulargraph_to_string(d.graph))
    return String(take!(io))
end

"""
    circulardecorated_from_string(s::AbstractString) -> CircularDecorated

Inverse of [`circulardecorated_to_string`](@ref). Throws if the number of labels
does not match the region count of the graph read (`CircularDecorated` checks that).
"""
function circulardecorated_from_string(s::AbstractString)
    labels = SoergelPoly[]
    outer  = one(SoergelPoly)
    body   = IOBuffer()
    for raw in split(s, '\n')
        line = strip(raw)
        if startswith(line, "#fd")
            v = parse(Int, split(line)[2])
            v == FD_VERSION || error("unsupported .cd version $v (expected $FD_VERSION)")
        elseif startswith(line, "labels ")
            labels = [soergelpoly_from_string(t) for t in split(line[8:end], '|')]
        elseif startswith(line, "outer ")
            outer = soergelpoly_from_string(line[7:end])
        else
            println(body, raw)
        end
    end
    g = circulargraph_from_string(String(take!(body)))
    return CircularDecorated(g, labels, outer)
end

"Writes `d` as a `.cd` text file to `path`."
save_circulardecorated(path::AbstractString, d::CircularDecorated) =
    open(io -> write(io, circulardecorated_to_string(d)), path, "w")

"Reads a `CircularDecorated` from a file written by `save_circulardecorated`."
load_circulardecorated(path::AbstractString) = circulardecorated_from_string(read(path, String))

# ---- the basis file `.cdb` ----------------------------------------------------

"""
    CLBasis

The circular-L "basis" of ONE word pair, as it stands in a `.cdb` file: whether
it is actually a basis is still open — more precisely a `leaf_set`.

* `x`, `y`      — the word pair (`leaves` are morphisms `x → y`);
* `leaves`      — one REPRESENTATIVE per leaf, in stable order;
* `degrees`     — `circular_degree` per leaf (same order);
* `zamo`        — the switch setting the run was computed under (`CIRCULAR_ZAMO_ENABLED[]`);
* `pkgversion`  — the package version of the run (informational only, not enforced);
* `ndl`         — number of double leaves worked through (provenance).

The canonical key is NOT stored, but recomputed from the representative
([`cl_basis_index`](@ref)).
"""
struct CLBasis
    x::Vector{Int}
    y::Vector{Int}
    leaves::Vector{CircularDecorated}
    degrees::Vector{Int}
    zamo::Bool
    pkgversion::String
    ndl::Int
end

"""
    cl_basis_index(b::CLBasis) -> Dict{Any,Int}

`circular_canonical_key(leaf) => position` — the coordinate assignment used to read
off a reduced sum in this basis. Freshly computed from the representatives,
never read from the file.
"""
cl_basis_index(b::CLBasis) =
    Dict{Any,Int}(circular_canonical_key(d) => i for (i, d) in enumerate(b.leaves))

Base.length(b::CLBasis) = length(b.leaves)

Base.show(io::IO, b::CLBasis) =
    print(io, "CLBasis(", _dlb_tok(b.x), " → ", _dlb_tok(b.y), ", ",
          length(b.leaves), " leaf/leaves, zamo ", b.zamo ? "on" : "off", ")")

"""
    cl_basis_to_string(b::CLBasis) -> String

`b` in the `.cdb` format (see file header).
"""
function cl_basis_to_string(b::CLBasis)
    io = IOBuffer()
    println(io, "#fdb $FDB_VERSION")
    println(io, "pkgversion ", b.pkgversion)
    println(io, "zamo ", b.zamo ? "on" : "off")
    println(io, "x ", _fdb_word_token(b.x))
    println(io, "y ", _fdb_word_token(b.y))
    println(io, "ndl ", b.ndl)
    for (i, d) in enumerate(b.leaves)
        println(io, "---")
        println(io, "deg ", b.degrees[i])
        print(io, circulardecorated_to_string(d))
    end
    return String(take!(io))
end

"""
    cl_basis_from_string(s::AbstractString) -> CLBasis

Inverse of [`cl_basis_to_string`](@ref).
"""
function cl_basis_from_string(s::AbstractString)
    chunks = split(s, "\n---\n")
    head = chunks[1]
    x = Int[]; y = Int[]; zamo = true; pkgver = "?"; ndl = 0
    for raw in split(head, '\n')
        line = strip(raw)
        if startswith(line, "#fdb")
            v = parse(Int, split(line)[2])
            v == FDB_VERSION || error("unsupported .cdb version $v (expected $FDB_VERSION)")
        elseif startswith(line, "pkgversion ")
            pkgver = strip(line[12:end])
        elseif startswith(line, "zamo ")
            zamo = strip(line[6:end]) == "on"
        elseif startswith(line, "x ")
            x = _fdb_token_word(strip(line[3:end]))
        elseif startswith(line, "y ")
            y = _fdb_token_word(strip(line[3:end]))
        elseif startswith(line, "ndl ")
            ndl = parse(Int, strip(line[5:end]))
        end
    end
    leaves = CircularDecorated[]; degrees = Int[]
    for ch in chunks[2:end]
        strip(ch) == "" && continue
        deg = nothing
        body = IOBuffer()
        for raw in split(ch, '\n')
            line = strip(raw)
            if startswith(line, "deg ")
                deg = parse(Int, strip(line[5:end]))
            else
                println(body, raw)
            end
        end
        push!(leaves, circulardecorated_from_string(String(take!(body))))
        push!(degrees, deg === nothing ? circular_degree(leaves[end]) : deg)
    end
    return CLBasis(x, y, leaves, degrees, zamo, pkgver, ndl)
end

"Writes `b` as a `.cdb` text file to `path` (creates the directory)."
function save_cl_basis(path::AbstractString, b::CLBasis)
    mkpath(dirname(path))
    open(io -> write(io, cl_basis_to_string(b)), path, "w")
    return path
end

"""
    load_cl_basis(path; zamo = nothing) -> Union{CLBasis,Nothing}

Reads a `.cdb` file. Returns `nothing` if the file is missing OR (when `zamo`
is set) its switch setting does not match the requested one — the caller then
recomputes. A format error, in contrast, throws.
"""
function load_cl_basis(path::AbstractString; zamo::Union{Nothing,Bool} = nothing)
    isfile(path) || return nothing
    b = cl_basis_from_string(read(path, String))
    (zamo === nothing || b.zamo == zamo) || return nothing
    return b
end
