# CircularLightLeaves.jl — light leaves & double leaves as circular diagrams
#
# Built on the existing LL/DL constructors in morphism/LightLeaves.jl, but
# immediately convert the result to CircularMorphismGraph: build the diagram
# with the normal tooling, then interpret it as a circular diagram.

"""
    circular_light_leaf(word, e) -> CircularMorphismGraph

Build [`light_leaf`](@ref) and interpret it as `CircularMorphismGraph`.
"""
function circular_light_leaf(word::Vector{Int}, e::Vector{Int})
    return circular_morphism(light_leaf(word, e))
end

"""
    circular_double_leaf(xword, e, yword, f) -> CircularMorphismGraph

Build [`double_leaf`](@ref) and interpret it as `CircularMorphismGraph`.
"""
function circular_double_leaf(xword::Vector{Int}, e::Vector{Int},
                         yword::Vector{Int}, f::Vector{Int})
    return circular_morphism(double_leaf(xword, e, yword, f))
end

"""
    circular_light_leaves(xword, zword; degree = nothing)
        -> Vector{Tuple{Vector{Int}, CircularMorphismGraph}}

Like [`light_leaves`](@ref), but with `CircularMorphismGraph` instead of `MorphismGraph`.
"""
function circular_light_leaves(xword::Vector{Int}, zword::Vector{Int};
                          degree::Union{Nothing,Int} = nothing)
    return [(e, circular_morphism(m)) for (e, m) in light_leaves(xword, zword; degree = degree)]
end

"""
    circular_double_leaves(xword, yword; degree = nothing)
        -> Vector{@NamedTuple{e::Vector{Int}, f::Vector{Int}, z::Vector{Int},
                              degree::Int, morphism::CircularMorphismGraph}}

Like [`double_leaves`](@ref), but with `CircularMorphismGraph` instead of `MorphismGraph`.
"""
function circular_double_leaves(xword::Vector{Int}, yword::Vector{Int};
                           degree::Union{Nothing,Int} = nothing)
    raw = double_leaves(xword, yword; degree = degree)
    T = @NamedTuple{e::Vector{Int}, f::Vector{Int}, z::Vector{Int},
                    degree::Int, morphism::CircularMorphismGraph}
    out = T[]
    for r in raw
        push!(out, (e = r.e, f = r.f, z = r.z, degree = r.degree,
                    morphism = circular_morphism(r.morphism)))
    end
    return out
end

"""
    _parse_meta(meta::AbstractString) -> @NamedTuple{e::Vector{Int}, f::Vector{Int},
                                                     z::Vector{Int}, degree::Int}

Parses a `.wgm` meta line of the form
`# e=1,1,0 f=1,0,1 z=1 degree=2` (or `z=eps`).
"""
function _parse_meta(meta::AbstractString)
    s = strip(meta)
    startswith(s, "#") && (s = strip(s[2:end]))
    e = Int[]; f = Int[]; z = Int[]; degree = 0
    for part in split(s)
        kv = split(part, '='; limit = 2)
        length(kv) != 2 && continue
        key, val = kv[1], kv[2]
        if key == "e" || key == "f"
            vals = val == "" ? Int[] : parse.(Int, split(val, ','))
            key == "e" ? (e = vals) : (f = vals)
        elseif key == "z"
            z = val == "eps" ? Int[] : parse.(Int, collect(val))
        elseif key == "degree"
            degree = parse(Int, val)
        end
    end
    return (; e, f, z, degree)
end

"""
    circular_load_double_leaves(path::AbstractString)
        -> Vector{@NamedTuple{e::Vector{Int}, f::Vector{Int}, z::Vector{Int},
                              degree::Int, morphism::CircularMorphismGraph}}

Loads saved double leaves from a `.wgm` file (MorphismEntry format) and
converts each morphism to `CircularMorphismGraph`.

The metadata `e`, `f`, `z`, `degree` is read from the file's comment lines
(`# e=... f=... z=... degree=...`). If no metadata line is present, these
fields stay empty resp. 0.
"""
function circular_load_double_leaves(path::AbstractString)
    entries = load_morphisms(path)
    T = @NamedTuple{e::Vector{Int}, f::Vector{Int}, z::Vector{Int},
                    degree::Int, morphism::CircularMorphismGraph}
    out = T[]
    for entry in entries
        meta = _parse_meta(entry.meta)
        push!(out, (e = meta.e, f = meta.f, z = meta.z, degree = meta.degree,
                    morphism = circular_morphism(entry.morphism)))
    end
    return out
end

"""
    circular_reduce_double_leaf(dl::CircularMorphismGraph; rules = CIRCULAR_RULES,
                           all_matches = false) -> CircularComboR

Fully reduce a single double leaf (as a circular diagram).
"""
function circular_reduce_double_leaf(dl::CircularMorphismGraph; rules = CIRCULAR_RULES,
                                all_matches::Bool = false)
    return reduce_circular_full(dl.graph; rules = rules, all_matches = all_matches)
end

"""
    circular_reduce_double_leaf(dl::MorphismGraph; kwargs...) -> CircularComboR

Convenience entry point: convert a normal morphism to circular first.
"""
function circular_reduce_double_leaf(dl::MorphismGraph; kwargs...)
    return circular_reduce_double_leaf(circular_morphism(dl); kwargs...)
end
