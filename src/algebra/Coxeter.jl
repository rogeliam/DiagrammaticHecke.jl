# algebra/Coxeter.jl — the Coxeter group of type A₃, from CoxeterGroups.jl.
#
# Elements are `CoxEltMin`s of the minimal-roots implementation; the primitives
# the package uses are `one`, `right_multiply`, `coxeter_length`,
# `is_right_descent`, `short_lex` (the canonical reduced word) and
# `coxeter_matrix`. The methods below extend CoxeterGroups.jl: a constructor
# from a type name, printing, indexing a group by a word, and the enumeration
# of a finite group.

import CoxeterGroups
using CoxeterGroups: CoxGrp, CoxGrpMin, CoxEltMin, coxeter_group_min,
                     coxeter_matrix_from_group_type, generators, rank,
                     is_right_descent, length as coxeter_length,
                     short_lex, right_multiply, coxeter_matrix

"""
    coxeter_group_min(group_type::AbstractString) -> (W, generators)

The Coxeter group of the named type (`"A3"`, `"B2"`, `"I2(5)"`, products such
as `"A3 x B4"`), via `coxeter_matrix_from_group_type`.
"""
CoxeterGroups.coxeter_group_min(group_type::AbstractString) =
    coxeter_group_min(coxeter_matrix_from_group_type(group_type))

function Base.show(io::IO, W::CoxGrpMin)
    n = size(W.coxeter_matrix, 1)
    print(io, "Coxeter group of rank ", n, " (Coxeter matrix ",
          replace(string(W.coxeter_matrix), "\n" => " "), ")")
end

"An element prints as its short-lex word, `s1s2s1`, and the identity as `id`."
Base.show(io::IO, x::CoxEltMin) =
    isempty(x.word) ? print(io, "id") : print(io, join(("s$(Int(i))" for i in x.word)))

"""
    W[i₁, i₂, …] -> CoxEltMin

The element of `W` spelled by the word in the generators, `W[1, 2, 1] = s₁s₂s₁`;
the word need not be reduced. `W[]` is the identity.
"""
function Base.getindex(W::CoxGrpMin, word::Integer...)
    w = one(W)
    for i in word
        1 <= i <= rank(W) || error("generator index $i out of range 1:$(rank(W))")
        w = right_multiply(w, Int(i))
    end
    return w
end

"""
    enumerate_whole_group(W::CoxGrp) -> Vector

All elements of the finite group `W`, the identity first and in non-decreasing
Coxeter length: one breadth-first pass over the right Cayley graph, generators
in increasing order.
"""
function enumerate_whole_group(W::CoxGrp)
    e = one(W)
    queue = [e]
    seen = Set([e])
    gens = generators(W)
    i = 1
    while i <= length(queue)
        for s in gens
            ws = queue[i] * s
            ws in seen && continue
            push!(seen, ws)
            push!(queue, ws)
        end
        i += 1
    end
    return queue
end

const _A3_GROUP = Ref{Any}(nothing)

"""
    coxeter_group_a3() -> (W, generators)

The Coxeter group of type A₃ together with its three simple reflections, built
once and cached.
"""
function coxeter_group_a3()
    _A3_GROUP[] === nothing && (_A3_GROUP[] = coxeter_group_min("A3"))
    return _A3_GROUP[]
end
