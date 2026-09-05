# words layer — CircularWord.jl
#
# A CircularWord is a word over the alphabet {1,2,3}, considered up to CYCLIC
# rotation:
#
#     (i1, i2, …, in)  ≡  (i2, …, in, i1).
#
# This layer is purely combinatorial: the letters are the symbols 1, 2, 3 and no
# group theory enters. The connection to type-A₃ Coxeter combinatorics and Soergel
# diagrams is made in the layers above.
#
# NORMAL FORM. Since all n rotations denote the same element, we pick a canonical
# representative: the lexicographically smallest rotation. Ties (rotations that
# are equal as tuples, e.g. the two even rotations of (1,2,1,2)) collapse to the
# same tuple, so `min` over all rotations is well defined. The empty word () is a
# valid element, length 0, and is its own normal form.
#
# LENGTH is simply the number of letters. Note: with the rules in Rules.jl the
# length is NOT a rotation invariant only — it is honestly the letter count of
# the normal form (rotation does not change length).

const Letter = Int   # ∈ {1,2,3}

"""
    CircularWord(letters)

A word over {1,2,3} taken up to cyclic rotation, stored in its canonical
(lexicographically minimal rotation) normal form.

Construct from any iterable of integers in 1:3, e.g. `CircularWord([1,3,1,2])`
or `CircularWord((1,3,1,2))`. The empty word is `CircularWord(Int[])`.

Two circular words are `==` (and hash equal) iff they are the same element, i.e.
rotations of each other — this is guaranteed because both are stored in normal
form.
"""
struct CircularWord
    letters::Vector{Letter}   # ALWAYS the normal form (canonical rotation)
    len::Int                  # cached length == length(letters)
    str::String               # cached compact form, e.g. "1213" (or "ε" for the
                              # empty word). Handy for printing and as a dict key.

    function CircularWord(raw::AbstractVector{<:Integer})
        for x in raw
            (1 <= x <= 3) || throw(ArgumentError("letters must be in 1:3, got $x"))
        end
        canon = _canonical_rotation(collect(Letter, raw))
        new(canon, length(canon), isempty(canon) ? "ε" : join(canon, ""))
    end
end

CircularWord(raw::Tuple) = CircularWord(collect(Letter, raw))
CircularWord(raw::Integer...) = CircularWord(collect(Letter, raw))

# ---------------------------------------------------------------------------
# Normal form: lexicographically minimal rotation.
# ---------------------------------------------------------------------------

"""
    _canonical_rotation(v) -> Vector

Return the lexicographically smallest cyclic rotation of `v`. `v` is not
mutated. The empty vector is its own canonical form.

(Straightforward O(n²) scan — n is tiny here. Could be Booth's O(n) algorithm
later if it ever matters.)
"""
function _canonical_rotation(v::Vector{Letter})
    n = length(v)
    n <= 1 && return copy(v)
    best = copy(v)
    for k in 1:(n-1)
        rot = vcat(v[k+1:end], v[1:k])
        if _lex_less(rot, best)
            best = rot
        end
    end
    return best
end

# lexicographic <, on equal-length vectors (all rotations share length n)
function _lex_less(a::AbstractVector, b::AbstractVector)
    @inbounds for i in eachindex(a)
        a[i] != b[i] && return a[i] < b[i]
    end
    return false   # equal
end

# ---------------------------------------------------------------------------
# Basic interface
# ---------------------------------------------------------------------------

Base.length(w::CircularWord) = w.len
Base.isempty(w::CircularWord) = w.len == 0

# element identity == rotation-equivalence; both are stored normalized already
Base.:(==)(a::CircularWord, b::CircularWord) = a.letters == b.letters
Base.hash(w::CircularWord, h::UInt) = hash(w.letters, hash(:CircularWord, h))

# a stable total order (by length, then lexicographically on the normal form) —
# handy for "take the lowest length not yet visited" style loops and for sorting.
function Base.isless(a::CircularWord, b::CircularWord)
    la, lb = length(a), length(b)
    la != lb && return la < lb
    return _lex_less(a.letters, b.letters)
end

# the underlying letters of the canonical representative (do not mutate)
letters(w::CircularWord) = w.letters

# the empty circular word
const EMPTY = CircularWord(Letter[])

"""
    rotations(w) -> Vector{Vector{Int}}

All `n` cyclic rotations of the normal-form representative (with repeats if the
word has cyclic symmetry). Mostly for tests / inspection.
"""
function rotations(w::CircularWord)
    v = w.letters
    n = length(v)
    n == 0 && return [Letter[]]
    return [vcat(@view(v[k+1:end]), @view(v[1:k])) for k in 0:(n-1)]
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

# show as e.g. CircularWord(1213) or CircularWord(ε), using the cached string.
Base.show(io::IO, w::CircularWord) = print(io, "CircularWord(", w.str, ")")

# the cached compact string form, e.g. "1213" (or "ε" for the empty word) — used
# by the graph/data code for labels / keys when a short human-readable name is wanted.
compact(w::CircularWord) = w.str
