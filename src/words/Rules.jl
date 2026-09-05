# words layer — Rules.jl
#
# Local rewrite rules on `CircularWord`s. Everything here is applied CYCLICALLY:
# a rule matches a window of consecutive positions that may WRAP AROUND the end
# of the word. E.g. in (1,3,1,2) the braid window 121 can be read at positions
# (3,4,1) — letters 1,2,1 — and rewriting gives (…) → (2,3,2,1)-as-a-cycle, whose
# normal form is then taken: 1312 --121→212--> 2321.
#
# THE RULES (each listed as lhs <-> rhs, both directions available). The names
# are the pictures the rules draw, i ∈ {1,2,3}:
#
#   cap:
#       ii <-> ε        length -2 / +2      (reducing / increasing)
#                       the two i-strands join into one arc and close off
#   merge:
#       ii <-> i        length -1 / +1      (reducing / increasing)
#                       the two i-strands run into one trivalent
#   commutation:
#       13 <-> 31       length  0           (neutral)
#   braid:
#       121 <-> 212     length  0           (neutral)
#       232 <-> 323     length  0           (neutral)
#
# `dot` is the name of `i <-> ε`, one strand capped by a dot, and it is NOT in
# this set: with it every word would reach ε one letter at a time and the
# reachability question of `SequenceGraph.jl` would be vacuous.
#
# STEP SIGN (label on a graph edge for the move x -> y):
#       -1  strictly reducing  (length decreases)
#       +1  strictly increasing (length increases)
#        0  neutral            (length unchanged: commutation or braid)
#
# A `Move` records one concrete rewrite: which rule, at which cyclic start
# position, in which direction, and the resulting CircularWord together with the
# sign. `neighbors(w)` returns every distinct resulting word reachable in one
# step, each with the sign of the move that produced it. (If two different moves
# land on the same normal-form word, one representative per (word) is kept;
# `moves(w)` returns every raw move.)

# A single rewrite rule as an ordered pair of letter-patterns. Applying it
# forward replaces an occurrence of `from` by `to`; backward does the reverse.
struct Rule
    from::Vector{Letter}
    to::Vector{Letter}
    name::Symbol
end

# The base rules, listed in ONE direction. Both directions are generated in
# `all_directed_rules()`. `ε` is the empty pattern `Letter[]`.
const BASE_RULES = Rule[
    Rule([1,1], Letter[],   :cap1),
    Rule([2,2], Letter[],   :cap2),
    Rule([3,3], Letter[],   :cap3),
    Rule([1,1], [1],        :merge1),
    Rule([2,2], [2],        :merge2),
    Rule([3,3], [3],        :merge3),
    Rule([1,3], [3,1],      :comm_13),
    Rule([1,2,1], [2,1,2],  :braid_12),
    Rule([2,3,2], [3,2,3],  :braid_23),
]

# ---------------------------------------------------------------------------
# A SECOND, extended rule set: the base rules PLUS extra strictly-reducing
# rules (derivable from the base ones, but handy as direct shortcuts). These
# extra rules are listed ONE-DIRECTIONAL and are meant to be used as REDUCERS
# only — we do NOT run them backwards as increasing moves (see `down_moves`).
#
#   121212 -> ε     12121 -> 2     1212 -> 21      (lengths -6, -3, -2)
#   232323 -> ε     23232 -> 3     2323 -> 32      (same for {2,3})
#   131 -> 3        313 -> 1                       (each -2, = 131->311->3)
#
# NOTE: with cyclic normalization 1212->21 etc. are stated on the linear pattern;
# applied cyclically at every start position they cover all rotations.
const EXTRA_REDUCERS = Rule[
    Rule([1,2,1,2,1,2], Letter[], :red_121212),
    Rule([1,2,1,2,1],   [2],      :red_12121),
    Rule([1,2,1,2],     [2,1],    :red_1212),
    Rule([2,3,2,3,2,3], Letter[], :red_232323),
    Rule([2,3,2,3,2],   [3],      :red_23232),
    Rule([2,3,2,3],     [3,2],    :red_2323),
    Rule([1,3,1],       [3],      :red_131),
    Rule([3,1,3],       [1],      :red_313),
]

# the full extended set (base ∪ extra), used when we want every rule available
const EXTENDED_RULES = vcat(BASE_RULES, EXTRA_REDUCERS)

# ---------------------------------------------------------------------------
# The UPWARD move set: the strictly-increasing rules with a NONEMPTY source
# (no ε-exit). Listed in their INCREASING direction (from -> to, |to|>|from|),
# so `up_neighbors` keeps only the forward application. These are EXACTLY the
# nonempty-source increasers derived from EXTENDED_RULES:
#
#     1→11  2→22  3→33            (+1)
#     21→1212  32→2323            (+2)
#     2→12121  3→23232            (+4)
#     3→131  1→313                (+2)
#
# The ε-exits (ε→11, ε→22, ε→121212, …) are absent here — they are only allowed
# as the first step of the walk, injected via the seeds in `up_component`.
const UP_RULES = Rule[
    Rule([1],     [1,1],         :up_1_11),
    Rule([2],     [2,2],         :up_2_22),
    Rule([3],     [3,3],         :up_3_33),
    Rule([2,1],   [1,2,1,2],     :up_21_1212),
    Rule([3,2],   [2,3,2,3],     :up_32_2323),
    Rule([2],     [1,2,1,2,1],   :up_2_12121),
    Rule([3],     [2,3,2,3,2],   :up_3_23232),
    Rule([3],     [1,3,1],       :up_3_131),
    Rule([1],     [3,1,3],       :up_1_313),
]

"""
    sign_of(from_len, to_len) -> Int

The edge label for a move that changes length `from_len -> to_len`:
`-1` reducing, `+1` increasing, `0` neutral.
"""
sign_of(from_len::Int, to_len::Int) = sign(to_len - from_len)

# ---------------------------------------------------------------------------
# One recorded move
# ---------------------------------------------------------------------------

"""
    Move

One concrete application of a rule to a source `CircularWord`, producing `result`.
Fields:
  - `result` :: CircularWord   (already in normal form)
  - `sign`   :: Int            (-1, 0, +1; the graph edge label)
  - `rule`   :: Symbol         (which rule fired, e.g. :braid_12)
  - `pos`    :: Int            (cyclic start position 1..n of the matched window)
  - `forward`:: Bool           (rule applied in its listed direction?)
"""
struct Move
    result::CircularWord
    sign::Int
    rule::Symbol
    pos::Int
    forward::Bool
end

# ---------------------------------------------------------------------------
# Cyclic pattern matching + replacement
# ---------------------------------------------------------------------------

# Read `len` letters of `v` starting at (1-based) position `start`, wrapping
# cyclically. Returns the window as a Vector.
function _cyclic_window(v::Vector{Letter}, start::Int, len::Int)
    n = length(v)
    return Letter[ v[mod1(start + k, n)] for k in 0:(len-1) ]
end

# Replace the cyclic window of length `patlen` starting at `start` by `repl`.
# Returns a plain Vector (NOT yet normalized); caller wraps in CircularWord.
# We build the "rest" of the word (the letters NOT in the window), starting just
# after the window, then prepend the replacement. Because CircularWord normalizes
# by rotation, where exactly we cut is irrelevant to the resulting element.
function _cyclic_replace(v::Vector{Letter}, start::Int, patlen::Int, repl::Vector{Letter})
    n = length(v)
    rest = Letter[ v[mod1(start + patlen + k, n)] for k in 0:(n - patlen - 1) ]
    return vcat(repl, rest)
end

# ---------------------------------------------------------------------------
# Enumerate all one-step moves
# ---------------------------------------------------------------------------

"""
    moves(w) -> Vector{Move}

Every one-step rewrite applicable to `w`, in every rule, every direction, at
every cyclic position. May contain duplicates in `result` (different positions /
rules landing on the same element); use `neighbors` for the deduplicated view.

The empty word ε: only increasing moves apply (ε -> ii via the reverse of
`ii -> ε`, i.e. inserting `11`, `22`, `33`).
"""
function moves(w::CircularWord, rules::AbstractVector{Rule} = BASE_RULES)
    out = Move[]
    v = w.letters
    n = length(v)
    for rule in rules
        for (from, to, forward) in ((rule.from, rule.to, true), (rule.to, rule.from, false))
            patlen = length(from)
            if patlen == 0
                # inserting into the empty pattern: this is ε -> `to`. Only
                # meaningful when the whole word is empty (there is a single
                # "position"); for nonempty words an ε-window would let us insert
                # `to` between ANY two cyclic positions.
                _apply_insertions!(out, w, to, rule.name, forward)
                continue
            end
            n == 0 && continue          # nonempty pattern can't match ε
            # a length-n word has n cyclic start positions; but if patlen > n the
            # window would re-read letters — that is a genuine cyclic match only
            # when patlen <= n. (patlen<=3 here, so only tiny words are limited.)
            patlen > n && continue
            for start in 1:n
                _cyclic_window(v, start, patlen) == from || continue
                newv = _cyclic_replace(v, start, patlen, to)
                res = CircularWord(newv)
                push!(out, Move(res, sign_of(n, length(res)), rule.name, start, forward))
            end
        end
    end
    return out
end

# Insertion of pattern `to` at every gap (the reverse of an `ii -> ε` rule).
# For a word of length n there are n cyclic gaps to insert at; for ε there is one.
function _apply_insertions!(out::Vector{Move}, w::CircularWord, to::Vector{Letter},
                            name::Symbol, forward::Bool)
    v = w.letters
    n = length(v)
    if n == 0
        res = CircularWord(copy(to))
        push!(out, Move(res, sign_of(0, length(res)), name, 1, forward))
        return
    end
    for gap in 1:n
        # insert `to` after position `gap`
        newv = vcat(v[1:gap], to, v[gap+1:end])
        res = CircularWord(newv)
        push!(out, Move(res, sign_of(n, length(res)), name, gap, forward))
    end
end

# ---------------------------------------------------------------------------
# Full-word increasing expansions (for the "expand then reduce" search)
# ---------------------------------------------------------------------------
#
# The BRAID-EXPANDERS: the ONLY increasing moves used for full-word expansion.
# These are NOT the dot increasers (no ε->ii, no i->ii). Each replaces a short
# alternating block by a longer alternating block within {1,2} or within {2,3}:
#
#   {1,2}:  12 -> 2121   21 -> 1212   1 -> 21212   2 -> 12121
#   {2,3}:  23 -> 3232   32 -> 2323   3 -> 23232   2 -> 32323
#
# (So the letter 2 has two expanders — one on each side; 1 and 3 have one each.)
# A single letter grows to a 5-letter alternating word (length +4); a 2-letter
# block grows to 4 letters (length +2).
# The four +2 braid-expanders (2-letter block → 4-letter block). These alone can
# only tile EVEN-length words.
const BRAID_EXPANDERS_2 = Tuple{Vector{Letter},Vector{Letter}}[
    ([1,2], [2,1,2,1]),  ([2,1], [1,2,1,2]),
    ([2,3], [3,2,3,2]),  ([3,2], [2,3,2,3]),
]

# The four single-letter braid-expanders (letter → 5-letter block, +4). Including
# these lets odd-length words tile too.
const BRAID_EXPANDERS_1 = Tuple{Vector{Letter},Vector{Letter}}[
    ([1], [2,1,2,1,2]), ([2], [1,2,1,2,1]),
    ([3], [2,3,2,3,2]), ([2], [3,2,3,2,3]),
]

# The full set of 8 (used by default in increasing_expansions et al.).
const BRAID_EXPANDERS = vcat(BRAID_EXPANDERS_2, BRAID_EXPANDERS_1)
const _INCREASERS = BRAID_EXPANDERS

# All increasing expansions of the LINEAR word `v` tiled left-to-right into
# consecutive blocks, each block = the LHS of some increaser, expanded to its RHS.
# EVERY position is covered (disjoint tiling). Recursive: peel a matching
# increaser off the front, recurse on the rest.
function _linear_expansions(v::Vector{Letter}, expanders = _INCREASERS)
    isempty(v) && return [Letter[]]
    out = Vector{Letter}[]
    for (from, to) in expanders
        L = length(from)
        L <= length(v) || continue
        v[1:L] == from || continue
        for rest in _linear_expansions(v[L+1:end], expanders)
            push!(out, vcat(to, rest))
        end
    end
    return out
end

"""
    increasing_expansions(w) -> Vector{CircularWord}

All circular words obtainable from `w` by choosing a **disjoint cyclic tiling** of
its positions into consecutive blocks, each block matching the LHS of an
increasing extended rule, and expanding every block (so EVERY letter takes part
in an increasing move). All tilings and all per-block rule choices are tried;
cyclic (wrapping) tilings are realised by tiling each rotation of `w` linearly.
Results are normal-form `CircularWord`s (deduplicated). No valid full tiling ⇒
empty result.

Example: `increasing_expansions(CircularWord([1,2,2,1]))` contains `21211212`
(from 12->2121, 21->1212), which reduces to ε.
"""
function increasing_expansions(w::CircularWord, expanders = _INCREASERS)
    isempty(w) && return CircularWord[]
    v = w.letters
    n = length(v)
    seen = Set{CircularWord}()
    for k in 0:(n-1)                       # each rotation, to realise wrapping tilings
        rot = vcat(v[k+1:end], v[1:k])
        for ex in _linear_expansions(rot, expanders)
            push!(seen, CircularWord(ex))
        end
    end
    return collect(seen)
end

# Like `_linear_expansions` but records the tiling as a list of (block, expansion)
# pairs so we can SHOW how the word was expanded.
function _linear_tilings(v::Vector{Letter}, expanders = _INCREASERS)
    isempty(v) && return [Tuple{Vector{Letter},Vector{Letter}}[]]
    out = Vector{Tuple{Vector{Letter},Vector{Letter}}}[]
    for (from, to) in expanders
        L = length(from)
        L <= length(v) || continue
        v[1:L] == from || continue
        for rest in _linear_tilings(v[L+1:end], expanders)
            push!(out, vcat([(copy(from), copy(to))], rest))
        end
    end
    return out
end

"""
    has_full_expansion(w, expanders = BRAID_EXPANDERS) -> Bool

True iff `w` admits at least one full disjoint cyclic tiling by `expanders`
(every position covered) — i.e. `w` "can be expanded all around". Pass
`BRAID_EXPANDERS_2` to allow only the four +2 (2-letter-block) expanders.
"""
function has_full_expansion(w::CircularWord, expanders = _INCREASERS)
    isempty(w) && return false
    v = w.letters
    for k in 0:(length(v)-1)
        rot = vcat(v[k+1:end], v[1:k])
        isempty(_linear_expansions(rot, expanders)) || return true
    end
    return false
end

"""
    a_full_expansion(w, expanders = BRAID_EXPANDERS) -> Union{Nothing, NamedTuple}

One witnessing full expansion of `w`, or `nothing` if none exists. The returned
named tuple has `rotation` (the rotated linear word that was tiled), `blocks` (the
list of (block, expansion) pairs, in order), and `result` (the expanded
CircularWord). Use for the "print how one would increase this" option.
"""
function a_full_expansion(w::CircularWord, expanders = _INCREASERS)
    isempty(w) && return nothing
    v = w.letters
    for k in 0:(length(v)-1)
        rot = vcat(v[k+1:end], v[1:k])
        tilings = _linear_tilings(rot, expanders)
        if !isempty(tilings)
            blocks = tilings[1]
            expanded = reduce(vcat, (to for (_, to) in blocks); init = Letter[])
            return (rotation = rot, blocks = blocks, result = CircularWord(expanded))
        end
    end
    return nothing
end

"""
    all_full_expansions(w, expanders = BRAID_EXPANDERS) -> Vector{NamedTuple}

Every witnessing full expansion of `w` (over all rotations and tilings), each as a
`(rotation, blocks, result)` named tuple — deduplicated by `result`. Empty if `w`
has no full tiling. Used when we must try each expansion to find one that then
reduces to ε.
"""
function all_full_expansions(w::CircularWord, expanders = _INCREASERS)
    isempty(w) && return NamedTuple[]
    v = w.letters
    seen = Set{CircularWord}()
    out = NamedTuple[]
    for k in 0:(length(v)-1)
        rot = vcat(v[k+1:end], v[1:k])
        for blocks in _linear_tilings(rot, expanders)
            expanded = reduce(vcat, (to for (_, to) in blocks); init = Letter[])
            res = CircularWord(expanded)
            if !(res in seen)
                push!(seen, res)
                push!(out, (rotation = rot, blocks = blocks, result = res))
            end
        end
    end
    return out
end

"""
    uses_all_letters(w) -> Bool

True iff `w` contains at least one 1, one 2, and one 3.
"""
uses_all_letters(w::CircularWord) =
    (1 in w.letters) && (2 in w.letters) && (3 in w.letters)

"""
    neighbors(w) -> Vector{Tuple{CircularWord,Int}}

The distinct one-step neighbors of `w`, each paired with the sign of a move that
reaches it. If several moves reach the same neighbor with different signs (can
happen at the boundary), the sign of least length-change magnitude is kept
(prefers 0 over ±1); ties keep the first seen. Self-loops (result == w) are
dropped.
"""
function neighbors(w::CircularWord, rules::AbstractVector{Rule} = BASE_RULES)
    best = Dict{CircularWord,Int}()
    for m in moves(w, rules)
        m.result == w && continue
        if haskey(best, m.result)
            abs(m.sign) < abs(best[m.result]) && (best[m.result] = m.sign)
        else
            best[m.result] = m.sign
        end
    end
    return [(k, v) for (k, v) in best]
end

"""
    down_neighbors(w) -> Vector{Tuple{CircularWord,Int}}

The distinct one-step neighbors reachable from `w` using only **non-increasing**
moves under the EXTENDED rule set:
  - every rule applied in the length-**reducing** direction (all of BASE ∪ EXTRA:
    ii→ε, ii→i, 131→3, 1212→21, 12121→2, 121212→ε, …), and
  - the length-**neutral** braid/commutation moves (13↔31, 121↔212, 232↔323),
    available in both directions (they keep length).

Increasing moves are excluded. This is the move set for the "can we get DOWN to
ε?" graph (see `reaches_empty_down` in SequenceGraph.jl). Each neighbor is paired
with its sign (−1 reducing, 0 neutral). Self-loops dropped.
"""
function down_neighbors(w::CircularWord)
    best = Dict{CircularWord,Int}()
    for m in moves(w, EXTENDED_RULES)
        m.result == w && continue
        m.sign > 0 && continue          # drop strictly increasing moves
        if haskey(best, m.result)
            abs(m.sign) < abs(best[m.result]) && (best[m.result] = m.sign)
        else
            best[m.result] = m.sign
        end
    end
    return [(k, v) for (k, v) in best]
end

"""
    up_neighbors(w) -> Vector{CircularWord}

The distinct one-step neighbors reachable from `w` using only strictly
**increasing** moves of the EXTENDED rule set, EXCLUDING any ε-exit — i.e. every
increasing extended-rule direction whose source pattern is nonempty:

    1→11  2→22  3→33            (dot_single reversed, +1)
    21→1212  32→2323           (red_1212 / red_2323 reversed, +2)
    2→12121  3→23232           (red_12121 / red_23232 reversed, +4)
    3→131  1→313               (red_131 / red_313 reversed, +2)

The ε-exits (ε→11, ε→22, ε→121212, and the dropped ε→33, ε→232323) are NOT here —
they are only ever used as the FIRST step of the upward walk (see
`up_component`). Neutral braid/commutation moves are also excluded (they don't
increase length). Self-loops dropped. Result is just the words (all moves +1).

Implemented from the EXPLICIT list `UP_RULES` of nonempty-source increasers (so no
ε-pattern insertion can sneak in), applied cyclically in both start positions.
"""
function up_neighbors(w::CircularWord)
    isempty(w) && return CircularWord[]     # ε-exits are handled by the seeds only
    out = Set{CircularWord}()
    for m in moves(w, UP_RULES)             # only nonempty-source increasers
        m.forward || continue               # listed in the increasing direction
        m.result == w && continue
        push!(out, m.result)
    end
    return collect(out)
end
