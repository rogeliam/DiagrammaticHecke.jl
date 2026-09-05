# CORE — the shared diagram/word layer under the circular pipeline.
# ColoredReduction.jl
#
# The "no simplification within one expanded block" condition.
#
# When a candidate word is expanded all-around into BLOCKS (e.g. 32->2323 |
# 21->1212 | …), each output letter belongs to exactly one block. The extra rule:
# on the way back DOWN to ε, a reducing/neutral move may fire ONLY IF
# its matched window spans ≥ 2 DISTINCT blocks. A move whose whole window sits
# inside a single block is FORBIDDEN — you may not undo (part of) a block using
# only that block's own letters. Simplifications must cross block boundaries.
#
# To decide that, every letter must carry its block tag through the rewrites. So
# we work with COLORED words: a cyclic word plus a per-letter block id, kept in the
# SAME normal form as CircularWord (lex-min rotation) with the tag vector rotated
# alongside. `down_neighbors_colored` returns the legal (cross-block) non-increasing
# moves; `colored_down_path` gives a witnessing path down to ε (or `nothing`).
#
# TAG PROPAGATION. We only ever use NON-INCREASING moves here (reducers + neutral
# braid/commutation). For each such rule we record, for every OUTPUT letter, WHICH
# input-window position it inherits its tag from (`_TAG_MAP`). Letters outside the
# window keep their tags unchanged. This is only consulted for surviving letters;
# the legality test itself just looks at the window's input tags.

# A colored circular word: letters in {1,2,3} plus a block tag per letter, stored
# in the canonical (lex-min) rotation of the LETTERS, with tags permuted to match.
# Ties in the letter rotation are broken to keep a deterministic (letters,tags).
struct ColoredWord
    letters::Vector{Letter}
    blocks::Vector{Int}
    str::String                     # e.g. "12·12" not needed; we cache plain letters
    function ColoredWord(raw_letters::AbstractVector{<:Integer},
                         raw_blocks::AbstractVector{<:Integer})
        length(raw_letters) == length(raw_blocks) ||
            throw(ArgumentError("letters and blocks must have equal length"))
        L = collect(Letter, raw_letters); B = collect(Int, raw_blocks)
        cl, cb = _canonical_colored(L, B)
        new(cl, cb, isempty(cl) ? "ε" : join(cl, ""))
    end
end

Base.length(w::ColoredWord) = length(w.letters)
Base.isempty(w::ColoredWord) = isempty(w.letters)
Base.:(==)(a::ColoredWord, b::ColoredWord) = a.letters == b.letters && a.blocks == b.blocks
Base.hash(w::ColoredWord, h::UInt) = hash(w.blocks, hash(w.letters, hash(:ColoredWord, h)))

"underlying circular word, forgetting the coloring"
forget_color(w::ColoredWord) = CircularWord(w.letters)
compact(w::ColoredWord) = w.str

# Canonical rotation: pick the rotation minimizing the LETTER sequence; among ties
# (rotations equal as letters) pick the lexicographically smallest TAG sequence, so
# the representative is deterministic.
function _canonical_colored(L::Vector{Letter}, B::Vector{Int})
    n = length(L)
    n <= 1 && return (copy(L), copy(B))
    bestL = copy(L); bestB = copy(B)
    for k in 1:(n-1)
        rl = vcat(@view(L[k+1:end]), @view(L[1:k]))
        rb = vcat(@view(B[k+1:end]), @view(B[1:k]))
        if _lex_less(rl, bestL) || (rl == bestL && _lex_less(rb, bestB))
            bestL = rl; bestB = rb
        end
    end
    return (bestL, bestB)
end

"""
    colored_from_expansion(ex) -> ColoredWord

Build the colored word for a full expansion `ex` (a `(rotation, blocks, result)`
named tuple from `all_full_expansions`): letters = the expanded word, each block's
letters tagged with that block's index (1,2,3,…) in tiling order.
"""
function colored_from_expansion(ex)
    letters = Letter[]; blocks = Int[]
    for (bi, (_, to)) in enumerate(ex.blocks)
        for c in to
            push!(letters, c); push!(blocks, bi)
        end
    end
    return ColoredWord(letters, blocks)
end

# ---------------------------------------------------------------------------
# Tag maps: for each NON-INCREASING rule direction, output-letter i inherits the
# tag of input-window position _TAG_MAP[name_dir][i]. Keyed by (rule name, forward).
# Only reducers and neutral moves appear (we never increase here).
# ---------------------------------------------------------------------------

# helper: identity map of length k
_idmap(k) = collect(1:k)

# Built from EXTENDED_RULES; for a rule from->to we consider BOTH directions but
# keep only the non-increasing one(s). Output length = |to| (forward) or |from|.
# For same-length (neutral) rules both directions are non-increasing.
const _TAG_MAP = let d = Dict{Tuple{Symbol,Bool},Vector{Int}}()
    for r in EXTENDED_RULES
        lf, lt = length(r.from), length(r.to)
        # forward: from -> to, non-increasing iff lt <= lf
        if lt <= lf
            # each output position maps to the SAME index in the window (prefix);
            # for ii->i, ii->ε, 131->3 etc. the survivors are a prefix of the window
            d[(r.name, true)] = _idmap(lt)
        end
        # backward: to -> from, non-increasing iff lf <= lt
        if lf <= lt
            d[(r.name, false)] = _idmap(lf)
        end
    end
    d
end

# ---------------------------------------------------------------------------
# One colored non-increasing move set
# ---------------------------------------------------------------------------

# All legal (cross-block) non-increasing colored neighbors of `w`. A move is legal
# iff its matched window's tags are NOT all identical (spans ≥2 blocks). Output
# letters inherit tags via _TAG_MAP; letters outside the window keep their tags.
function down_neighbors_colored(w::ColoredWord)
    isempty(w) && return ColoredWord[]
    L = w.letters; B = w.blocks; n = length(L)
    seen = Set{ColoredWord}()
    out = ColoredWord[]
    for r in EXTENDED_RULES
        for (from, to, forward) in ((r.from, r.to, true), (r.to, r.from, false))
            length(to) <= length(from) || continue      # non-increasing only
            patlen = length(from)
            patlen == 0 && continue                      # no ε-source here
            patlen > n && continue
            tagmap = get(_TAG_MAP, (r.name, forward), nothing)
            tagmap === nothing && continue
            for start in 1:n
                # read window letters + tags
                win = Letter[ L[mod1(start + k, n)] for k in 0:(patlen-1) ]
                win == from || continue
                wtags = Int[ B[mod1(start + k, n)] for k in 0:(patlen-1) ]
                # LEGALITY: window must span ≥2 distinct blocks
                all(==(wtags[1]), wtags) && continue
                # build result letters (as in _cyclic_replace) + tags
                newL = vcat(to, Letter[ L[mod1(start + patlen + k, n)] for k in 0:(n - patlen - 1) ])
                newtags = vcat(Int[ wtags[tagmap[j]] for j in 1:length(to) ],
                               Int[ B[mod1(start + patlen + k, n)] for k in 0:(n - patlen - 1) ])
                cw = ColoredWord(newL, newtags)
                cw == w && continue
                if !(cw in seen)
                    push!(seen, cw); push!(out, cw)
                end
            end
        end
    end
    return out
end

"""
    colored_down_path(w; maxsteps = 100_000) -> Union{Nothing, Vector{ColoredWord}}

A concrete legal (cross-block) non-increasing path of colored words from `w` down
to ε, or `nothing` if none exists. BFS over `down_neighbors_colored`.
"""
function colored_down_path(w::ColoredWord; maxsteps::Int = 100_000)
    isempty(w) && return [w]
    EMPTYC = ColoredWord(Letter[], Int[])
    prev = Dict{ColoredWord,ColoredWord}()
    seen = Set{ColoredWord}([w]); q = ColoredWord[w]; steps = 0
    while !isempty(q) && steps < maxsteps
        x = popfirst!(q); steps += 1
        isempty(x) && break
        for nb in down_neighbors_colored(x)
            if !(nb in seen)
                push!(seen, nb); prev[nb] = x; push!(q, nb)
            end
        end
    end
    target = nothing
    for s in seen
        isempty(s) && (target = s; break)
    end
    target === nothing && return nothing
    path = ColoredWord[target]
    cur = target
    while cur != w
        haskey(prev, cur) || return nothing
        cur = prev[cur]; pushfirst!(path, cur)
    end
    return path
end
