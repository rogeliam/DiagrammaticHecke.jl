# BraidMoves.jl  —  word combinatorics: braid moves that reorder a reduced word.
#
# A "braid move" on a reduced word replaces a length-m alternating block
# [s,t,s,…] by [t,s,t,…] (the two sides of the braid relation for m = m_{st}).
# Applying a sequence of braid moves turns one reduced word for an element w into
# another. Light leaves need this: before a "down" step on a letter s, the current
# TOP word must be reordered to END in s (so merge/cap can act on the last strand).
#
# These functions are PURELY combinatorial (words = Vector{Int} over A₃). The
# morphism-level realisation (inserting the actual braid VERTEX at a move) is a
# separate step.
#
# The braid-move ORCHESTRATION (`braid_to_end_with`) is the standard descent
# recursion: peel off a right descent, braid the shorter word, put the letter
# back. The group-theoretic primitives it leans on — element equality, `is_right_descent`, `coxeter_matrix`
# and `short_lex` (the canonical reduced word) — all come from `algebra/Coxeter.jl`,
# so "are two words the same element" is not reimplemented here.

# ---- A₃ element helpers (local; algebra/Coxeter.jl does the real work) ------

const _BM_A3 = Ref{Any}(nothing)
function _bm_a3()
    _BM_A3[] === nothing && (_BM_A3[] = coxeter_group_a3())
    return _BM_A3[]                                    # (W, gens)
end

"The A₃ group element of the word `word` (right-multiply from the identity)."
function _bm_element(word)
    W, _ = _bm_a3()
    w = one(W)
    for s in word
        w = right_multiply(w, Int(s))
    end
    return w
end

"""
    same_element(u, v) -> Bool

Do the two words `u`, `v` denote the SAME element of A₃? (Element equality in
`algebra/Coxeter.jl`, i.e. equality of permutations — no word comparison.)
"""
same_element(u, v) = _bm_element(u) == _bm_element(v)

"""
    is_reduced(word) -> Bool

Is `word` a REDUCED expression (its length equals the Coxeter length of its
element)? Uses `short_lex` (the canonical reduced word) to read off the length.
"""
is_reduced(word) = length(word) == length(short_lex(_bm_element(word)))

"""
    canonical_word(word) -> Vector{Int}

The canonical (short-lex) reduced word of `word`'s element (`algebra/Coxeter.jl`).
Two words share this iff they are the same element.
"""
canonical_word(word) = Int.(short_lex(_bm_element(word)))

# ---- braid moves ------------------------------------------------------------

"""
    apply_braid_move(word, i, j) -> Vector{Int}

Apply the braid move on the block `word[i..j]`: the block must be alternating
`[s,t,s,…]` of length `m = j−i+1`, and is replaced by `[t,s,t,…]` (same length).
Does NOT check that `m == m_{st}` — the caller guarantees that (it is a genuine
relation only when the block length equals `m_{st}`).
"""
function apply_braid_move(word::AbstractVector{<:Integer}, i::Integer, j::Integer)
    @assert 1 ≤ i ≤ j ≤ length(word) "braid-move block [$i,$j] out of range 1:$(length(word))"
    s = word[i]
    t = i < j ? word[i + 1] : word[i]
    @assert s != t || i == j "block [$i,$j] is not alternating (starts $s,$s)"
    new = collect(Int, word)
    for k in i:j
        new[k] = iseven(k - i) ? t : s          # [s,t,s,…] ↦ [t,s,t,…]
    end
    return new
end

"""
    replay_braid_moves(word, moves) -> Vector{Int}

Apply braid `moves` (a list of `(i,j)` blocks, as returned by `braid_to_end_with`)
to `word` in order, returning the resulting word.
"""
function replay_braid_moves(word::AbstractVector{<:Integer}, moves)
    w = collect(Int, word)
    for (i, j) in moves
        w = apply_braid_move(w, i, j)
    end
    return w
end

"""
    braid_to_end_with(word, t) -> (Vector{Int}, Vector{Tuple{Int,Int}})

Given a REDUCED `word` for an A₃ element `w` and a right descent `t` of `w`, return

  * `newword` : another reduced word for the SAME `w` that ENDS in `t`, and
  * `moves`   : braid moves `(i,j)` that, applied in order to `word`, produce it.

Errors if `t` is not a right descent of `w`, or if `word` is not reduced.
"""
function braid_to_end_with(word::AbstractVector{<:Integer}, t::Integer)
    isempty(word) && error("cannot braid an empty word")
    W, _ = _bm_a3()
    w = _bm_element(word)
    is_reduced(word) || error("word $(collect(Int, word)) is not reduced")
    is_right_descent(w, Int(t)) ||
        error("letter $t is not a right descent of $(collect(Int, word))")

    word = collect(Int, word)
    s = word[end]
    s == t && return (word, Tuple{Int,Int}[])

    M = coxeter_matrix(W)
    mst = Int(M[s, t])

    braided = word
    moves = Tuple{Int,Int}[]
    # Bring the prefix into a form ending in the alternating tail t,s,… recursively.
    for i in 2:mst
        endletter = iseven(i) ? Int(t) : Int(s)
        prefix = braided[1:(length(braided) - i + 1)]
        newprefix, bmoves = braid_to_end_with(prefix, endletter)
        braided = vcat(newprefix, braided[(length(braided) - i + 2):end])
        append!(moves, bmoves)
    end

    # Now the word ends in the length-mst alternating block; braid it to swap ends.
    oldend = reverse([isodd(k) ? s : t for k in 1:mst])
    newend = reverse([isodd(k) ? t : s for k in 1:mst])
    @assert braided[(end - mst + 1):end] == oldend "braid_to_end_with: unexpected tail $(braided[(end-mst+1):end]) (want $oldend)"
    braided = vcat(braided[1:(end - mst)], newend)
    push!(moves, (length(braided) - mst + 1, length(braided)))

    return braided, moves
end

"""
    braid_to_word(v, vprime) -> Vector{Tuple{Int,Int}}

Find a sequence of braid moves that turns the reduced word `v` into the reduced
word `vprime`, when both are reduced words for the SAME A₃ element. Returns the
move list (replay it on `v` to get `vprime`); errors if `v`, `vprime` are not the
same reduced element.

Thin wrapper over `braid_to_end_with`: peel matching letters off the RIGHT end,
and whenever `v` and `vprime` disagree at the current last position, braid `v` to
end with the letter `vprime` needs there, then recurse on the shortened words.
"""
function braid_to_word(v::AbstractVector{<:Integer}, vprime::AbstractVector{<:Integer})
    is_reduced(v) || error("v = $(collect(Int, v)) is not reduced")
    is_reduced(vprime) || error("vprime = $(collect(Int, vprime)) is not reduced")
    same_element(v, vprime) ||
        error("$(collect(Int, v)) and $(collect(Int, vprime)) are different elements")

    cur = collect(Int, v)
    target = collect(Int, vprime)
    moves = Tuple{Int,Int}[]
    # Fix the word from the right end inward. At each step the suffixes already
    # match; make position p agree by braiding `cur[1:p]` to end with target[p].
    for p in length(cur):-1:1
        if cur[p] != target[p]
            prefix = cur[1:p]
            newprefix, bmoves = braid_to_end_with(prefix, target[p])
            # shift move positions? no — braid_to_end_with works on prefix 1:p, and
            # those indices coincide with cur's indices (prefix is a leading slice).
            append!(moves, bmoves)
            cur = vcat(newprefix, cur[(p + 1):end])
        end
    end
    @assert cur == target "braid_to_word failed: got $cur, want $target"
    return moves
end

# ---- the reduced-expression (reex) graph ------------------------------------

"""
    braid_move_sites(word) -> Vector{Tuple{Int,Int}}

All positions `(i, j)` where a single braid move applies to `word`: a maximal-ish
alternating block `[s,t,s,…]` of length EXACTLY `m_{st}` (so the move is a genuine
relation). In A₃ that means length-3 blocks `sts` with `m_{st}=3` (the pairs
{1,2},{2,3}) and length-2 blocks for the commuting pair {1,3} (`m_{13}=2`, i.e. the
commutation `13↔31`). Every returned move keeps the word reduced and the element
fixed.
"""
function braid_move_sites(word::AbstractVector{<:Integer})
    W, _ = _bm_a3()
    M = coxeter_matrix(W)
    w = collect(Int, word)
    n = length(w)
    sites = Tuple{Int,Int}[]
    for i in 1:(n - 1)
        s = w[i]
        t = w[i + 1]
        s == t && continue                        # block must alternate two colours
        mst = Int(M[s, t])
        j = i + mst - 1                            # the only block length that is a relation
        j > n && continue
        # is w[i:j] the alternating block s,t,s,t,… of length mst?
        if all(w[i + k] == (iseven(k) ? s : t) for k in 0:(mst - 1))
            push!(sites, (i, j))
        end
    end
    return sites
end

"""
    reduced_words(word) -> Vector{Vector{Int}}

Every reduced word for the element of `word`, reached by braid moves. Breadth-first
closure under `braid_move_sites` (Matsumoto: all reduced words of one element are
connected by braid moves). Returns them in discovery order (starting from `word`).
"""
function reduced_words(word::AbstractVector{<:Integer})
    is_reduced(word) || error("word $(collect(Int, word)) is not reduced")
    start = collect(Int, word)
    seen = Set([start])
    order = [start]
    queue = [start]
    while !isempty(queue)
        cur = popfirst!(queue)
        for (i, j) in braid_move_sites(cur)
            nxt = apply_braid_move(cur, i, j)
            if !(nxt in seen)
                push!(seen, nxt); push!(order, nxt); push!(queue, nxt)
            end
        end
    end
    return order
end

"""
    reex_graph(word) -> (words, edges)

The reduced-expression graph of the element of `word`: `words` is the list of all
its reduced words (`reduced_words`), and `edges` is the list of single-braid-move
connections `(a, b, (i,j))` where `words[a]` becomes `words[b]` by the braid move
at block `(i,j)`. Undirected in effect (each move is its own inverse), but each
unordered pair is listed once with `a < b`. For w₀ = 121321 this is the graph of all
reduced words joined step by step.
"""
function reex_graph(word::AbstractVector{<:Integer})
    words = reduced_words(word)
    index = Dict(w => k for (k, w) in enumerate(words))
    edges = Tuple{Int,Int,Tuple{Int,Int}}[]
    for (a, w) in enumerate(words)
        for (i, j) in braid_move_sites(w)
            b = index[apply_braid_move(w, i, j)]
            a < b && push!(edges, (a, b, (i, j)))
        end
    end
    return words, edges
end
