# LightLeaves.jl  —  light leaves & double leaves on MorphismGraphs.
#
# The morphism generators (identity strand, dot) and the light leaf for the PURE-UP
# subexpression (all e_k = 1) are group-free: they need no Bruhat decision. The full
# U0/U1/D0/D1 construction does — whether appending a letter goes UP or DOWN in the
# Coxeter group — and takes it from `algebra/Coxeter.jl`; see `_updown`.

# ---- morphism generators (each a MorphismGraph bottom → top) ---------------

"Identity strand of colour `s`: one leaf below, one above, joined by an edge.
bottom = [s], top = [s]."
function identity_strand(s::Integer)
    c = Int(s)
    g = WordGraph(CircularWord([c, c]), Node[], Edge[Edge(c, Leaf(1), Leaf(2))])
    MorphismGraph(g, 0, 1)                    # cut after leaf2 (=0) and after leaf1
end

"The dot/counit morphism on colour `s`: bottom = [s], top = ε. A single dot node."
function dot_morphism(s::Integer)
    c = Int(s)
    g = WordGraph(CircularWord([c]), Node[Node(:dot, [c], 0)], Edge[Edge(c, Leaf(1), NodePort(1, 1))])
    # bottom = whole boundary (leaf1), top = ε: this NEEDS cut1 == cut2 == 0 (the
    # "whole boundary" special case in _bottom_leaves/_top_leaves), NOT cut2 = 1 — with
    # a 1-leaf boundary there is only ONE gap, so cut1=0, cut2=1 both name that SAME
    # gap but the raw (non-modular) `cut1 == cut2` check in MorphismGraph.jl misses it,
    # making `top` come out as `[s]` instead of `ε`. The light_leaf induction relies
    # on dot_morphism's top being really empty.
    MorphismGraph(g, 0, 0)
end

# unit/startdot (bottom = ε, top = [s]) is NOT expressible in the two-cut model for a
# 1-leaf boundary: with one gap you cannot have BOTH a non-empty top and an empty
# bottom distinct from the reverse. It needs the tensor-based construction (a strand
# tensored with a dot), deferred with the full LL build. Do not add a broken n=1
# version — the cut arithmetic is degenerate there.

# ---- tensor of morphisms (bottoms concatenated, tops concatenated) ---------
#
# For MorphismGraphs we cannot just `tensor` the underlying graphs (diagram/Join.jl) —
# that only concatenates boundaries; the cut positions must be recomputed so the
# result's bottom = bottom(f)·bottom(g) and top = top(f)·top(g). We build it directly
# (planar-word-graph plumbing): the cuts of f and g must be re-derived from their
# (possibly non-trivial) bottom/top splits.

"""
    identity_morphism() -> MorphismGraph

The empty morphism `ε → ε`: no leaves, no nodes, no edges. The tensor unit and the
starting point of the `light_leaf` induction.
"""
identity_morphism() = MorphismGraph(WordGraph(CircularWord(Int[]), Node[], Edge[]), 0, 0)

"""
    tensor(f::MorphismGraph, g::MorphismGraph) -> MorphismGraph

Horizontal juxtaposition: `bottom(result) = bottom(f)·bottom(g)`,
`top(result) = top(f)·top(g)`. Built by hand (not via `WordGraph` `tensor`, which just
concatenates boundaries) because the two cuts have to be recomputed:

New PHYSICAL leaf order (counter-clockwise, keeps the picture planar): `bottom(f)`
leaves, `bottom(g)` leaves, then `top(g)`'s leaves (forward), then `top(f)`'s leaves
(forward). `cut1 = 0`, `cut2 = |bottom(f)| + |bottom(g)|` — reduced mod the new total
leaf count so that a fully-bottom (empty-top) result correctly hits the `cut1==cut2`
"whole boundary" special case in `_bottom_leaves`/`_top_leaves` (MorphismGraph.jl),
which compares cuts RAW, not modulo `n` (the same degeneracy fixed in `dot_morphism`
above — here it recurs whenever `top(f)` and `top(g)` are BOTH empty, since then
`|bottom(f)|+|bottom(g)|` equals the new total leaf count exactly).

NB (chaining ≥ 2 tensors): `top()` (MorphismGraph.jl) reads its
physical arc REVERSED, so getting `top(result) == top(f)·top(g)` out of a SINGLE
overall reversal of the post-cut2 arc needs the PRE-reversal physical order to be
`top(g)`'s leaves forward THEN `top(f)`'s leaves forward (own internal order kept,
piece order swapped) — reversing each piece individually first (as one might guess
from "planar ⇒ mirror both") instead reverses the piece order TWICE and silently
reproduces `top(f)·top(g)` only when both tops have length ≤ 1. The naive
double-reverse already fails on a 3-fold chained tensor
(`light_leaf([1,2,1],[1,1,1])`, i.e. plain `light_leaf_up`), giving `[2,1,1]`
instead of `[1,2,1]`.

Node indices of `g` are shifted by `length(f.graph.nodes)`; `f`'s are untouched.
Handles `f`/`g` being the empty morphism (`identity_morphism()`) correctly (no
leaves/nodes contributed from that side).
"""
function tensor(f::MorphismGraph, g::MorphismGraph)
    bf = _bottom_leaves(f); tf = _top_leaves(f)
    bg = _bottom_leaves(g); tg = _top_leaves(g)

    order = Tuple{Symbol,Int}[]
    for k in bf; push!(order, (:f, k)); end
    for k in bg; push!(order, (:g, k)); end
    for k in tg; push!(order, (:g, k)); end
    for k in tf; push!(order, (:f, k)); end

    leafmap = Dict{Tuple{Symbol,Int},Int}()
    cols = Int[]
    for (idx, (src, k)) in enumerate(order)
        leafmap[(src, k)] = idx
        push!(cols, src === :f ? leaf_colour(f.graph, k) : leaf_colour(g.graph, k))
    end
    ntot = length(order)

    offset = length(f.graph.nodes)
    nodes = vcat(f.graph.nodes, g.graph.nodes)

    trport(src::Symbol, p::Leaf) = Leaf(leafmap[(src, p.k)])
    trport(src::Symbol, p::NodePort) = src === :f ? p : NodePort(p.node + offset, p.slot)
    trport(::Symbol, p::Circle) = p

    edges = Edge[]
    for e in f.graph.edges
        push!(edges, Edge(e.colour, trport(:f, e.a), trport(:f, e.b)))
    end
    for e in g.graph.edges
        push!(edges, Edge(e.colour, trport(:g, e.a), trport(:g, e.b)))
    end

    word = CircularWord(cols)
    graph = WordGraph(word, nodes, edges)
    cut1 = 0
    cut2raw = length(bf) + length(bg)
    cut2 = ntot == 0 ? 0 : mod(cut2raw, ntot)
    return MorphismGraph(graph, cut1, cut2)
end

# ---- merge / cap generators (the D0 / D1 down-step morphisms) --------------

"""
    merge_morphism(i) -> MorphismGraph

The merge (trivalent) morphism on colour `i`: bottom = `[i,i]`, top = `[i]`. One
trivalent node, boundary word `(i i i)` (bottom leaves 1,2; top leaf 3).

WINDING ORIENTATION. The bottom arc runs DESCENDING into the slots, exactly
as in `braid_move_morphism` (see there, `slot = mod1(c − position)`):

    leaf 1 → slot 2,  leaf 2 → slot 1,  leaf 3 (top) → slot 3

Slot 3 (the upper output) stays fixed — only the two lower legs swap, so the
role "slot 1,2 = legs, slot 3 = output" is preserved. (Any renderer that reads
slot order as a rotation system would otherwise draw this node with the
opposite winding from every braid node.)
"""
function merge_morphism(i::Integer)
    c = Int(i)
    g = WordGraph(CircularWord([c, c, c]), Node[Node(:trivalent, [c], 0)],
                  Edge[Edge(c, Leaf(1), NodePort(1, 2)),
                       Edge(c, Leaf(2), NodePort(1, 1)),
                       Edge(c, Leaf(3), NodePort(1, 3))])
    MorphismGraph(g, 0, 2)          # bottom = leaves 1,2 (=[i,i]); top = leaf 3 (=[i])
end

"""
    split_morphism(i) -> MorphismGraph

The split (dual trivalent) morphism on colour `i`: bottom = `[i]`, top = `[i,i]`. Same
physical wiring as `merge_morphism` (one trivalent node, boundary word `(i i i)`, stem
= slot 3, legs = slots 2,1 descending) but with bottom/top swapped: bottom = leaf 1
(the stem), top = leaves 2,3 (the legs).
"""
function split_morphism(i::Integer)
    c = Int(i)
    g = WordGraph(CircularWord([c, c, c]), Node[Node(:trivalent, [c], 0)],
                  Edge[Edge(c, Leaf(1), NodePort(1, 3)),
                       Edge(c, Leaf(2), NodePort(1, 2)),
                       Edge(c, Leaf(3), NodePort(1, 1))])
    MorphismGraph(g, 0, 1)          # bottom = leaf 1 (=[i]); top = leaves 2,3 (=[i,i])
end

"""
    cap_morphism(i) -> MorphismGraph

The cap morphism on colour `i`: bottom = `[i,i]`, top = ε. A single edge joining the
two boundary leaves directly (no node). Bottom = the WHOLE boundary here, so — same
degeneracy as `dot_morphism` above — the cuts must be `cut1 == cut2 == 0` (raw
equality) to hit the "whole boundary, top = ε" special case, NOT `cut2 = 2` (which
names the identical gap only modulo `n = 2`, and is missed by the raw compare).
"""
function cap_morphism(i::Integer)
    c = Int(i)
    g = WordGraph(CircularWord([c, c]), Node[], Edge[Edge(c, Leaf(1), Leaf(2))])
    MorphismGraph(g, 0, 0)
end

# ---- pure-up light leaf (all e_k = 1): no group decision needed -------------

"""
    light_leaf_up(word) -> MorphismGraph

The light leaf for the all-ones subexpression of `word`: every strand is kept (U1),
so the diagram is just parallel identity strands and top == bottom == word. This is
the base case that needs no Bruhat/up-down decision; the general
`light_leaf(word, e)` (U0/U1/D0/D1) takes that decision from `_updown`.
"""
function light_leaf_up(word::Vector{Int})
    n = length(word)
    # bottom leaves 1..n, top leaves n+1..2n (reversed arc gives top == word).
    nodes = Node[]
    edges = Edge[]
    for (i, c) in enumerate(word)
        # strand i: bottom leaf i ↔ top leaf (2n - i + 1) so that reading the top arc
        # reversed reproduces `word`.
        push!(edges, Edge(c, Leaf(i), Leaf(2n - i + 1)))
    end
    wcols = vcat(word, reverse(word))          # boundary: bottom then top-arc (reversed)
    g = WordGraph(CircularWord(wcols), nodes, edges)
    return MorphismGraph(g, 0, n)              # cut after leaf 2n (=0) and after leaf n
end

# ---- braid moves as morphisms (join a braid vertex onto a map's top) --------
#
# A light-leaf top is reordered by BRAID MOVES to end in the needed
# letter (or to re-canonicalise). Combinatorially that is `braid_move_sites` /
# `braid_to_end_with` (BraidMoves.jl). Here we realise ONE braid move as an actual
# morphism `v → v'`: the strands outside the block [i,j] run straight through
# (identities); the block itself is one braid VERTEX. Post-composing this onto a map
# `f : w → v` (via `compose`) gives `w → v'` with the braid drawn.

"""
    braid_move_morphism(v, i, j) -> MorphismGraph

The morphism `v → v'` that applies the single braid move on block `v[i..j]`: outside
the block every strand is an identity; inside, the alternating block `[s,t,s,…]` is
carried by one 2m-valent braid vertex to `[t,s,t,…]` (`v' = apply_braid_move(v,i,j)`).
Bottom = `v`, top = `v'`. The block length `m = j−i+1` must be `m_{st}` (checked by
`apply_braid_move`; a genuine relation only then).
"""
function braid_move_morphism(v::AbstractVector{<:Integer}, i::Integer, j::Integer)
    vv = collect(Int, v)
    n = length(vv)
    vprime = apply_braid_move(vv, i, j)            # also validates the block
    m = j - i + 1
    s, t = vv[i], vv[i + 1 <= n ? i + 1 : i]
    # boundary: bottom leaves 1..n (=v), top leaves n+1..2n whose REVERSED reading = v'.
    # top position k (1..n) of v' sits on leaf 2n-k+1.
    topleaf(k) = 2n - k + 1
    # Slot colour follows slot PARITY (`_slot_colour`, diagram/Graph.jl). The
    # reversed wiring below (`mod1(m − k)`) shifts parity exactly when `m` is
    # EVEN — in that case the two node colours must be swapped so edge colour
    # and slot colour line up. For odd `m` parity is preserved and there is
    # nothing to do.
    ncols = iseven(m) ? [t, s] : [s, t]
    nodes = Node[Node(:braid, ncols, m)]           # the single braid vertex, node index 1
    edges = Edge[]
    # identity strands outside the block
    for p in 1:n
        (i <= p <= j) && continue
        push!(edges, Edge(vv[p], Leaf(p), Leaf(topleaf(p))))
    end
    # block: braid slots 1..2m. Convention: slots are the cyclic order.
    #
    # WINDING ORIENTATION. The boundary
    # is read counter-clockwise at the BOTTOM and clockwise at the TOP (`top()`
    # reverses the top arc, MorphismGraph.jl). Both halves therefore run
    # DESCENDING into the slots, not ascending: with `c1 = m` for the bottom and
    # `c2 = mod1(c1 + m, 2m)` for the top (the top sits directly opposite the
    # bottom, `opposite_slot`), `slot = mod1(c − position)`. Result:
    #
    #     m = 2:  leaves 1..4 → slots [2, 1, 4, 3]
    #     m = 3:  leaves 1..6 → slots [3, 2, 1, 6, 5, 4]
    #
    # The reversal hits slots of EQUAL parity each time, and since a braid
    # slot's colour depends only on parity (`_slot_colour`, diagram/Graph.jl),
    # the colour assignment stays consistent — the naive "just mirror the top"
    # variant `[1,2,4,3]` would not, and is therefore colour-impossible for
    # even `m`.
    c1 = m
    c2 = mod1(c1 + m, 2m)
    for k in 0:(m - 1)
        push!(edges, Edge(vv[i + k], Leaf(i + k), NodePort(1, mod1(c1 - k, 2m))))
    end
    for k in 0:(m - 1)
        # `topleaf` runs BACKWARD over the boundary (`topleaf(i+k)` decreases as
        # k grows), so the position here is counted from the end: `m−1−k`.
        push!(edges, Edge(vprime[i + k], Leaf(topleaf(i + k)), NodePort(1, mod1(c2 - (m - 1 - k), 2m))))
    end
    wcols = vcat(vv, reverse(vprime))              # boundary colours: bottom then top-arc reversed
    g = WordGraph(CircularWord(wcols), nodes, edges)
    return MorphismGraph(g, 0, n)                  # bottom = v (leaves 1..n), top = v'
end

"""
    braid_top(f, i, j) -> MorphismGraph

Post-compose the single braid move on block `[i,j]` of `f`'s TOP word onto `f`: if
`f : w → v` then the result is `w → v'` with `v' = apply_braid_move(top(f), i, j)`,
the braid vertex joined onto v. Just `compose(f, braid_move_morphism(top(f), i, j))`.
"""
braid_top(f::MorphismGraph, i::Integer, j::Integer) =
    compose(f, braid_move_morphism(top(f), i, j))

"""
    braid_top_to(f, vprime) -> MorphismGraph

Post-compose a whole SEQUENCE of braid moves onto `f`'s top, reordering `top(f)` into
the reduced word `vprime` of the same element (via `braid_to_word`). Result: `w → vprime`.
"""
function braid_top_to(f::MorphismGraph, vprime::AbstractVector{<:Integer})
    cur = f
    moves = braid_to_word(top(f), vprime)
    for (i, j) in moves
        cur = braid_top(cur, i, j)
    end
    return cur
end

# ---- the up/down decision --------------------------------------------------
#
# `_updown(t, s)` answers: does appending letter `s` to the element with reduced word
# `t` go UP (length +1) or DOWN in the A₃ Coxeter group? Appending `s` goes DOWN iff
# `s` is a RIGHT DESCENT of the element. We build the A₃ group once and cache it.

const _A3 = Ref{Any}(nothing)
function _a3()
    _A3[] === nothing && (_A3[] = coxeter_group_a3())
    return _A3[]                                   # (W, gens)
end

"The A₃ group element with reduced word `t` (letters in 1:3)."
function _element(t::Vector{Int})
    W, gens = _a3()
    el = one(W)
    for i in t; el = el * gens[i]; end
    return el
end

"`:up` if appending letter `s` to the element of reduced word `t` increases length,
else `:down` (s is a right descent)."
function _updown(t::Vector{Int}, s::Int)
    is_right_descent(_element(t), s) ? :down : :up
end

# ---- subexpression combinatorics ------------------------------------------
# A₃-specialised) ----------------------------------------------------------

"""
    subexpressions(n) -> Vector{Vector{Int}}

All `2^n` 01-subexpressions `e ∈ {0,1}ⁿ` of a length-`n` word, as `Vector{Int}`s of
0s and 1s (`e[k]` decides whether letter `k` is KEPT (`1`) or DROPPED (`0`)).
`subexpressions(0) == [Int[]]` (the unique subexpression of the empty word).
"""
function subexpressions(n::Integer)
    n == 0 && return [Int[]]
    n >= 0 || error("subexpressions: n must be ≥ 0")
    return [[((mask >> (k - 1)) & 1) for k in 1:n] for mask in 0:(2^n - 1)]
end

"""
    decorations(word, e) -> Vector{Symbol}

The Bruhat-stroll decoration of each step of `(word, e)`: walk the CANONICAL
short-lex word `cur` of the element expressed up to that point (start `Int[]`); at step `k`
with letter `s = word[k]`, `up = _updown(cur, s) == :up`:

    e[k] == 1 && up    -> :U1   (keep the strand, element grows)
    e[k] == 0 && up    -> :U0   (drop the strand, element unchanged — it dies)
    e[k] == 0 && !up   -> :D0   (element unchanged — MERGE two `s`-strands)
    e[k] == 1 && !up   -> :D1   (element shrinks — CAP both `s`-strands)

After a `U1`/`D1` step (`e[k] == 1`) `cur` is updated to `canonical_word(cur · s)`
(this single call correctly handles BOTH the up-extension and the down-cancellation,
since `canonical_word` just re-reads off the short-lex word of the resulting group
element, which is where the cancelling happens). `cur` is left unchanged on a `0`-step
(the element does not move). D0 = merge, D1 = cap.
"""
function decorations(word::Vector{Int}, e::Vector{Int})
    n = length(word)
    length(e) == n || error("decorations: e (length $(length(e))) and word (length $n) size mismatch")
    cur = Int[]
    decs = Vector{Symbol}(undef, n)
    for k in 1:n
        s = word[k]
        up = _updown(cur, s) == :up
        if e[k] == 1
            decs[k] = up ? :U1 : :D1
            cur = canonical_word(vcat(cur, [s]))
        else
            decs[k] = up ? :U0 : :D0
        end
    end
    return decs
end

"""
    defect(word, e) -> Int

`#{:U0} − #{:D0}` over `decorations(word, e)` — the grading/defect of the light leaf
`light_leaf(word, e)`.
"""
function defect(word::Vector{Int}, e::Vector{Int})
    decs = decorations(word, e)
    count(==(:U0), decs) - count(==(:D0), decs)
end

"""
    expressed_word(word, e) -> Vector{Int}

The canonical short-lex reduced word of the element expressed by `(word, e)`
(walking the same `cur` update as `decorations`, keeping only the final value).
"""
function expressed_word(word::Vector{Int}, e::Vector{Int})
    n = length(word)
    length(e) == n || error("expressed_word: e (length $(length(e))) and word (length $n) size mismatch")
    cur = Int[]
    for k in 1:n
        e[k] == 1 && (cur = canonical_word(vcat(cur, [word[k]])))
    end
    return cur
end

# ---- the full light leaf (U0/U1/D0/D1 induction) ---------------------------

"""
    light_leaf(word, e) -> MorphismGraph

The light leaf `LL(word, e)`: bottom = `word`, top = `expressed_word(word, e)` (the
canonical short-lex word of the expressed element). Built letter by letter,
tensoring the next strand onto `cur` and, on a down step, first re-canonicalising
`top(cur)` (BEFORE tensoring the new strand) to end in the needed letter via
`braid_to_end_with`/`braid_top_to`, so the merge/cap always acts on the last two
strands:

    U1: cur = tensor(cur, identity_strand(i))                    # strand kept
    U0: cur = tensor(cur, dot_morphism(i))                       # strand dies (dot)
    D0: braid top(cur) to end in i, THEN tensor identity_strand(i)
        (now top ends `…, i, i`), then compose with a MERGE on the last two strands
    D1: same braid + tensor, then compose with a CAP on the last two strands

Asserts `bottom(result) == word` and `top(result) == expressed_word(word, e)` as a
consistency check.
"""
function light_leaf(word::Vector{Int}, e::Vector{Int})
    n = length(word)
    length(e) == n || error("light_leaf: e (length $(length(e))) and word (length $n) size mismatch")
    decs = decorations(word, e)
    cur = identity_morphism()
    for k in 1:n
        i = word[k]
        dec = decs[k]
        if dec === :U1
            cur = tensor(cur, identity_strand(i))
        elseif dec === :U0
            cur = tensor(cur, dot_morphism(i))
        elseif dec === :D0 || dec === :D1
            told = top(cur)
            newtop, _ = braid_to_end_with(told, i)
            cur = braid_top_to(cur, newtop)
            cur = tensor(cur, identity_strand(i))       # top now = newtop · i (ends i,i)
            prefix = newtop[1:(end - 1)]
            gen = dec === :D0 ? merge_morphism(i) : cap_morphism(i)
            cur = compose(cur, tensor(light_leaf_up(prefix), gen))
        else
            error("light_leaf: unknown decoration $dec")
        end
        # Re-canonicalise the top to short-lex (the top word is ALWAYS the
        # canonical short-lex reduced word, so double leaves compose without ambiguity).
        # A U1 append (or a merge/cap) can leave the top a NON-short-lex reduced word of
        # the same element — e.g. appending 1 after …3 gives …3,1 where short-lex wants
        # …1,3 (the 1↔3 commutation). Braid it back to canonical form; braid_top_to is a
        # no-op when the top is already short-lex.
        canon = canonical_word(top(cur))
        top(cur) == canon || (cur = braid_top_to(cur, canon))
    end
    @assert bottom(cur) == word "light_leaf($word, $e): bottom mismatch, got $(bottom(cur))"
    ew = expressed_word(word, e)
    @assert top(cur) == ew "light_leaf($word, $e): top mismatch, got $(top(cur)), want $ew"
    return cur
end

# ---- double leaves + enumeration -------------------------------------------

"""
    double_leaf(xword, e, yword, f) -> MorphismGraph

`compose(light_leaf(xword, e), flip(light_leaf(yword, f)))`: a morphism `xword →
yword`, factored through the common expressed element `z = expressed_word(xword, e)
== expressed_word(yword, f)` (checked; errors if they differ).
"""
function double_leaf(xword::Vector{Int}, e::Vector{Int}, yword::Vector{Int}, f::Vector{Int})
    zx = expressed_word(xword, e)
    zy = expressed_word(yword, f)
    zx == zy || error("double_leaf: expressed elements differ: $zx (from x,e) vs $zy (from y,f)")
    return compose(light_leaf(xword, e), flip(light_leaf(yword, f)))
end

"""
    light_leaves(xword, zword; degree = nothing) -> Vector{Tuple{Vector{Int},MorphismGraph}}

All `(e, light_leaf(xword, e))` with `expressed_word(xword, e) == zword`, optionally
filtered to a fixed `defect(xword, e) == degree`.
"""
function light_leaves(xword::Vector{Int}, zword::Vector{Int}; degree::Union{Nothing,Int} = nothing)
    n = length(xword)
    out = Tuple{Vector{Int},MorphismGraph}[]
    for e in subexpressions(n)
        expressed_word(xword, e) == zword || continue
        (degree === nothing || defect(xword, e) == degree) || continue
        push!(out, (e, light_leaf(xword, e)))
    end
    return out
end

"""
    double_leaves(xword, yword; degree = nothing)
        -> Vector{@NamedTuple{e::Vector{Int}, f::Vector{Int}, z::Vector{Int}, degree::Int, morphism::MorphismGraph}}

All double leaves `xword → yword`: for every element `z` reachable as BOTH
`expressed_word(xword, ·)` and `expressed_word(yword, ·)` (i.e. `z ≤ x` and `z ≤ y`
in Bruhat order, read off combinatorially), every pair `(e, f)` from
`light_leaves(xword, z) × light_leaves(yword, z)` gives one entry
`double_leaf(xword, e, yword, f)`, with `degree = defect(xword,e) + defect(yword,f)`
(optionally filtered).
"""
function double_leaves(xword::Vector{Int}, yword::Vector{Int}; degree::Union{Nothing,Int} = nothing)
    nx = length(xword); ny = length(yword)
    zsx = Set(expressed_word(xword, e) for e in subexpressions(nx))
    zsy = Set(expressed_word(yword, f) for f in subexpressions(ny))
    zs = intersect(zsx, zsy)
    T = @NamedTuple{e::Vector{Int}, f::Vector{Int}, z::Vector{Int}, degree::Int, morphism::MorphismGraph}
    out = T[]
    for z in zs
        lx = light_leaves(xword, z)
        ly = light_leaves(yword, z)
        for (e, _) in lx, (f, _) in ly
            d = defect(xword, e) + defect(yword, f)
            (degree === nothing || d == degree) || continue
            dl = double_leaf(xword, e, yword, f)
            push!(out, (e = e, f = f, z = z, degree = d, morphism = dl))
        end
    end
    return out
end

# ---- A₃ enumeration helper -------------------------------------------------

"""
    a3_words() -> Vector{Vector{Int}}

The canonical short-lex reduced word of every element of A₃ (24 words, from
`enumerate_whole_group`); picks word triples for double-leaf compositions.
"""
function a3_words()
    W, _ = _a3()
    return [Int.(short_lex(w)) for w in enumerate_whole_group(W)]
end
