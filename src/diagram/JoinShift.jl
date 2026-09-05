# diagram/JoinShift.jl — shift_right/shift_left: reinterpreting a MorphismGraph's
# cut2 to move a boundary letter between bottom and top, without touching the graph.
# Included after Join.jl.
#
# ---- shift_right / shift_left: reinterpret the SAME graph with a moved cut -----
#
# Both functions move ONLY `cut2` (never `cut1`, never a single node/edge of the
# graph) — the letter that "moves" between bottom and top is never touched at the
# data level, only which arc it is read as part of. Derivation (see the
# shift_right/shift_left docstrings): `bottom(m)`
# reads leaves `cut1+1 … cut2` FORWARD (MorphismGraph.jl `_bottom_leaves`), so its
# LAST letter is always `Leaf(cut2)`. `top(m)` reads leaves `cut2+1 … cut1` forward
# and then REVERSES the result (`_top_leaves` + the `reverse(...)` in `top`), so its
# LAST letter is always `Leaf(cut2+1)` (the FIRST leaf of the forward top-arc walk,
# since reversing puts it last). Moving `Leaf(cut2)` from bottom's end to top's end
# is therefore exactly `cut2 -= 1` (bottom loses its last leaf, and that same leaf
# becomes the new `cut2+1`, i.e. top's forward-walk start, i.e. top's new last
# letter after reversal). `shift_left` undoes this: `cut2 += 1`.

"""
    _empty_bottom_state(m, c1, newc2) -> MorphismGraph

Build the `MorphismGraph` for "bottom has just become empty" at effective cut1
`c1` (already normalised: `_EMPTY_BOTTOM_CUT ↦ 0`) and effective new cut2
`newc2 == c1`. The plain pair `(c1, c1)` is read by `_bottom_leaves`/`_top_leaves`
as "bottom = WHOLE boundary, top = ε" for **every** value of `c1` (not just `0`) —
so it is the wrong reading here; only the dedicated `_EMPTY_BOTTOM_CUT` sentinel
gives "bottom = ε, top = whole". That sentinel is however wired (in
`MorphismGraph.jl`) to always read the whole boundary starting at `Leaf(1)`,
regardless of `cut1`'s previous value — so it is only a faithful re-reading of THIS
graph when `c1 == 0` (mod `n`). For `c1 != 0` there is no `(cut1,cut2)` pair in the
current `MorphismGraph` representation that expresses "bottom = ε starting at gap
`c1`" (every non-sentinel `cut1==cut2` reads as the OTHER degenerate case, for any
`cut1` value) — a structural gap in the two-cut model, so this function raises a
clear error instead of returning a silently-wrong morphism.
"""
function _empty_bottom_state(m::MorphismGraph, c1::Int, newc2::Int)
    n = length(m.graph.word)
    if mod(c1, n) == 0
        return MorphismGraph(m.graph, _EMPTY_BOTTOM_CUT, _EMPTY_BOTTOM_CUT)
    end
    error("shift_right: bottom just became empty at a non-zero cut1=$(c1) — " *
          "the two-cut MorphismGraph representation has no sentinel for \"bottom=ε\" " *
          "except when cut1≡0 (mod n); this shift cannot be represented without " *
          "changing cut1 (which shift_right never does). Not supported.")
end

"""
    shift_right(m::MorphismGraph, k::Integer = 1) -> MorphismGraph

Reinterpret the boundary cut: move the LAST `k` letters of `bottom(m)` to the END
of `top(m)`, one leaf at a time, WITHOUT touching the graph (same nodes, same
edges, same boundary word) — only `cut2` moves, `cut1` never changes. If
`m : x → y` and `x = x'·s` (`s` the last letter of `x`), `shift_right(m)` is the
morphism `x' → y·s`; `k`-fold is `k` single-letter shifts (`shift_right(m,2) ==
shift_right(shift_right(m))`).

**Derivation** (see the module-level note above `_empty_bottom_state`): `bottom`
reads its arc forward, so its last letter is `Leaf(cut2)`; `top` reads its arc
forward then reverses, so ITS last letter is `Leaf(cut2+1)` — the same leaf that
`cut2 -= 1` hands to it. Hence `shift_right(m, k) = MorphismGraph(m.graph, cut1,
cut2 - k)` (mod `n`), modulo the degenerate bookkeeping below.

**Degenerate / edge cases:**
- `k > length(bottom(m))`: no such letter exists to move — raises `ArgumentError`.
  A full trip around is `shift_right(m, length(bottom(m)))`, landing on the
  empty-bottom sentinel below, not a return to `m`.
- `n == 0` (empty boundary): `bottom(m) == top(m) == []`, so `length(bottom(m)) ==
  0 < k` for any `k ≥ 1` — caught by the same `ArgumentError` above, no special
  case needed.
- The shift makes bottom become empty for the first time (`k == length(bottom(m))`
  and bottom was nonempty): the plain pair `(cut1, cut1)` would be misread as
  "bottom=whole/top=ε" (the OTHER degenerate case), so the result switches to the
  `_EMPTY_BOTTOM_CUT` sentinel instead — but ONLY valid when `cut1 ≡ 0 (mod n)`
  (see `_empty_bottom_state`); this is the only case ever produced by the
  library's own constructors (`identity_morphism`, `dot_morphism`, `cap_morphism`,
  `merge_morphism`, `identity_strand`, `tensor`, `compose`, `double_leaf` all fix
  `cut1 = 0`; only `flip` can produce a nonzero `cut1`, in which case this function
  errors rather than silently mis-shifting).
- Starting from an ALREADY-empty bottom (`m.cut1 == m.cut2` as the sentinel, or
  any `cut1==cut2` "whole/ε" reading with `length(bottom(m)) == 0`): `k ≥ 1 >
  length(bottom(m)) == 0`, caught by the `ArgumentError` above.

The graph is untouched
(`canonical_key(shift_right(m).graph) == canonical_key(m.graph)`), and
`shift_left ∘ shift_right == id` (and vice versa) wherever both sides are defined.
"""
function shift_right(m::MorphismGraph, k::Integer = 1)
    k >= 1 || error("shift_right: k must be ≥ 1, got $k")
    n = length(m.graph.word)
    blen = length(bottom(m))
    k <= blen || throw(ArgumentError(
        "shift_right: k=$k exceeds length(bottom(m))=$blen — no letter to shift"))
    c1 = m.cut1 == _EMPTY_BOTTOM_CUT ? 0 : m.cut1
    c2 = m.cut2 == _EMPTY_BOTTOM_CUT ? 0 : m.cut2
    newc2 = mod(c2 - k, n)
    if newc2 == mod(c1, n) && k == blen
        return _empty_bottom_state(m, c1, newc2)
    end
    return MorphismGraph(m.graph, c1, newc2)
end

"""
    shift_left(m::MorphismGraph, k::Integer = 1) -> MorphismGraph

Reinterpret the boundary cut: move the LAST `k` letters of `top(m)` to the END of
`bottom(m)`, one leaf at a time, WITHOUT touching the graph — only `cut2` moves,
`cut1` never changes. If `m : x → y` and `y = y'·t` (`t` the last letter of `y`),
`shift_left(m)` is the morphism `x·t → y'`; `k`-fold is `k` single-letter shifts.
Exact inverse of `shift_right` wherever both are defined:
`shift_left(shift_right(m)) == m == shift_right(shift_left(m))`.

**Derivation:** `top`'s last letter is `Leaf(cut2+1)` (see `shift_right`'s note);
moving it onto the end of `bottom` (whose forward walk ends at `Leaf(cut2)`) is
exactly `cut2 += 1`. Hence `shift_left(m, k) = MorphismGraph(m.graph, cut1, cut2 +
k)` (mod `n`), modulo the degenerate bookkeeping below.

**Degenerate / edge cases:**
- `k > length(top(m))`: raises `ArgumentError`, matching `shift_right` (no cyclic
  wrap-around).
- `n == 0`: `top(m) == []`, caught by the same `ArgumentError`.
- The shift makes TOP become empty for the first time (`k == length(top(m))` and
  top was nonempty): unlike `shift_right`'s empty-BOTTOM case, this does **not**
  need the `_EMPTY_BOTTOM_CUT` sentinel — the plain pair `(cut1, cut1)` (mod `n`)
  is READ AS "bottom = whole, top = ε", which is exactly the state wanted here (no
  ambiguity, since this is the case the plain pair already means).
- Starting from an already-empty top: `k ≥ 1 > length(top(m)) == 0`, caught by the
  `ArgumentError` above.
- Leaving the `_EMPTY_BOTTOM_CUT` sentinel (bottom was ε, `cut1==cut2==
  _EMPTY_BOTTOM_CUT`): the effective `c1` is normalised to `0` and `c2` to `0`, so
  `shift_left` produces `(0, k mod n)` — a plain (non-sentinel) cut pair, which is
  correct since bottom is no longer empty.

The graph is untouched by `shift_left`, matching `shift_right`.
"""
function shift_left(m::MorphismGraph, k::Integer = 1)
    k >= 1 || error("shift_left: k must be ≥ 1, got $k")
    n = length(m.graph.word)
    tlen = length(top(m))
    k <= tlen || throw(ArgumentError(
        "shift_left: k=$k exceeds length(top(m))=$tlen — no letter to shift"))
    c1 = m.cut1 == _EMPTY_BOTTOM_CUT ? 0 : m.cut1
    c2 = m.cut2 == _EMPTY_BOTTOM_CUT ? 0 : m.cut2
    newc2 = mod(c2 + k, n)
    return MorphismGraph(m.graph, c1, newc2)
end

