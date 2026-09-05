# Zamolodchikov.jl  —  the A₃ Zamolodchikov relation as a 14-step braid cycle.
#
# The longest element w₀ of A₃ has 16 reduced words (14 "essentially", because the
# 13-commutation gives 132132 two ways). The Zamolodchikov relation lives on the
# HEXAGON connecting the two reduced words 321323 and 121321 of w₀: there are two
# braid-move paths between them, and the relation says the two composite morphisms are
# EQUAL.
#
# The two sides are stored as paths of SINGLE braid moves. Each walks
# 321323 → 121321 in 7 moves, every move a genuine braid block with no mismatch:
#
#   LHS path (7 moves):  321323 → 323123 → 232123 → 231213 → 213213 → 213231
#                        → 212321 → 121321
#   RHS path (7 moves):  321323 → 321232 → 312132 → 132132 → 132312 → 123212
#                        → 123121 → 121321
#
# ORIENTATION: the walk is built BOTTOM→TOP starting from the CANONICAL
# short-lex word 121321 at the bottom — like a light-leaf target sits at the bottom and
# strands build upward. So the 14-step CYCLE starts and ends at 121321:
#   steps 1..7  = LHS path read BACKWARD  (121321 → … → 321323),  midpoint = 321323;
#   steps 8..14 = RHS path FORWARD        (321323 → … → 121321),  back to start.
# `Zamo(1,8)` then has bottom = 121321, top = 321323. The relation "LHS = RHS" is
# "steps 1..7 = steps 8..14 reversed" (both are morphisms 121321 → 321323 once you read
# the RHS half backward).
#
# INDEXING: the cycle is stored as 14 DISTINCT words (no repeated
# 121321). Word indices are CYCLIC (word 14 → word 1 is a step). `Zamo(i,j)` walks
# forward i→j (i>j wraps), and `Zamo(i,i)` is the WHOLE 14-step cycle, not the identity.
#
# ⚠ ONE VERTEX PER STEP. The relation has a term that performs TWO DISJOINT braids
# at once. Here it is two consecutive steps instead, so that every step is exactly
# one braid vertex — the shape `braid_top` builds. The two orderings describe the
# same morphism, but anything that assumes the two braids sit at the SAME HEIGHT
# (a relation stated on the pair, a renderer aligning them) sees a difference.

# Each step is a braid move (block start i, block end j) applied to the current word,
# taking word → apply_braid_move(word, i, j). Stored as the sequence of words so a step
# is fully self-describing; the move is recovered by comparing consecutive words.

const _ZAMO_LHS_PATH = [
    [3,2,1,3,2,3],  # start
    [3,2,3,1,2,3],  # 13 @ 3:4
    [2,3,2,1,2,3],  # 32 @ 1:3
    [2,3,1,2,1,3],  # 21 @ 3:5
    [2,1,3,2,1,3],  # 31 @ 2:3
    [2,1,3,2,3,1],  # 13 @ 5:6
    [2,1,2,3,2,1],  # 32 @ 3:5
    [1,2,1,3,2,1],  # 21 @ 1:3  (= 121321 = target)
]

const _ZAMO_RHS_PATH = [
    [3,2,1,3,2,3],  # start
    [3,2,1,2,3,2],  # 32 @ 4:6
    [3,1,2,1,3,2],  # 21 @ 2:4
    [1,3,2,1,3,2],  # 31 @ 1:2
    [1,3,2,3,1,2],  # 13 @ 4:5
    [1,2,3,2,1,2],  # 32 @ 2:4
    [1,2,3,1,2,1],  # 21 @ 4:6
    [1,2,1,3,2,1],  # 31 @ 3:4  (= 121321 = target)
]

"""
    zamo_words() -> Vector{Vector{Int}}

The 14 DISTINCT words of the Zamolodchikov cycle, as a CYCLE (word 14 → word 1 is also
one braid move). Oriented to start at the canonical short-lex word `121321`
(bottom-to-top): words 1..8 are the LHS path read BACKWARD
`121321 → 321323`, words 8..14 continue the RHS path FORWARD and then close back to
word 1. So `zamo_words()[1] == [1,2,1,3,2,1]`, `zamo_words()[8] == [3,2,1,3,2,3]`, and
every consecutive pair (including 14→1, cyclically) differs by ONE braid move.
"""
function zamo_words()
    lhs_back = reverse(_ZAMO_LHS_PATH)          # 121321 → … → 321323   (8 words)
    rhs_fwd = _ZAMO_RHS_PATH                     # 321323 → … → 121321   (8 words)
    # drop the shared 321323 at the seam AND the final 121321 (== word 1, closes the
    # cycle): 8 + 7 − 1 = 14 distinct words.
    return vcat(lhs_back, rhs_fwd[2:(end - 1)])  # 8 + 6 = 14 words
end

# the braid move (block) turning word `k` into word `k+1`, CYCLICALLY (14 → 1).
function _zamo_move(k::Int)
    ws = zamo_words()
    n = length(ws)                               # 14
    a = ws[mod1(k, n)]; b = ws[mod1(k + 1, n)]
    for (i, j) in braid_move_sites(a)
        apply_braid_move(a, i, j) == b && return (i, j)
    end
    error("no single braid move takes word $(mod1(k,n)) ($(join(a))) to $(join(b))")
end

"""
    Zamo(i, j) -> MorphismGraph

Build the morphism that walks the Zamolodchikov cycle FORWARD from word `i` to word `j`,
starting bottom-to-top from `zamo_words()[i]` (as identity strands) and joining one
braid vertex per step (via `braid_top`). Result: bottom = `zamo_words()[i]`, top =
`zamo_words()[j]`.

The cycle has **14 words**, indexed cyclically (word 14 → word 1 is a step too):

* `i` and `j` may be any integers; they are read mod 14 onto `1..14`.
* the walk always goes FORWARD; `i > j` simply wraps around (e.g. `Zamo(10, 6)` takes
  `mod(6-10, 14) = 10` steps).
* **`Zamo(i, i)` is the WHOLE cycle** (14 steps, back to word `i`) — NOT the identity.

Examples: `Zamo(1, 2)` is one step; `Zamo(1, 8)` the LHS path `121321 → 321323`;
`Zamo(8, 1)` the RHS path back `321323 → 121321`; `Zamo(1, 1)` the full 14-step cycle
(bottom = top = 121321).
"""
function Zamo(i::Int, j::Int)
    ws = zamo_words()
    n = length(ws)                               # 14
    nsteps = mod(j - i, n)
    nsteps == 0 && (nsteps = n)                  # Zamo(i,i) = the whole cycle
    f = light_leaf_up(copy(ws[mod1(i, n)]))      # identity strands on the start word
    for s in 0:(nsteps - 1)
        (bi, bj) = _zamo_move(i + s)             # _zamo_move is cyclic
        f = braid_top(f, bi, bj)
    end
    return f
end

# build a morphism bottom→top from a list of words `path` (consecutive single moves),
# joining one braid vertex per step. bottom = path[1], top = path[end].
function _build_path(path)
    f = light_leaf_up(copy(path[1]))
    for k in 1:(length(path) - 1)
        a, b = path[k], path[k + 1]
        mv = nothing
        for (i, j) in braid_move_sites(a)
            apply_braid_move(a, i, j) == b && (mv = (i, j); break)
        end
        mv === nothing && error("_build_path: no single move $(join(a)) → $(join(b))")
        f = braid_top(f, mv[1], mv[2])
    end
    return f
end

"""
    zamo_lhs() -> MorphismGraph
    zamo_rhs() -> MorphismGraph

The two sides of the Zamolodchikov relation, both morphisms `121321 → 321323` (built
bottom-to-top from the canonical `121321`): `zamo_lhs()` is the LHS path
(`= Zamo(1, 8)`), `zamo_rhs()` is the RHS path over the OTHER half of the hexagon. The
Zamolodchikov relation asserts these two morphisms are EQUAL. Here we can only check
they have the same bottom and top (`zamo_relation_endpoints_agree`); proving the
MORPHISMS equal needs diagram simplification (a later step).
"""
zamo_lhs() = Zamo(1, 8)

# RHS half: 121321 → … → 321323, i.e. the RHS path (_ZAMO_RHS_PATH: 321323→…→121321)
# read BACKWARD so it starts at 121321.
zamo_rhs() = _build_path(reverse(_ZAMO_RHS_PATH))

"""
    zamo_relation_endpoints_agree() -> Bool

Combinatorial sanity check of the Zamolodchikov relation: LHS and RHS have the same
bottom (`121321`) and the same top (`321323`). This is NECESSARY but not sufficient —
morphism equality is proved later by diagram simplification.
"""
zamo_relation_endpoints_agree() =
    bottom(zamo_lhs()) == bottom(zamo_rhs()) && top(zamo_lhs()) == top(zamo_rhs())
