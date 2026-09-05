# circular/CircularGraph.jl — the CircularNode/CircularGraph type.
#
# WHY. When reducing, diagrams do not end up in DL form. Two causes:
#   (a) `111 → 1` can be wired from trivalents in TWO ways (the two
#       bracketings). By associativity they are equal, but PROVING that is
#       expensive and fragile.
#   (b) 1/3 combinations are isomorphic in a way the plain node kinds
#       (`:dot`/`:trivalent`/`:braid` with fixed degrees 1/3/2m) cannot see —
#       they only know ONE colour resp. a fixed colour pair, not "these arms
#       in this order".
#
# THE SOLUTION: a node may have ARBITRARILY many arms, and its only
# information is the cyclic arm-colour sequence. A cluster of several nodes
# becomes ONE node by merging (`merge_nodes`) — and the two bracketings of
# `111 → 1` thereby become literally the same node, with nothing left to
# check.
#
# Merging is NOT forced in the constructor (it is a reduction rule instead) —
# unmerged CircularGraphs must stay representable.
#
# REUSE: `Port`, `Leaf`, `NodePort`, `Circle`, `Edge`, `Cells` and
# `CircularWord` come unchanged from diagram/Graph.jl — none of these types
# knows about `Node`, so there is nothing to duplicate. `CircularNode` and
# `CircularGraph` are the only new types.

# ---- the node type -------------------------------------------------------------

"""
    CircularNode

A node with ARBITRARILY many arms. Its only information is `arms`: the
colours of slots `1..deg` in CLOCKWISE order — the same convention as the plain
`Node` (diagram/Graph.jl §PLANAR SLOT CONVENTION, diagram/Faces.jl
§ROTATIONSSYSTEM). `kind` is derivable from `arms` (see `circular_node`) and serves
only as a fast dispatch/validation tag:

- `:mono`  — arbitrarily many arms, all colour 2. **Only** colour 2: colours 1
  and 3 belong to the 1/3 world and are therefore `:mixed`, not `:mono`.
- `:mixed` — the 1/3 world: an arbitrary sequence of 1s and 3s. A pure 1-node
  `[1,1,1]` is the SPECIAL CASE here, not its own kind — that is exactly what
  turns `111 → 1` into one node kind instead of two.
- `:braid` — the 12-/23-node: `2k` arms with `k ≥ 3`, alternating `s,t,s,t,…`
  with `{s,t} = {1,2}` or `{2,3}`. The **ARM COUNT IS AN ATTRIBUTE**, not a
  separate kind: `2k = 6` is the ordinary m=3 braid, `2k = 8` the
  "general-12-node" (the cluster of 2 braids + 2 trivalents), larger `2k` the
  tower family `π_k`. What depends on the arm count lives in FUNCTIONS instead
  of the type: `circular_degree = 6 − 2k`, `circular_opposite_slot = slot + k`,
  and C6–C8 explicitly check `arm_count == 6`. C15/C16 take any arm count —
  their surgery derives everything from the two counts. What still stops at 8
  is the PREIMAGE side: `GBRAID_PREIMAGES` and `expand_gbraid` know only the
  8-armed node, and C14 is missing for `k ≥ 5`. The m=2 commutation `braid(1,3)`, by
  contrast, is JUST a 4-armed `:mixed` node, not `:braid` — see CircularConvert.jl.

NO provenance/history field: that would destroy exactly the identification
this type is built for.

EQUALITY holds only up to ROTATION (oriented, clockwise), NOT up to
reflection — see `circular_key`.
"""
struct CircularNode
    kind::Symbol
    arms::Vector{Int}

    function CircularNode(kind::Symbol, arms::Vector{Int})
        isempty(arms) && throw(ArgumentError("CircularNode needs at least one arm"))
        if kind === :mono
            all(==(2), arms) ||
                throw(ArgumentError(":mono nodes may only carry colour 2, got $arms"))
        elseif kind === :mixed
            all(c -> c == 1 || c == 3, arms) ||
                throw(ArgumentError(":mixed nodes may only carry colours 1/3, got $arms"))
            # If both colors appear, require at least two occurrences of each to avoid
            # spurious too-small mixed nodes like [1,3,1] or [1,3]. Single-colour
            # mixed nodes (e.g. [1] or [1,1,1]) remain allowed HERE — the separate
            # "no arm_count==2" ban (2-armed "beads") lives only in
            # `circular_node`, NOT in this inner constructor;
            # [1,1] still passes this constructor (test fixtures/rule internals
            # need that transiently), even though `circular_node([1,1])` throws.
            c1 = count(==(1), arms)
            c3 = count(==(3), arms)
            if c1 > 0 && c3 > 0
                (c1 >= 2 && c3 >= 2) ||
                    throw(ArgumentError(":mixed nodes with two colours need each colour at least twice, got $arms"))
            end
        elseif kind === :braid
            # ONE type covers 6 and 2k arms: the same three invariants hold
            # regardless of arm count, with `k ≥ 3`. The arm count is an
            # ATTRIBUTE, not its own type.
            (iseven(length(arms)) && length(arms) >= 6) ||
                throw(ArgumentError(":braid nodes need 2k arms with k >= 3, got $(length(arms))"))
            s, t = arms[1], arms[2]
            (Set((s, t)) == Set((1, 2)) || Set((s, t)) == Set((2, 3))) ||
                throw(ArgumentError(":braid colour pair must be {1,2} or {2,3}, got {$s,$t}"))
            arms == [isodd(k) ? s : t for k in 1:length(arms)] ||
                throw(ArgumentError(":braid arms must alternate, got $arms"))
        else
            throw(ArgumentError("unknown CircularNode kind :$kind"))
        end
        return new(kind, arms)
    end

    # Unchecked back door for the reduction rule hot path: declared as `global`
    # INSIDE the struct block so it can reach the implicit `new`. Only for code
    # that guarantees the invariants itself.
    global _unchecked_circularnode(k::Symbol, a::Vector{Int}) = new(k, a)
end

"""
    circular_node(arms) -> CircularNode

Derives `kind` from the arm colours instead of requiring it separately: all
arms `2` ⇒ `:mono`; `2k ≥ 6` arms alternating with colour pair `{1,2}`/`{2,3}`
⇒ `:braid` (any arm count, see `CircularNode`); all arms in `{1,3}` ⇒ `:mixed`;
otherwise `ArgumentError`.

Additionally checks the arm-count invariant: a
CircularNode has EXACTLY 1 arm (dot) or AT LEAST 3 arms — `arm_count == 2` (a
"bead") is invalid regardless of colours (so both `[1,1]` and `[1,3]` are
forbidden). This check sits only here, in the convenient public constructor,
and NOT in the inner `CircularNode(kind, arms)` resp. `_unchecked_circularnode`:
rule internals (e.g. `_circular_drop_adjacent_selfloop`,
`src/circular/rules/CircularGen12Merge.jl`) and test fixtures build transient 2-armed
intermediate nodes that disappear only in the next step (collapse resp.
wiring) — `circular_node` throws so a bead never appears as an END state;
the internal paths remain unhindered.
"""
function circular_node(arms::Vector{Int})
    isempty(arms) && throw(ArgumentError("CircularNode needs at least one arm"))
    length(arms) == 2 &&
        throw(ArgumentError("circular_node: a CircularNode has 1 arm (dot) or at least 3 arms, got 2"))
    if all(==(2), arms)
        return CircularNode(:mono, arms)
    elseif length(arms) >= 6 && iseven(length(arms)) &&
           (Set((arms[1], arms[2])) == Set((1, 2)) || Set((arms[1], arms[2])) == Set((2, 3))) &&
           arms == [isodd(k) ? arms[1] : arms[2] for k in 1:length(arms)]
        # ONE branch covers 6 and 2k arms.
        return CircularNode(:braid, arms)
    elseif all(c -> c == 1 || c == 3, arms)
        # If both colors appear, require at least two occurrences of each to avoid
        # spurious too-small mixed nodes like [1,3,1] or [1,3]. Single-colour
        # mixed nodes (e.g. [1] or [1,1,1]) remain allowed.
        c1 = count(==(1), arms)
        c3 = count(==(3), arms)
        if c1 > 0 && c3 > 0
            (c1 >= 2 && c3 >= 2) || throw(ArgumentError(":mixed nodes with two colours need each colour at least twice, got $arms"))
        end
        return CircularNode(:mixed, arms)
    else
        throw(ArgumentError("arm sequence $arms fits no CircularNode kind"))
    end
end

# ---- the graph type ------------------------------------------------------------
#
# Placed before the accessor functions because `arm_count(g::CircularGraph, i)`
# needs to know the type.

"""
    CircularGraph

Second, PARALLEL type to `WordGraph` — the existing `WordGraph` stack stays
untouched. `word`/`edges` are literally the same types as there; only the
nodes are `CircularNode` instead of `Node`. The 3-argument constructor in
circular/CircularFaces.jl fills `cells` via the same generic tracer as `WordGraph`
(diagram/Faces.jl) — and, as there, it is a pure function of
`(word, nodes, edges)`, so it plays no role in `==`/`hash`.
"""
struct CircularGraph
    word::CircularWord
    nodes::Vector{CircularNode}
    edges::Vector{Edge}
    cells::Cells
    # Lazy cache for `circular_canonical_key`: `==`/`hash` pull the key on every
    # comparison/lookup, which without caching costs ~80% of warm reduction
    # time. The field is a `Ref` so the immutable graph can memoise it on first
    # access; it carries NO information (a pure function of the other fields),
    # so it stays out of `==`/`hash`/`show`.
    _ckey::Base.RefValue{Any}

    CircularGraph(word::CircularWord, nodes::Vector{CircularNode}, edges::Vector{Edge}, cells::Cells) =
        new(word, nodes, edges, cells, Ref{Any}(nothing))
end

Base.length(g::CircularGraph) = length(g.word)

"The boundary word of `g` (analogue of `boundary` for `WordGraph`)."
circular_boundary(g::CircularGraph) = g.word

"Colour of the edge attached to boundary `Leaf(k)` in a `CircularGraph` (0 if none)."
function leaf_colour(g::CircularGraph, k::Int)
    for e in g.edges
        e.a isa Leaf && e.a.k == k && return e.colour
        e.b isa Leaf && e.b.k == k && return e.colour
    end
    return 0
end

# ---- accessors ------------------------------------------------------------------

"The arm sequence of `nd` (slots `1..deg` clockwise). Do NOT mutate."
arms(nd::CircularNode) = nd.arms

"Degree of `nd` — the function the generic tracer needs."
_circular_degree(nd::CircularNode) = length(nd.arms)

"""
    arm_count(nd::CircularNode) -> Int
    arm_count(g::CircularGraph, i::Int) -> Int

The ARM COUNT of `nd` (resp. of node `i` in `g`) — needed structurally for the
slots (tracer, rotation system, renderer). Arm count and `circular_degree`
(see below) are two different numbers (e.g. dot: 1 arm, degree +1;
trivalent: 3 arms, degree −1).
"""
arm_count(nd::CircularNode) = length(nd.arms)
arm_count(g::CircularGraph, i::Int) = arm_count(g.nodes[i])

"""
    _circular_braidlike(nd::CircularNode) -> Bool

The nodes the merge machinery must NOT touch: edge contraction with a
`:mono`/`:mixed` neighbour would destroy the alternating arm sequence that is
their invariant.

This is literally `kind === :braid`. The function stays under its own name
because it names the INTENT ("braid-like, hence locked for merging") and
because it is the one place a possible third braid-like kind would need to be
added.
"""
_circular_braidlike(nd::CircularNode) = nd.kind === :braid

"""
    circular_degree(nd::CircularNode) -> Int
    circular_degree(g::CircularGraph) -> Int

The DEGREE of `nd`:

    circular_degree(nd) = 2·(number of DISTINCT colours in `arms(nd)`) − arm_count(nd)

for `:mono` and `:mixed`. The braid-like nodes (`_circular_braidlike`) have their
OWN formula, `6 − arm_count`: for the 6-armed braid this is the fixed value
`0`, from 8 arms on it is the tower computation (see below), stated ONCE
instead of per case. Not a special case of the formula above, but an
exception: circular nodes only arise from COMMUTING colours (the 1/3
world and mono-2 — differently-coloured arms of the same node only come about
because their colours commute), whereas the `(1,2)`/`(2,3)` braid (m=3) is a
fixed, non-mergeable shape — nodes like `[1,2,1,2,2,1]` do not structurally
exist (see `CircularNode`'s constructor invariants). The braid is thus the ONLY
node kind the formula does not have to handle, because its degree never arises
from merging.

Exactly reproduces the plain `degree(::WordGraph)` (`dot_count − trivalent_count`,
`MorphismGraph.jl:139`) on every node kind it shares with that type:

| node | `circular_degree` | plain `degree` |
|---|---|---|
| dot `[1]`, `[2]` | +1 (2·1−1) | +1 |
| trivalent `[1,1,1]`, `[2,2,2]` | −1 (2·1−3) | −1 |
| crossing `[1,3,1,3]` (m=2 commutation) | 0 (2·2−4) | 0 |
| braid `[1,2,1,2,1,2]` (m=3) | 0 (6−6) | 0 |
| general-12 `[1,2,1,2,1,2,1,2]` (k=4) | −2 (6−8) | — (cluster: 2 braids + 2 trivalents = 0+0−1−1) |
| `[1,1,1,1]` (merged `111→1`) | −2 (2·1−4) | — (no plain equivalent) |
| `[1,3,3,1,3,3]` | −2 (2·2−6) | — (no plain equivalent) |

**The tower family `π_k` has degree `6 − #arms = −(2k−6)`**: `π_k` is
`BS(x) ↠ B_{w₀}(extreme shift) ↪ BS(y)`, and this degree space is
one-dimensional. The `2k`-node IS the tower, so its degree is the sum of its
parts. The formula agrees with `2 − arm_count/2` exactly at 8 arms and is
independently confirmed for every other arm count, see below.

Both formulas agree at **8** arms (`6−8 == 2−8/2 == −2`). From 10 arms on they
diverge, and there `2 − arm_count/2` is independently refuted: for **every**
circular diagram, `circular_degree ≡ #boundary arms (mod 2)` holds (dot 1 arm/`+1`,
trivalent 3/`−1`, crossing 4/`0`, braid 6/`0`; every inner edge consumes two
arms) — so a 10-armed node of degree `−3` cannot exist.

Only with this formula does the **fusion balance** of C13 work out
(`_circular_gen12_clusters`): `2k` and `2l`, triply connected with two
trivalents outside, give `2(k+l−2)`, and
`(6−2k) + (6−2l) + (−1) + (−1) = 6 − 2(k+l−2)` ✓.

Important for `test/pairing.jl`: the Soergel grading there hangs on the Hecke
algebra ground truth — the formula must not let anything slip there, and the
table above is exactly the check for that.
"""
circular_degree(nd::CircularNode) =
    _circular_braidlike(nd) ? 6 - length(nd.arms) :
    2 * length(unique(nd.arms)) - length(nd.arms)

circular_degree(g::CircularGraph) = sum(circular_degree, g.nodes; init = 0)

"""
    CIRCULAR_DOT_MIN_ARMS

The MINIMAL arm count per colour pair — `{1,2}` and `{2,3}` have 6, `{1,3}`
has 4. Used by [`circular_dot_may_pass`](@ref).
"""
const CIRCULAR_DOT_MIN_ARMS = Dict(Set([1, 2]) => 6, Set([2, 3]) => 6, Set([1, 3]) => 4)

"""
    CIRCULAR_DOT_POLICY

Is the dot policy ([`circular_dot_may_pass`](@ref)) in force? Default
**`true`**. ONE switch for both places it acts:

* the RULES — C14 `_fr_dot_on_gen12` and the generic preimage stage decline to
  push a dot through a node it may not pass (circular/rules/CircularDotOnGen12.jl,
  circular/rules/CircularDriver.jl);
* the CRITERION — `circular_unreduced_transition` drops the `:dot` transitions
  such a dot would raise (circular/rules/CircularRegionWord.jl).

Both call sites take a `policy` keyword resp. read this `Ref`, so a
measurement can switch it off without editing the source: calling C14 with
`policy = false` shows what pushing through WOULD give.
"""
const CIRCULAR_DOT_POLICY = Ref(true)

"""
    circular_dots_reducible(g::CircularGraph, v::Int) -> Bool

May the dots on node `v` be pushed through it? **This is the gate the rules
use** ([`CIRCULAR_DOT_POLICY`](@ref)); [`circular_dot_may_pass`](@ref) is only
its first clause.

Dots on big braid nodes are allowed in general; a simplification exists
exactly when two dots sit next to each other, or when so many arms carry
dots that fewer than three strands of one colour remain.

So on a braid-like node with `2k ≥ 8` arms the dots stay put UNLESS

* **(a)** two dotted arms are ADJACENT — then C20 `_fr_two_adjacent_dots`
  collapses the node to `2(k−1)` arms in one step, or
* **(b)** after removing the dotted arms, one colour is left with **fewer than
  three** arms — the node is then no longer a genuine `2k`-braid (the minimal
  one, 6 arms, has 3 of each colour) and collapses.

Everything else — a single dot, or dots that leave both colours at ≥ 3 and are
never neighbours — is a NORMAL FORM: the dot sits on the big node and stays
there. That is what makes "10 arms with two opposite dots" irreducible, while
"8 arms with two adjacent dots" is not.

At a MINIMAL node (`circular_dot_may_pass`: 6 arms for `{1,2}`/`{2,3}`, 4 for
`{1,3}`) and at every non-braid node the answer is always `true` — those are
C5/C6 territory and unchanged.
"""
function circular_dots_reducible(g::CircularGraph, v::Int)
    nd = g.nodes[v]
    circular_dot_may_pass(nd) && return true
    (nd.kind === :braid && arm_count(nd) >= 8) || return true
    n = arm_count(nd)
    dotted = falses(n)
    for e in g.edges
        (e.a isa NodePort && e.b isa NodePort) || continue
        e.a.node == v && arm_count(g.nodes[e.b.node]) == 1 && (dotted[e.a.slot] = true)
        e.b.node == v && arm_count(g.nodes[e.a.node]) == 1 && (dotted[e.b.slot] = true)
    end
    any(dotted) || return false                                   # no dot: nothing to pull
    any(dotted[j] && dotted[mod1(j + 1, n)] for j in 1:n) && return true          # (a)
    for c in unique(nd.arms)                                                      # (b)
        count(s -> nd.arms[s] == c && !dotted[s], 1:n) < 3 && return true
    end
    return false
end

"""
    circular_dot_may_pass(nd::CircularNode) -> Bool

May a dot be pushed THROUGH the node `nd`? Only through the **minimal** node of
its colour pair: 6 arms for `{1,2}`/`{2,3}`, 4 arms for `{1,3}`.

Dots are pulled only through 12-, 23-, 13-nodes with the minimal number of
arms, i.e. 6, 6, 4.

⚠ **This predicate alone is NOT the gate.** Dots on big braid nodes are also
allowed, with the two cases in which a bigger node DOES collapse named
explicitly. The gate the rules read is therefore
[`circular_dots_reducible`](@ref), whose first clause is this predicate. Kept
as its own name because "minimal node of its colour pair" is the notion both
clauses are written against.
"""
function circular_dot_may_pass(nd::CircularNode)
    colours = Set(unique(nd.arms))
    haskey(CIRCULAR_DOT_MIN_ARMS, colours) || return false
    return arm_count(nd) == CIRCULAR_DOT_MIN_ARMS[colours]
end

"""
    arm_colour(nd::CircularNode, slot::Int) -> Int

Colour at `slot`, read cyclically (`mod1`): slot `deg+1` is slot 1 again.
"""
arm_colour(nd::CircularNode, slot::Int) = nd.arms[mod1(slot, length(nd.arms))]

"""
    rotate_arms(nd::CircularNode, r::Int) -> CircularNode

Rotates the arm sequence by `r`: slot `s` of the result is slot `s+r` of `nd`.
Since equality holds only up to rotation anyway, `rotate_arms(nd, r) == nd`
for every `r` — useful for testing exactly that.
"""
function rotate_arms(nd::CircularNode, r::Int)
    d = length(nd.arms)
    return CircularNode(nd.kind, [arm_colour(nd, s + r) for s in 1:d])
end

"""
    rotate_arms(g::CircularGraph, node::Int, r::Int) -> CircularGraph

Graph-level version of `rotate_arms(nd::CircularNode, r)`. Rotates the arm
sequence of `g.nodes[node]` by `r` AND rewrites the edges hanging off it
accordingly, so the same physical wiring is preserved: an edge hanging at
slot `p` before the rotation hangs at slot `mod1(p - r, deg)` afterwards —
the slot whose rotated colour carries `p`'s original colour again. Since
`CircularNode` equality holds only up to rotation, the result represents the
SAME diagram as `g` — the tool for testing rotation invariance (e.g. of
`circular_canonical_key`).
"""
function rotate_arms(g::CircularGraph, node::Int, r::Int)
    nd = g.nodes[node]
    d = arm_count(nd)
    new_nodes = copy(g.nodes)
    new_nodes[node] = rotate_arms(nd, r)
    remap(p::Int) = mod1(p - r, d)
    new_edges = Edge[]
    for e in g.edges
        a = e.a isa NodePort && e.a.node == node ? NodePort(node, remap(e.a.slot)) : e.a
        b = e.b isa NodePort && e.b.node == node ? NodePort(node, remap(e.b.slot)) : e.b
        push!(new_edges, Edge(e.colour, a, b))
    end
    return CircularGraph(g.word, new_nodes, new_edges)
end

"""
    circular_opposite_slot(nd::CircularNode, slot::Int) -> Int

The slot OPPOSITE `slot` at a braid-like node with `2k` arms: `mod1(slot + k,
2k)`; for the 6-armed braid this specialises to `mod1(slot + 3, 6)`. Analogue
of `opposite_slot` in diagram/Graph.jl.

Only defined for braid-like nodes (`_circular_braidlike`): a `:mono`/`:mixed` node
has no distinguished opposite — its arms don't alternate, so "opposite" would
have no meaning there. Other kinds throw.
"""
function circular_opposite_slot(nd::CircularNode, slot::Int)
    _circular_braidlike(nd) ||
        error("circular_opposite_slot only defined for braid-like nodes")
    n = length(nd.arms)
    return mod1(slot + n ÷ 2, n)
end

