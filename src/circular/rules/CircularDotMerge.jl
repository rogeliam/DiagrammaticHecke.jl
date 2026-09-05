# circular/rules/CircularDotMerge.jl — three rules
# about a DOT on a braid-like `2k`-armed node.
#
# THE BACKGROUND. The dot policy `circular_dot_may_pass` (circular/CircularGraph.jl)
# says a dot may only be pushed through the MINIMAL node of its colour pair
# (6 arms for {1,2}/{2,3}, 4 for {1,3}); on a bigger node it stays and counts as
# reduced. That leaves two holes, and this file closes both:
#
#   C19 `trivalent_into_gen12` — the INVERSE of C14 `dot_on_gen12`:
#       trivalent (colour c') on an arm of a `2k`-node  ⟹  `2(k+1)`-node whose
#       new middle arm (the other colour c) carries a DOT: a
#       trivalent plus a general-12 node with `2k` arms merges into a
#       `2(k+1)`-armed node with a dot of the other colour.
#       It RAISES the BRAID arm sum, so it does not fall under
#       `circular_arm_weight`. It falls under `circular_sep_weight`, where the
#       second component counts separating edges and swallowing the trivalent
#       turns one of them into a dot edge — which is why it IS registered, and
#       why C14 is not (circular/rules/CircularWeight.jl, second header block).
#       Verified: two C19 steps turn
#       „trivalent + 6-armed + trivalent" into „10-armed + two opposite dots",
#       bit-identical (`circular_canonical_key`) to the C14 round trip backwards.
#
#   C20 `two_adjacent_dots` — two dots on ADJACENT arms of a `2k`-node give the
#       naked `2(k-1)`-node. This is C14 followed by C5 `merge` in one step
#       (measured for 8 and 10 arms), and WITHOUT it the dot
#       policy is not confluent: „8-armed + two adjacent dots" and „naked
#       6-armed" would both be normal forms of the same morphism, and neither
#       the region count nor the number of separating edges tells them apart
#       (6/6 and 6/6). C20 lowers arm sum AND node count, so it is
#       safe for the present weight and IS registered.
#
#   C23 `mixed_into_gen12` — a mixed {1,3} node with at least three arms of
#       the braid colour, hanging on a `2k`-node, is broken open and the freed
#       trivalent merged in: a `2(k+1)`-node with a dot, plus the two pieces of
#       the mixed node. This is C19 applied to a trivalent the rule has to
#       manufacture first, and it is ONE step so that the intermediate figure
#       never has to be a normal form. It raises the arm sum, and unlike C19
#       no reading of the weight catches it, so it carries a termination
#       argument of its own and stands in `CIRCULAR_WEIGHT_ASSERT_EXEMPT`.
#       Full statement in the box before `CIRCULAR_MIXED_MERGE` below.

"""
    _fr_two_adjacent_dots(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C20 `two_adjacent_dots`.** Two dots on adjacent arms of a braid-like node with
`2k ≥ 6` arms ⟹ the naked `2(k−1)`-armed node (for `2k = 6`: the `{1,3}`/mono
case cannot occur, `:braid` nodes always alternate two colours). One term,
coefficient 1. Degree check: `(6−2k) + 1 + 1 = 6 − 2(k−1)`.

See the file header for why this rule has to exist next to the dot policy.
"""
function _fr_two_adjacent_dots(g::CircularGraph)
    for (i, nd) in enumerate(g.nodes)
        (nd.kind === :braid && arm_count(nd) >= 8) || continue
        n = arm_count(nd)
        # dot node hanging off slot s?
        dot_at = Dict{Int,Int}()
        for e in g.edges
            (e.a isa NodePort && e.b isa NodePort) || continue
            e.a.node == i && arm_count(g.nodes[e.b.node]) == 1 && (dot_at[e.a.slot] = e.b.node)
            e.b.node == i && arm_count(g.nodes[e.a.node]) == 1 && (dot_at[e.b.slot] = e.a.node)
        end
        length(dot_at) >= 2 || continue
        for j in 1:n
            jp = mod1(j + 1, n)
            (haskey(dot_at, j) && haskey(dot_at, jp)) || continue

            a    = collect(arms(nd))
            keep = [mod1(jp + t, n) for t in 1:(n - 2)]        # cyclic, order preserved
            newnode = circular_node([a[s] for s in keep])

            weg = Set([i, dot_at[j], dot_at[jp]])
            old = [v for v in 1:length(g.nodes) if !(v in weg)]
            idx = Dict(v => t for (t, v) in enumerate(old))
            nodes = vcat(CircularNode[g.nodes[v] for v in old], CircularNode[newnode])
            NI = length(old) + 1
            slotmap = Dict(keep[t] => t for t in 1:(n - 2))

            function shift(p)
                p isa Leaf && return p
                p.node in weg || return NodePort(idx[p.node], p.slot)
                p.node == i && haskey(slotmap, p.slot) && return NodePort(NI, slotmap[p.slot])
                return nothing                                  # the two dot arms drop out
            end

            edges = Edge[]
            for ed in g.edges
                pa = shift(ed.a); pb = shift(ed.b)
                (pa === nothing || pb === nothing) && continue
                push!(edges, Edge(ed.colour, pa, pb))
            end
            return CircularComboR(CircularGraph(g.word, nodes, edges))
        end
    end
    return nothing
end

"""
    _fr_trivalent_into_gen12(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C19 `trivalent_into_gen12`**, the inverse of C14 `_fr_dot_on_gen12`
(circular/rules/CircularDotOnGen12.jl): a trivalent of colour `c'` sitting on an
arm of a braid-like `2k`-node is absorbed — the arm becomes three arms
`(c', c, c')`, the two outer ones take over the trivalent's free legs, and the
middle one (the OTHER colour `c`) gets a new DOT. Result: a `2(k+1)`-armed node.
One term, coefficient 1. Degree check: `−1 + (6−2k) = (6−2(k+1)) + 1`.

Registered in `CIRCULAR_RULES` behind `CIRCULAR_TRIVALENT_MERGE`, which is also
what takes C14 out: the two are inverse, so at most one may run. C19 raises the
BRAID arm sum, so it does not fall under `circular_arm_weight`; it falls under
`circular_sep_weight`, because swallowing the trivalent turns a separating edge
into a dot edge (CircularWeight.jl, second header block).

This is the matcher only: it looks for the pair (braid-like node, trivalent
joined to it) and hands it to [`_circular_mono_into_gen12_at`](@ref), which
does the surgery and owns the cyclic-order and loop-guard conditions.
"""
function _fr_trivalent_into_gen12(g::CircularGraph)
    for (i, nd) in enumerate(g.nodes)
        nd.kind === :braid || continue
        for (j, td) in enumerate(g.nodes)
            # A one-coloured node with `m ≥ 3` arms. `:braid` nodes alternate
            # two colours, so `unique == 1` already excludes them. `m == 3` is
            # the trivalent the rule is named after; the `m ≥ 4` half is the
            # `m`-valent merge, see `CIRCULAR_MVALENT_MERGE`.
            m = arm_count(td)
            (m >= 3 && length(unique(td.arms)) == 1) || continue
            (m == 3 || CIRCULAR_MVALENT_MERGE[]) || continue
            hit = nothing
            for e in g.edges
                (e.a isa NodePort && e.b isa NodePort) || continue
                (e.a.node == i && e.b.node == j) && (hit = (e.a.slot, e.b.slot))
                (e.b.node == i && e.a.node == j) && (hit = (e.b.slot, e.a.slot))
            end
            hit === nothing && continue
            gg = _circular_mono_into_gen12_at(g, i, j, hit[1], hit[2])
            gg === nothing && continue
            return CircularComboR(gg)
        end
    end
    return nothing
end


"""
    _circular_mono_into_gen12_at(g, i, j, conn, t0) -> Union{Nothing, CircularGraph}

The surgery behind C19 for ONE given pair: the braid-like node `i` swallows the
single-colour node `j` (`m ≥ 3` arms, colour `c'`) they share the edge
`(i, conn) — (j, t0)` on.

The one arm `conn` becomes `2m−3` arms `c' c c' … c c'`: the `m−1` arms of
colour `c'` take over `j`'s free legs in cyclic order, and each of the `m−2`
arms of the other colour `c` carries a new DOT. The result is a
`2(k+m−2)`-armed braid-like node. For `m = 3` that is one trivalent, one dot.

Degree check: `(2−m) + (6−2k) = (6−2(k+m−2)) + (m−2)`, both sides `8−2k−m`.
Arm check: `(2k−1) + (m−1)` free arms before, `2(k+m−2) − (m−2)` after.

`nothing` on a colour mismatch at `conn`, and `nothing` when the merged node
would hand its new dots straight back to C14 `_fr_dot_on_gen12`
(`circular_dots_reducible`) — building a node the dot policy immediately
undoes would let the two rules push each other in a circle.

The cyclic order is not free (as in C14): the legs of `j`, read clockwise from
the connector, must land on the new node's `c'` arms in that same order.
`check_wiring` is the test for it.

It is a function of its own, and not inlined into `_fr_trivalent_into_gen12`,
because C23 `_fr_mixed_into_gen12` needs exactly this step on the trivalent it
frees from a mixed node.
"""
function _circular_mono_into_gen12_at(g::CircularGraph, i::Int, j::Int, conn::Int, t0::Int)
    nd = g.nodes[i]; n = arm_count(nd)
    td = g.nodes[j]; m = arm_count(td); cp = td.arms[1]
    arm_colour(nd, conn) == cp || return nothing

    a = collect(arms(nd))
    c = a[mod1(conn + 1, n)]                      # the other colour
    # the one arm `conn` becomes 2m-3 arms: c' c c' c … c' — the m-1
    # `c'` arms carry the legs, the m-2 `c` arms carry the new dots.
    middle = Int[isodd(t) ? cp : c for t in 1:(2m - 3)]
    newnode = circular_node(vcat(a[1:conn-1], middle, a[conn+1:end]))

    old = [v for v in 1:length(g.nodes) if v != i && v != j]
    idx = Dict(v => t for (t, v) in enumerate(old))
    nodes = vcat(CircularNode[g.nodes[v] for v in old],
                 CircularNode[newnode],
                 CircularNode[circular_node([c]) for _ in 1:(m - 2)])
    NI = length(old) + 1                          # dots: NI+1 … NI+m-2
    new_slot(s) = s < conn ? s : s + (2m - 4)
    # leg r of `j` (clockwise from the connector) lands on the r-th `c'` arm
    # of the new node; for m = 3 that is `conn` and `conn+2`.
    leg = Dict(mod1(t0 + r, m) => conn + 2 * (r - 1) for r in 1:(m - 1))

    function shift(p)
        p isa Leaf && return p
        p.node == i && return NodePort(NI, new_slot(p.slot))
        if p.node == j
            haskey(leg, p.slot) && return NodePort(NI, leg[p.slot])
            return nothing                        # the connector itself
        end
        return NodePort(idx[p.node], p.slot)
    end

    edges = Edge[]
    for ed in g.edges
        (ed.a isa NodePort && ed.b isa NodePort &&
         ((ed.a.node == i && ed.b.node == j) || (ed.b.node == i && ed.a.node == j))) && continue
        pa = shift(ed.a); pb = shift(ed.b)
        (pa === nothing || pb === nothing) && continue
        push!(edges, Edge(ed.colour, pa, pb))
    end
    for t in 1:(m - 2)
        push!(edges, Edge(c, NodePort(NI, conn + 2t - 1), NodePort(NI + t, 1)))
    end
    gg = CircularGraph(g.word, nodes, edges)
    # The C14 loop guard: a node whose dots the dot policy pulls straight back
    # out is not a step forward, it is the other half of a cycle.
    circular_dots_reducible(gg, NI) && return nothing
    return gg
end


# ---- C23 `mixed_into_gen12`: the mixed {1,3} node at a braid ---------------
#
# THE EQUATION (one term, coefficient 1)
#
#   braid `{c, c'}` (2k arms) —c— mixed `{1,3}` node with ≥ 3 arms of colour c
#     ↦  braid with 2(k+1) arms (the connector arm becomes `c, c', c`, the new
#        `c'` arm carries a DOT) + TWO mixed nodes A, B hanging on the two new
#        `c` arms, joined to each other by an edge of the other commuting
#        colour `c̄` (the {1,3} partner of `c`).
#
# WHY IT LOOKS LIKE THAT. It is two moves in one. Pull a `c`-trivalent out of
# the mixed node: its legs run to the two remaining pieces A and B, which stay
# joined by the `c̄`-strand that used to cross the node. Then C19 merges the
# freed trivalent into the braid. The rule performs exactly that pair —
# `_circular_split_mixed_at`, then `_circular_mono_into_gen12_at` — as a single
# step, so that the intermediate figure never has to be a normal form.
#
# ≥ 3 ARMS OF COLOUR `c` is what makes the split legal: the trivalent takes one
# `c`-arm away and gives one `c`-leg to each piece, and a two-coloured mixed
# node needs every colour at least twice (`circular_node`). Below three there
# is nothing to break open.
#
# DEGREE. A mixed node has degree `2·2 − arms` (`13113` → −1, `1313` → 0). The
# split gives trivalent (−1) plus A and B with `a + b = m + 3` arms and one
# inner `c̄` edge: `(−1) + (8 − (m+3)) = 4 − m`, the degree of the node that
# went in. C19 is degree-neutral as well, so the balance closes with no extra
# polynomial.
#
# THE CUT. Reading the mixed node's arms cyclically after the connector,
# `w = (w₁ … w_{m−1})`: A takes `w₁ … w_t`, B the rest, with `t` the FIRST cut
# at which both pieces carry both colours. For `13113` that gives A = (3,1) and
# B = (1,3), i.e. two `1313` crossings on the legs of the trivalent. From six
# arms on, several cuts balance the degree equally well; the rule commits to
# the first admissible one so that its result is a function of its input.
#
# TERMINATION. `circular_weight` RISES here — the braid node grows by two arms
# — so C23 stands in `CIRCULAR_WEIGHT_ASSERT_EXEMPT` (CircularDriver.jl) and
# carries its own argument instead: every piece has strictly fewer arms than
# the node it came from (A takes at least two of them), so a piece with ≥ 3
# `c`-arms is broken open only finitely often and the chain ends at four arms.

"""
    CIRCULAR_MIXED_MERGE

Does C23 `_fr_mixed_into_gen12` run? Default **`true`**. Off, a mixed `{1,3}`
node stays whole next to a braid node, and the diagram keeps whichever of the
two shapes it was built in.

It is a switch rather than a plain rule because C23 is the one registered rule
whose termination does not come from `circular_weight`
(see `CIRCULAR_WEIGHT_ASSERT_EXEMPT`), so a measurement wants to be able to
take it out of the registry without editing the source.
"""
const CIRCULAR_MIXED_MERGE = Ref(true)

"""
    _circular_split_mixed_at(g, j, t0) -> Union{Nothing, CircularGraph}

Step 1 of C23: break the mixed node `j` open at its arm `t0` (colour `c`). A
`c`-trivalent takes over arm `t0` on its slot 1 and hangs the two pieces A
(slot 2) and B (slot 3) on its other two legs; A and B stay joined by one new
edge of the other commuting colour. The trivalent is placed at index `j`, A and
B are appended, so every other node index is unchanged.

`nothing` unless `j` is a two-coloured mixed node with at least three arms of
colour `c` and a cut leaving both pieces two-coloured — see the C23 box above
for why those are the conditions.
"""
function _circular_split_mixed_at(g::CircularGraph, j::Int, t0::Int)
    nd = g.nodes[j]
    nd.kind === :mixed || return nothing
    a = collect(arms(nd)); m = length(a)
    c = a[t0]
    count(==(c), a) >= 3 || return nothing
    cols = unique(a); length(cols) == 2 || return nothing
    cbar = only(x for x in cols if x != c)
    w = [a[mod1(t0 + r, m)] for r in 1:(m - 1)]          # the arms after the connector
    slot_w = [mod1(t0 + r, m) for r in 1:(m - 1)]
    # The cut: the FIRST `t` such that both pieces carry both colours. Each
    # piece gains one more `c` (its leg) and one more `cbar` (the AB edge), and
    # a two-coloured mixed node needs every colour at least twice
    # (`circular_node`).
    both(v) = any(==(c), v) && any(==(cbar), v)
    t = findfirst(t -> both(w[1:t]) && both(w[t+1:end]), 1:(m - 2))
    t === nothing && return nothing
    # A: (leg from T, w₁ … w_t, AB edge);  B: (AB edge, w_{t+1} … w_{m−1}, leg from T)
    A = circular_node(vcat([c], w[1:t], [cbar]))
    B = circular_node(vcat([cbar], w[t+1:end], [c]))
    T = circular_node([c, c, c])
    nA = length(g.nodes) + 1; nB = nA + 1
    nodes = vcat(CircularNode[v == j ? T : g.nodes[v] for v in 1:length(g.nodes)],
                 CircularNode[A, B])
    target = Dict{Int,NodePort}()                        # old slot of j -> new port
    target[t0] = NodePort(j, 1)
    for r in 1:t;         target[slot_w[r]] = NodePort(nA, 1 + r); end
    for r in (t+1):(m-1); target[slot_w[r]] = NodePort(nB, 1 + (r - t)); end
    shift(p) = (p isa NodePort && p.node == j) ? target[p.slot] : p
    edges = Edge[Edge(e.colour, shift(e.a), shift(e.b)) for e in g.edges]
    push!(edges, Edge(c, NodePort(j, 2), NodePort(nA, 1)))
    push!(edges, Edge(c, NodePort(j, 3), NodePort(nB, arm_count(B))))
    push!(edges, Edge(cbar, NodePort(nA, arm_count(A)), NodePort(nB, 1)))
    return CircularGraph(g.word, nodes, edges)
end

"""
    _fr_mixed_into_gen12(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C23 `mixed_into_gen12`**: a mixed `{1,3}` node with **at least three arms of
colour `c`**, hanging by a `c`-arm on a braid-like `{c, c'}` node, is broken
open (`_circular_split_mixed_at`) and the freed `c`-trivalent merged into the
braid (`_circular_mono_into_gen12_at`, C19). One term, coefficient 1. It reads
the same for `c = 1` at a `{1,2}` braid and for `c = 3` at a `{2,3}` braid.

Behind `CIRCULAR_MIXED_MERGE`. Registered in `CIRCULAR_RULES` and exempt from
the driver's termination assert — the box above this function says why, and
what carries the termination instead.
"""
function _fr_mixed_into_gen12(g::CircularGraph)
    CIRCULAR_MIXED_MERGE[] || return nothing
    for (i, nd) in enumerate(g.nodes)
        nd.kind === :braid || continue
        for (j, td) in enumerate(g.nodes)
            td.kind === :mixed || continue
            for e in g.edges
                (e.a isa NodePort && e.b isa NodePort) || continue
                hit = e.a.node == i && e.b.node == j ? (e.a.slot, e.b.slot) :
                      e.b.node == i && e.a.node == j ? (e.b.slot, e.a.slot) : nothing
                hit === nothing && continue
                conn, t0 = hit
                c = arm_colour(td, t0)
                arm_colour(nd, conn) == c || continue
                count(==(c), arms(td)) >= 3 || continue
                gs = _circular_split_mixed_at(g, j, t0)
                gs === nothing && continue
                isempty(check_wiring(gs)) || continue
                gg = _circular_mono_into_gen12_at(gs, i, j, conn, 1)
                gg === nothing && continue
                isempty(check_wiring(gg)) || continue
                return CircularComboR(gg)
            end
        end
    end
    return nothing
end


