# circular/rules/CircularRexFusion.jl — the rex-cycle relation APPLIED LOCALLY,
# joining four pieces into one rule:
#
#   transition (`circular_unreduced_transition`) → path (`circular_path_edges`)
#   → relation (`circular_rex_relation`) → splice (`circular_splice`) → sum.
#
# THE TRIGGER IS THE REGION-WORD CRITERION: where a region word is not reduced,
# the diagram must simplify. Whether `circular_dot_reducible` should gate it
# instead is not yet resolved; there is no driver hook and the switch
# defaults to OFF.
#
# THE CURVE. The transition is a region `R` with reduced region word `w` and an
# `s`-edge on its boundary leading one level further, to region `to` with the
# NON-reduced word `w·s`. So the cut curve of the splice is the tree path to
# `to` — the path to `R` plus the transition edge — and its colours are exactly
# `w·s`, the word the relation of step 1 lives on. That equality is CHECKED, not
# assumed (`q.edges == vcat(p.edges, T.edge)`).
#
# ⚠ THE EXISTENCE OF THAT PARTNER IS NOT A GATE. At `:dihedral` (`id_1212`, no
# other 2-edge at `R4`) and `:long` (`id_1213213`, none at `R7`) the two
# `s`-strands are one resp. two regions apart, and the braid move of the rex cycle
# `1212 → 2122 → 1212` is what brings them together — inside the relation, not in
# `g`. Demanding the partner in `g` would assume the outcome and make the step
# decline on exactly the figures it is for.
#
# THE FUSION PAIR: the
# two edges that fuse are (transition edge, the OTHER `s`-coloured edge on the
# boundary of the SAME region) — NOT `prev_edge`, the tree parent. Fusion merges
# two edges of the SAME colour; at `case == :direct` the parent edge is
# `s`-coloured too and both coincide (which is why the `:regionword` version of
# the 2parallel rule applies there), at `:dihedral`/`:long` it does
# not. The pair is not used for surgery here — the surgery happens INSIDE the
# derived relation, where step 1's own search picks the carrying pair — but its
# EXISTENCE is the gate: without a second `s`-edge at `R` there is nothing to
# fuse and the step declines.
#
# THE ANCHOR (why the assignment is not searched per term). `circular_splice`
# searches which end of each cut edge goes to the bottom side. All terms of ONE
# relation must be inserted the SAME way round, or they do not describe one
# local replacement. So the assignment is fixed ONCE on the IDENTITY on `w·s`,
# which must reproduce `g` exactly (canonical key) — the step-3 self-test in
# situ — and every term is then spliced with that pinned assignment. If the
# identity does not reproduce `g`, the step declines instead of guessing.
#
# THE STEP IS LABEL-FREE (same policy as CircularDotSlide):
# `all(isone, region_labels) && isone(outer_label)` on the input AND on every term
# of the relation; anything else is REFUSED, not guessed. Transferring labels
# through the splice (`_circular_transfer_labels` with `overrides`) is a separate
# question.
#
# TWO ROUTES, ONE STEP. A trivial cycle — `wi` with `w` ending on `i` — needs only
# plain fusion:
#
#   `case == :direct`  →  route `:fusion`   — the two `s`-edges at the transition
#                         region ALREADY lie next to each other, so the plain
#                         2parallel surgery does it: a two-term sum with
#                         `α_s/2`, no relation, no cache, no splice.
#   otherwise          →  route `:relation` — the derived rex-cycle relation,
#                         spliced into the cut curve (below).
#
# ⚠ At `case == :direct` the relation is NOT derivable at all — route B of the
# trivial rex cycle (zero
# braid moves) keeps the identity summand, so `circular_rex_relation([2], 2)`
# raises. That is not a gap: with zero braid moves there is nothing to derive,
# the fusion IS the step. The derivation error stays uncaught for the non-trivial
# cases: an underivable relation there is a finding, not a silent `nothing`.

"""
    CIRCULAR_REX_FUSION_ENABLED

Switch for [`circular_rex_fusion_step`](@ref), default **`true`**.

The driver calls it in `reduce_to_circular_leave`, in the `:after_d4` group: the
structural rules run first, then Zamo, then 2parallel and D4, then this step, then
`circular_dot_fusion_step` — each one going back to the head of the round. It sits
that late because it is the expensive one (derived relation, disk cache) and
should only run where the cheap rules find nothing.

The trigger is the REGION WORD (`circular_unreduced_transition`), which no other
rule looks at. Termination rests on the lexicographic measure (diagram degree, max
unreduced region word) enforced by the α-guard `rex_alpha_blocked` in the driver,
with `_circular_leave_growth_guard` as the backstop. The step is callable
independently of the switch.
"""
const CIRCULAR_REX_FUSION_ENABLED = Ref(true)

"""
    CIRCULAR_MONO_BREAK

Whether the fusion step may BREAK A ONE-COLOURED NODE OPEN in front of a splice
that would otherwise decline ([`_crf_mono_break`](@ref)). Default **`true`**.

The criterion reads the SATURATED region word, the curve reads real edges — and
where the curve crosses two arms of ONE one-coloured node it reads `…cc…` where
the criterion reads `…c…`. That is not a different place, it is the same place
read twice. Breaking the node in two makes the curve cross only the new
connecting edge, and the two readings agree again. C5 merges the pieces back in
the next round, so breaking is never a rule of its own.

What it buys, measured on degree 0 of `321232/321232`: with the break the 12
double leaves give 12 circular leaves and the coordinate matrix has kernel
dimension 0; without it they give 13, and the kernel is one-dimensional. A
kernel of 0 is the statement one wants — the circular leaves are then a basis,
not one element too many.

⚠ The derived caches do not know about this switch, so a measurement of it
needs a FRESH PROCESS and a cache directory of its own per setting
(`cl_basis(...; recompute = true, dir = ...)`). Flipping it inside one session
and recomputing measures the cache, not the rule.
"""
const CIRCULAR_MONO_BREAK = Ref(true)

"""
    _circular_split_mono_at(g::CircularGraph, i::Int, s::Int)
        -> Union{Nothing, CircularGraph}

Break the ONE-COLOURED node `i` (`:mono`, `n >= 4` arms, all of colour `c`) open
along its two cyclically adjacent arms `s` and `s+1`: a TRIVALENT `A` takes those
two arms (its slots 1, 2) and hangs on the rest of the node by a new `c`-edge
(its slot 3); the remaining `n − 2` arms plus that edge become the `(n−1)`-armed
node `B` (slot 1 = the new edge, then the old arms `s+2 … s−1` in their cyclic
order). `A` is placed at index `i`, `B` is appended.

This is the INVERSE of C5 (`_fr_merge`): the degree is unchanged
(`2 − n = (2 − 3) + (2 − (n − 1))`) and so is the region count — one node and one
edge more leave `V − E + F` alone. For `n = 4` both pieces are trivalents.

`nothing` if the node is not `:mono`, has fewer than 4 arms, or `s` is not a slot.
A cell trace that refuses throws; that is the caller's `try`.
"""
function _circular_split_mono_at(g::CircularGraph, i::Int, s::Int)
    (1 <= i <= length(g.nodes)) || return nothing
    nd = g.nodes[i]
    nd.kind === :mono || return nothing
    n = arm_count(nd)
    n >= 4 || return nothing
    (1 <= s <= n) || return nothing
    c = arms(nd)[1]
    s2 = mod1(s + 1, n)
    A = circular_node([c, c, c])                     # (arm s, arm s+1, connector)
    B = circular_node(fill(c, n - 1))                # (connector, arm s+2, …, arm s−1)
    nB = length(g.nodes) + 1
    nodes = vcat(CircularNode[v == i ? A : g.nodes[v] for v in 1:length(g.nodes)],
                 CircularNode[B])
    target = Dict{Int,NodePort}(s => NodePort(i, 1), s2 => NodePort(i, 2))
    for r in 1:(n - 2)
        target[mod1(s2 + r, n)] = NodePort(nB, 1 + r)
    end
    shift(p) = (p isa NodePort && p.node == i) ? target[p.slot] : p
    edges = Edge[Edge(e.colour, shift(e.a), shift(e.b)) for e in g.edges]
    push!(edges, Edge(c, NodePort(i, 3), NodePort(nB, 1)))
    return CircularGraph(g.word, nodes, edges)
end

"""
    _crf_mono_break(m::CircularMorphismGraph, p) -> Union{Nothing, CircularMorphismGraph}

The pre-step of [`CIRCULAR_MONO_BREAK`](@ref): find two CONSECUTIVE edges of the
path `p` that have the same colour and meet at two cyclically adjacent arms of
one `:mono` node with at least 4 arms, and split that node
([`_circular_split_mono_at`](@ref)). `nothing` if the switch is off, if there is
no such place, or if the split does not give a planar, wired graph.
"""
function _crf_mono_break(m::CircularMorphismGraph, p)
    CIRCULAR_MONO_BREAK[] || return nothing
    g = m.graph
    for i in 1:(length(p.edges) - 1)
        e1 = g.edges[p.edges[i]]; e2 = g.edges[p.edges[i + 1]]
        e1.colour == e2.colour || continue
        for pa in (e1.a, e1.b), pb in (e2.a, e2.b)
            (pa isa NodePort && pb isa NodePort) || continue
            pa.node == pb.node || continue
            pa.slot == pb.slot && continue
            nd = g.nodes[pa.node]
            (nd.kind === :mono && arm_count(nd) >= 4) || continue
            n = arm_count(nd)
            s = mod1(pa.slot + 1, n) == pb.slot ? pa.slot :
                mod1(pb.slot + 1, n) == pa.slot ? pb.slot : 0
            s == 0 && continue                       # arms not cyclically adjacent
            # the cell trace may refuse (the constructor throws) — that is a
            # "not here", not an error, exactly as inside `circular_splice`
            gs = try
                _circular_split_mono_at(g, pa.node, s)
            catch
                nothing
            end
            gs === nothing && continue
            (euler(gs) == 2 && isempty(check_wiring(gs))) || continue
            return CircularMorphismGraph(gs, m.cut1, m.cut2)
        end
    end
    return nothing
end

"""
    _crf_common_sides(g, curve, rel, n) -> Union{Nothing, Vector{Bool}}

The side assignment for the WHOLE relation: the first assignment under which
EVERY term of `rel.rhs` splices into `curve` ([`circular_splice_sides`](@ref),
intersected in mask order). `nothing` if a term has no admissible assignment at
all, if the terms do not agree on one, or if a term is not label-free / does not
carry the right cuts.

`n` is the length of the curve word; each term is read as a morphism `w → w`
with cuts `(0, n)`.

FAST PATH: `false…false` is tried on every term FIRST — that is the assignment
the identity anchor produces, so wherever the step works without the scan it
costs one splice per term, and the `2^k` scan runs only where that gives up.
"""
function _crf_common_sides(g::CircularGraph, curve::Vector{Int}, rel, n::Int)
    w = [g.edges[e].colour for e in curve]
    terms = CircularMorphismGraph[]
    for (d, _) in pairs_of(rel.rhs)
        (all(isone, d.region_labels) && isone(d.outer_label)) || return nothing
        tm = CircularMorphismGraph(d.graph, 0, n)
        (bottom(tm) == w && top(tm) == w) || return nothing
        push!(terms, tm)
    end
    null = Bool[false for _ in curve]
    isempty(terms) && return null                              # empty relation
    all(tm -> circular_splice(g, curve, tm; sides = null) !== nothing, terms) && return null
    common = nothing
    for tm in terms
        sd = circular_splice_sides(g, curve, tm)
        isempty(sd) && return nothing
        common = common === nothing ? sd : [x for x in common if x in sd]
        isempty(common) && return nothing
    end
    return common[1]
end

"""
    _crf_simple(g, cut) -> Bool

Is the cut curve SIMPLE — pairwise distinct edges and no node carrying two of
them (else `:node_shared`, the unfolding case)?

⚠ REPORTED, NOT GATING: a non-simple curve splices FINE. On L4 (`[10,9]`), L7
(`[4,3]`) and the rex-cycle figure of `(121, 2)` (`[1,2,3,10]`) the curve shares a
node and the identity still splices back to `g` exactly. The binding check is the
ANCHOR (file head), not this predicate; the predicate is a classification and is
handed back in the result.
"""
function _crf_simple(g::CircularGraph, cut::Vector{Int})
    allunique(cut) || return false
    seen = Set{Int}()
    for e in cut, p in (g.edges[e].a, g.edges[e].b)
        p isa NodePort || continue
        p.node in seen && return false
        push!(seen, p.node)
    end
    return true
end

"""
    _crf_direct_fusion(fdm, T, adj) -> Union{Nothing, CircularComboR}

The `case == :direct` route (file head): fuse the transition edge with the
`s`-edge it already shares the transition region with, exactly as
[`circular_2parallel_step`](@ref) does it — [`circular_2parallel_apply`](@ref)
plus the `α_s/2` override on each new region, a two-term sum.

NOT routed through [`find_circular_2parallel`](@ref): that search
carries its own condition 4 (`:d4` distance pattern resp. the `:regionword`
pinning) and answers a different question. Here the region-word criterion already
says that this place must be simplified, and at `:direct` the pair it names is the
pair to fuse. `nothing` if there is no partner or the surgery is not
constructible.
"""
function _crf_direct_fusion(fdm::CircularDecoratedMorphism, T, adj)
    g = fdm.m.graph
    isone(fdm.region_labels[T.region]) || return nothing     # condition 5, as in 2parallel
    partners = unique(Int[ei for (_, ei) in adj[T.region]
                          if ei != T.edge && g.edges[ei].colour == T.colour &&
                             !_circular_edges_touch(g.edges[ei], g.edges[T.edge])])
    isempty(partners) && return nothing
    ei, ej = minmax(T.edge, partners[1])
    f = circular_2parallel_apply(g, ei, ej)
    f === nothing && return nothing
    half = (1//2) * alpha(T.colour)
    out = CircularComboR()
    for r in f.newregs
        labels = _circular_transfer_labels(g, f.graph, fdm.region_labels;
                                           overrides = Dict(r => half))
        out = out + CircularComboR(CircularDecorated(f.graph, labels))
    end
    return (combo = out, pair = (ei, ej), partner = partners[1])
end

"""
    circular_rex_fusion_step(fdm::CircularDecoratedMorphism; kwargs...)
        -> Union{Nothing, NamedTuple}

Rewrite `fdm` at the shortest unreduced transition — by plain fusion when the rex
cycle is trivial (`case == :direct`), else by
splicing the derived rex-cycle relation into the cut curve. Returns

* `combo`      — the sum that replaces the diagram (label-free terms),
* `route`      — `:fusion` (trivial cycle) or `:relation`,
* `transition` — the transition of [`circular_unreduced_transition`](@ref),
* `curve`      — the spliced edges `e_1..e_{k+1}` (colours `w·s`); `[T.edge]`
                 on route `:fusion`, where nothing is spliced,
* `partner`    — the other `s`-edge at the transition region if there is one,
                 else `nothing` (reported only — the partner is what the braid
                 move inside the relation creates, file head),
* `sides`      — the assignment anchored on the identity (`nothing` on route
                 `:fusion`),
* `pair`       — the i-pair that carried the fusion: inside the relation on
                 route `:relation`, in the graph itself on route `:fusion`,
* `simple_curve` — the step-3a classification of the curve (reported only, see
                 [`_crf_simple`](@ref)).

`nothing` whenever a gate declines: non-trivial labels, no transition, a `:dot`
transition (that is the dot-slide's case), on route `:fusion` no partner or a
surgery that is not constructible, on route `:relation` a path that is not the
transition path extended by the transition edge, an identity splice that does not
reproduce the graph, or a term that does not splice. Raises only
if the relation itself is not derivable (file head).

`kwargs` are passed to [`circular_rex_relation`](@ref) (`dir`, `recompute`,
`maxpairs`).
"""
function circular_rex_fusion_step(fdm::CircularDecoratedMorphism;
                                  break_budget::Int = 4, kwargs...)
    all(isone, fdm.region_labels) && isone(fdm.outer_label) || return nothing
    m = fdm.m
    g = m.graph

    T = circular_unreduced_transition(m)
    T === nothing && return nothing
    T.kind === :edge || return nothing          # `:dot` is CircularDotSlide's case

    adj, _ = circular_region_adjacency(g)

    # ROUTE `:fusion` — `:direct` = ZERO braid moves, the two `s`-edges already lie
    # at the same region, so plain fusion does it (file head). No relation is
    # derivable here and none is needed.
    if T.case === :direct
        d = _crf_direct_fusion(fdm, T, adj)
        d === nothing && return nothing
        return (combo = d.combo, route = :fusion, transition = T,
                curve = [T.edge], partner = d.partner, sides = nothing,
                pair = d.pair, simple_curve = true)
    end

    # ROUTE `:relation` — the fusion pair (file head): the OTHER s-edge at the
    # SAME region
    partners = unique(Int[ei for (_, ei) in adj[T.region]
                          if ei != T.edge && g.edges[ei].colour == T.colour &&
                             !_circular_edges_touch(g.edges[ei], g.edges[T.edge])])
    # ⚠ REPORTED, NOT GATING: the partner is what the braid move CREATES. At
    # `:dihedral`/`:long` the second `s`-edge lies one or
    # more regions away — `id_1212` has NO partner at `R4`, `id_1213213` none at
    # `R7` — yet the fusion happens INSIDE the derived relation (route
    # `1212 → 2122 → 1212`), where the braid move puts the two `s`-strands next to
    # each other first. Requiring the partner up front assumes the very situation
    # the relation is there to produce, so it is handed back, not demanded.

    p = circular_path_edges(m, T.region)
    p === nothing && return nothing
    q = circular_path_edges(m, T.to)
    # ⚠ THE PATH WORD IS NOT THE REGION WORD. `T.word` comes from the saturated
    # reading and is reduced, `p.word` from the level BFS and may not be — and
    # `circular_rex_relation` then THROWS (`w is not reduced`). The relation lives
    # on the region word; a tree path reading a different, unreduced word is not a
    # curve it can be spliced into. So this is a GATE, not an error: decline, and
    # the other rules carry on.
    ok = q !== nothing &&
         q.edges == vcat(p.edges, T.edge) &&                  # checked, not assumed
         q.word == vcat(p.word, T.colour) &&
         is_reduced(p.word) && !is_reduced(vcat(p.word, T.colour))
    if !ok
        # Where the curve crosses TWO arms of ONE one-coloured node, the two
        # readings differ by the contraction `cc → c` alone. Break that node in two
        # (`_crf_mono_break`) — then the curve crosses only the new connecting
        # edge, path word and saturated word agree, and the gate passes on the
        # retry. A PRE-STEP of this splice, not a rule of its own.
        break_budget > 0 || return nothing
        m2 = _crf_mono_break(m, p)
        m2 === nothing && return nothing
        return circular_rex_fusion_step(
            CircularDecoratedMorphism(m2, fill(one(SoergelPoly), region_count(m2.graph)),
                                      one(SoergelPoly));
            break_budget = break_budget - 1, kwargs...)
    end

    rel = circular_rex_relation(p.word, T.colour; kwargs...)

    # The assignment for the whole relation. The identity pins nothing — it splices
    # under every assignment — so it is taken from the terms themselves.
    n = length(q.word)
    sides = _crf_common_sides(g, q.edges, rel, n)
    sides === nothing && return nothing

    out = CircularComboR()
    for (d, coeff) in pairs_of(rel.rhs)
        tm = CircularMorphismGraph(d.graph, 0, n)             # cuts of the derivation
        r = circular_splice(g, q.edges, tm; sides = sides)
        r === nothing && return nothing
        out = out + coeff * CircularComboR(circular_decorated(r.graph))
    end
    return (combo = out, route = :relation, transition = T, curve = q.edges,
            partner = isempty(partners) ? nothing : partners[1],
            sides = sides, pair = rel.pair,
            simple_curve = _crf_simple(g, q.edges))
end
