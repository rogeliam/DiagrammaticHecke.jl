# circular/rules/CircularDotOnGen12.jl — C14 `dot_on_gen12`: a dot on one arm of a
# generalized braid node (a `:braid` node with `2k >= 8` arms).
#
# THE RULE. Let `nd` be a braid-like node with `n = 2k` arms `a₁…aₙ`
# (alternating colours `s,t,s,t,…`), and let a dot of colour `c = a_j` sit on
# arm `j`. Then:
#
#   * the three arms `a_{j-1}, a_j, a_{j+1}` (colours `c', c, c'`) disappear,
#   * one new arm of colour `c'` takes their place — the node now has `n − 2`
#     arms and still alternates (`…, c, c', c, …`); for `n = 8` this is an
#     ordinary `:braid`,
#   * that new arm connects to a trivalent node of colour `c'`, whose other two
#     arms take over the edges of `a_{j-1}` and `a_{j+1}` — i.e. exactly where
#     the dot sat.
#
# One term, coefficient 1. The boundary word loses exactly the letter `c` (the
# dot is an internal node, not a leaf — it disappears together with its edge).
#
# The rule gives the same `circular_canonical_key` as the full reduction of the
# preimage with the dot (dot first, then C1–C13). C6 `dot_into_braid` does not
# cover this spot: it matches only `kind === :braid` with six arms.
#
# THE ORIENTATION IS NOT FREE. The trivalent's three arms must sit
# cyclically as `(connector, old arm j−1, old arm j+1)`. The other order gives
# the same node kinds but `euler = 0` and two `check_wiring` violations —
# correct shape, wrong embedding. `check_wiring` and `circular_canonical_key`,
# not just node kinds, pin the result down.
#
# TERMINATION. This is why component 1 of `circular_weight` is the ARM SUM over
# braid-like nodes and not a count of nodes with exactly 6 arms: the input has
# `2k >= 8` arms and would not count under such a scheme, yet the `2k − 2`-armed
# result may have exactly 6 arms and would — a node count could rise from 0 to 1
# at unchanged total node count. The arm sum falls strictly, `8 → 6`. See the
# CircularWeight.jl header.
#
# The rule is written for general `n` (`2k ↦ 2k−2`); `k ≥ 5` is a CONJECTURE.
# A 10-armed braid-like node first arises from the chain `121212 → 212121`.
#
_circular_is_dot_node(nd::CircularNode) = arm_count(nd) == 1

# Looks for a dot node hanging off an arm of the gbraid-like node `i`.
# Returns `(slot, dot-node-index)` or `nothing`.
function _circular_dot_at_gbraid(g::CircularGraph, i::Int)
    for e in g.edges
        if e.a isa NodePort && e.a.node == i && e.b isa NodePort && _circular_is_dot_node(g.nodes[e.b.node])
            return (e.a.slot, e.b.node)
        elseif e.b isa NodePort && e.b.node == i && e.a isa NodePort && _circular_is_dot_node(g.nodes[e.a.node])
            return (e.b.slot, e.a.node)
        end
    end
    return nothing
end

"""
    _fr_dot_on_gen12(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C14 `dot_on_gen12`.** ⚠ Gated by the dot policy: with
`policy = true` (the default, [`CIRCULAR_DOT_POLICY`](@ref)) only nodes passing
[`circular_dots_reducible`](@ref) are touched — two adjacent dots, or a colour
left with fewer than three arms. `policy = false` measures the ungated behaviour.

A dot on one arm of a braid-like node with `2k` arms
caps three consecutive arms and replaces them with one arm to a newly created
**trivalent node of the other colour**; what remains is a braid-like node with
`2k − 2` arms (for `k = 4`: an ordinary `:braid`). One term, coefficient 1.

For the measurement, the orientation question, and termination see the file
header.
"""
function _fr_dot_on_gen12(g::CircularGraph; policy::Bool = CIRCULAR_DOT_POLICY[])
    for (i, nd) in enumerate(g.nodes)
        # THE DOT POLICY (`circular_dots_reducible`): the dots
        # of a `2k >= 8` node may only be pushed through when two of them are
        # ADJACENT or when a colour drops below three arms. Otherwise the dot
        # stays and the diagram counts as reduced. `policy = false` gets the
        # unconditional, ungated behaviour back for measurements.
        policy && !circular_dots_reducible(g, i) && continue
        # C14 wants only `:braid` nodes with `2k ≥ 8` arms (the 6-armed braid
        # belongs to C6).
        (nd.kind === :braid && arm_count(nd) >= 8) || continue
        hit = _circular_dot_at_gbraid(g, i)
        hit === nothing && continue
        (j, d) = hit

        n = arm_count(nd); a = collect(arms(nd))
        jm = mod1(j - 1, n); jp = mod1(j + 1, n)
        cp = a[jm]                                     # the OTHER colour
        a[jp] == cp || continue                        # cannot fail on a gbraid

        # The new node: the `n-3` untouched arms (from `j+2` cyclically to
        # `j-2`), followed by the connector to the trivalent node.
        keep    = [mod1(j + 1 + t, n) for t in 1:(n - 3)]
        newarms = vcat([a[s] for s in keep], [cp])
        conn    = length(newarms)
        newnode = circular_node(newarms)
        triv    = circular_node([cp, cp, cp])

        old   = [v for v in 1:length(g.nodes) if v != i && v != d]
        idx   = Dict(v => t for (t, v) in enumerate(old))
        nodes = vcat(CircularNode[g.nodes[v] for v in old], CircularNode[newnode, triv])
        NI = length(old) + 1; TI = length(old) + 2

        slotmap = Dict(keep[t] => t for t in 1:(n - 3))
        # Cyclic order at the trivalent node: 1 = connector, 2 = old arm
        # `j-1`, 3 = old arm `j+1`. NOT interchangeable (see file header).
        tslot = Dict(jm => 2, jp => 3)

        function shift(p)
            p isa Leaf && return p
            p.node == i || return NodePort(idx[p.node], p.slot)
            haskey(slotmap, p.slot) && return NodePort(NI, slotmap[p.slot])
            haskey(tslot, p.slot)   && return NodePort(TI, tslot[p.slot])
            return nothing                              # the dot's own arm drops out
        end

        edges = Edge[]
        for ed in g.edges
            (ed.a isa NodePort && ed.a.node == d) && continue
            (ed.b isa NodePort && ed.b.node == d) && continue
            pa = shift(ed.a); pb = shift(ed.b)
            (pa === nothing || pb === nothing) && continue
            push!(edges, Edge(ed.colour, pa, pb))
        end
        push!(edges, Edge(cp, NodePort(NI, conn), NodePort(TI, 1)))

        return CircularComboR(CircularGraph(g.word, nodes, edges))
    end
    return nothing
end

"""
    _fr_dot_on_gen12_collapse(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C14 as a COMPOSITE**, and this is the form the registry uses. A single C14
step RISES by exactly one separating edge; the descent only arrives once the
policy-degenerate node has collapsed completely. So the registered rule does
both in one go, after the C18 pattern: a policy-gated C14 step
(`_fr_dot_on_gen12`), then `reduce_circular` WITHOUT `:dot_on_gen12`, accepted
only if ALL terms are strictly lighter than the input and carry no foreign
labels. That keeps the driver's termination assert
(`_assert_circular_weight_drops`) intact.

A C14 ↔ C19 circle is excluded twice over: C19 refuses to BUILD a node the dot
policy would immediately collapse (the loop guard in
`_circular_mono_into_gen12_at`), and the intermediate reduction here runs
without C19.

C14 itself runs as a CHAIN in the work list: at `2k ≥ 10` arms with several
dots the node is still policy-degenerate after one step and needs C14 again.
The chain terminates because every C14 step lowers the node's arm count by two.

`_fr_dot_on_gen12` stays the single-step ground truth the tests are written
against.
"""
function _fr_dot_on_gen12_collapse(g::CircularGraph)
    step = _fr_dot_on_gen12(g)
    step === nothing && return nothing
    w0 = circular_weight(g)
    other_rules = [r for r in CIRCULAR_RULES
                   if r.name !== :dot_on_gen12 && r.name !== :trivalent_into_gen12]
    out  = CircularComboR()
    work = Vector{Tuple{CircularGraph,Any}}()
    for (t0, c0) in pairs_of(step)
        push!(work, (t0 isa CircularGraph ? t0 : t0.graph, c0))
    end
    while !isempty(work)
        (h, coeff) = pop!(work)
        c, _ = reduce_circular(circular_decorated(h); rules = other_rules)
        isempty(pairs_of(c)) && return nothing
        for (t, ct) in pairs_of(c)
            all(isone, t.region_labels) || return nothing
            nxt = _fr_dot_on_gen12(t.graph)
            if nxt === nothing
                # collapsed — the guard: strictly lighter than the input
                circular_weight(t.graph) < w0 || return nothing
                _add!(out, t.graph, coeff * ct)
            else
                for (u, cu) in pairs_of(nxt)
                    push!(work, (u isa CircularGraph ? u : u.graph, coeff * ct * cu))
                end
            end
        end
    end
    return out
end
