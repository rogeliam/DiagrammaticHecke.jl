# diagram/JoinCircularMirror.jl — hflip/vflip for CircularMorphismGraph: the
# circular-node mirror (_circular_mirror_node, _circular_mirror_slot,
# _circular_mirror_graph) and the horizon/vertical flip pair. Included after
# JoinShift.jl.
#
# ---- hflip / vflip for CircularMorphismGraph ---------------------------------------
#
# NAMING — DIFFERENT from the plain stack. A mirror is named after its
# AXIS:
#   `hflip` — axis = the HORIZON, the line through the two MARKERS (cuts). Swaps
#             bottom and top.
#   `vflip` — axis = the VERTICAL, perpendicular to it. Bottom stays bottom, top
#             stays top; only the letter order reverses (left↔right).
# In the plain `MorphismGraph` stack these two names are SWAPPED (`flip` = horizon,
# `hflip` = vertical). Those names stay as they are — `flip` is depended on by
# `double_leaf`, `LightLeaves.jl`, `Zamolodchikov.jl` and many tests.
# `flip(::CircularMorphismGraph)` is an alias of `hflip`.
#
# A SEPARATE IMPLEMENTATION is needed for `CircularMorphismGraph`:
# `light_leaf`/`double_leaf`/`compose`/`flip` all run on `MorphismGraph`, BEFORE any
# circular conversion, and `circular_morphism` comes after. Mirroring a FINISHED
# circular morphism therefore needs the operation directly on `CircularMorphismGraph`.
#
# THE DIFFERENCE TO `_mirror_graph`. A `Node` carries its slot colours implicitly
# (braid: slot parity; trivalent/dot: the one node colour), a `CircularNode` carries them
# EXPLICITLY in `arms`. Mirroring a circular node is therefore not just a slot
# permutation but the REVERSAL OF THE ARM SEQUENCE — and the two must agree, or a
# slot points at the wrong colour.
#
# Rule: `slot ↦ mod1(d − slot, d)` (slot `d` fixed, the rest reversed;
# `d = arm_count`) together with `newarms[s] = arms[mod1(d − s, d)]`. Then
# `arm_colour(newnd, mirror(s)) == arm_colour(nd, s)` for EVERY slot — colour
# fidelity is constructive here, not tied to a parity computation as for `:braid`. A
# 1-armed circular dot (`d = 1`) stays fixed automatically.
#
# The boundary is mirrored literally as in `_mirror_graph` (`Leaf(k) ↦
# Leaf(mod1(axis − k, n))`, `axis = cut1 + cut2 + 1`, `_EMPTY_BOTTOM_CUT` normalised
# to `0`) — same axis, same sentinel handling, so the contracts
# `bottom(flip(m)) == top(m)` / `top(flip(m)) == bottom(m)` hold as before.

"""
    _circular_mirror_node(nd::CircularNode) -> CircularNode

Mirrors a single `CircularNode`: the cyclic arm order is REVERSED, slot `d`
(`d = arm_count(nd)`) stays fixed. Counterpart to the slot rule
`_circular_mirror_slot`; together they are colour-faithful (see the section header).

Builds via `_unchecked_circularnode` with the UNCHANGED `kind`: reversing an arm sequence
cannot change the node kind (`:mono` stays monochrome 2, `:mixed` stays in {1,3} with
the same multiplicities, `:braid` stays alternating with the same colour pair), and
the regular constructor would needlessly reject the transient 2-armed intermediate
nodes of the rules.
"""
function _circular_mirror_node(nd::CircularNode)
    d = arm_count(nd)
    return _unchecked_circularnode(nd.kind, [nd.arms[mod1(d - s, d)] for s in 1:d])
end

"The slot rule matching the node mirror `_circular_mirror_node`: `slot ↦ mod1(d − slot, d)`."
_circular_mirror_slot(slot::Int, d::Int) = mod1(d - slot, d)

"""
    _circular_mirror_graph(m::CircularMorphismGraph) -> CircularGraph

Geometric core of `flip(::CircularMorphismGraph)`/`hflip(::CircularMorphismGraph)` — the
literal analogue of [`_mirror_graph`](@ref), with the `CircularNode` mirror
(`_circular_mirror_node`/`_circular_mirror_slot`) instead of the braid/trivalent slot rules.
Returns the graph unchanged if the boundary is empty.
"""
function _circular_mirror_graph(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    n == 0 && return g
    c1 = m.cut1 == _EMPTY_BOTTOM_CUT ? 0 : m.cut1
    c2 = m.cut2 == _EMPTY_BOTTOM_CUT ? 0 : m.cut2
    axis = c1 + c2 + 1
    newnodes = [_circular_mirror_node(nd) for nd in g.nodes]
    function remap(p::Port)
        if p isa Leaf
            return Leaf(mod1(axis - p.k, n))
        elseif p isa NodePort
            return NodePort(p.node, _circular_mirror_slot(p.slot, arm_count(g.nodes[p.node])))
        else
            return p                                # Circle
        end
    end
    newedges = [Edge(e.colour, remap(e.a), remap(e.b)) for e in g.edges]
    oldcols  = [leaf_colour(g, k) for k in 1:n]
    newword  = CircularWord([oldcols[mod1(axis - i, n)] for i in 1:n])
    return CircularGraph(newword, newnodes, newedges)
end

"""
    hflip(m::CircularMorphismGraph) -> CircularMorphismGraph
    flip(m::CircularMorphismGraph)  -> CircularMorphismGraph   (alias, see below)

Mirror in the **HORIZON** — the line through the two MARKERS (cuts). It swaps bottom
and top. For `m : a → b` we get `hflip(m) : b → a`, with the LETTERS of each boundary
word in their order —

    bottom(hflip(m)) == top(m)
    top(hflip(m))    == bottom(m)

An involution. Counterpart: [`vflip`](@ref) (mirror in the vertical).

**NAMING.** The axis of this mirror is the horizon, hence `hflip` — unlike the plain
`MorphismGraph` stack, where the same map is called `flip` and `hflip` denotes the
vertical mirror. The name `flip(::CircularMorphismGraph)` is kept as an alias;
NEW code should use `hflip`/`vflip`.

Technically: the same boundary mirror as `vflip`, the same special handling of the
`cut1 == cut2` case (sentinel toggle instead of a no-op swap), with the circular node
mirror `_circular_mirror_node` — derivation in the section header above.
"""
function hflip(m::CircularMorphismGraph)
    n = length(m.graph.word)
    n == 0 && return CircularMorphismGraph(m.graph, m.cut2, m.cut1)
    newgraph = _circular_mirror_graph(m)
    if m.cut1 == m.cut2
        newcut = m.cut1 == _EMPTY_BOTTOM_CUT ? 0 : _EMPTY_BOTTOM_CUT
        return CircularMorphismGraph(newgraph, newcut, newcut)
    end
    return CircularMorphismGraph(newgraph, m.cut2, m.cut1)
end

# Alias for `hflip` on a `CircularMorphismGraph`; new code uses `hflip`/`vflip`.
flip(m::CircularMorphismGraph) = hflip(m)

"""
    vflip(m::CircularMorphismGraph) -> CircularMorphismGraph

Mirror in the **VERTICAL** — the axis PERPENDICULAR to the horizon. Bottom stays
bottom and top stays top (the roles of the cuts do not change), only the LETTER
ORDER of each boundary word reverses —

    bottom(vflip(m)) == reverse(bottom(m))
    top(vflip(m))    == reverse(top(m))

An involution. Counterpart: [`hflip`](@ref) (mirror in the horizon, swaps bottom and
top).

**NAMING:** left↔right is mirrored in the vertical, hence `vflip`. In the plain
`MorphismGraph` stack this very map is called `hflip` — there the two names are
swapped relative to this convention (see [`hflip`](@ref)`(::CircularMorphismGraph)`).

Technically the same boundary mirror as `hflip` (`_circular_mirror_graph`), except the
cuts keep their roles: no swap, no sentinel toggle.
"""
function vflip(m::CircularMorphismGraph)
    n = length(m.graph.word)
    n == 0 && return CircularMorphismGraph(m.graph, m.cut1, m.cut2)
    return CircularMorphismGraph(_circular_mirror_graph(m), m.cut1, m.cut2)
end
