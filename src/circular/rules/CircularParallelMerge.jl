# circular/rules/CircularParallelMerge.jl — P1: merge parallel 1/3-edges.
#
# Contains find_circular_parallel_merge, circular_parallel_merge_step,
# _circular_parallel_connection_clean and their small helpers.
#
# ---- P1: merge parallel 1/3-edges -------------------------------------------
#
# WHY: a plain (non-circular) double-braid `13 -> 31 -> 13` is the identity, but a
# circular diagram may instead carry a `1133`-node, which `circular_canonical_key`
# distinguishes from the identity — the same morphism ends up in two different
# normal forms (e.g. on the boundary word `1331`, both planar and
# violation-free).
#
# THE RULE: scan all regions; wherever a 1-edge and a 3-edge sit at the same
# region and do not touch, merge them into ONE `[1,1,3,3]` node. Candidates are
# all edges EXCEPT leaf-to-dot edges. The node is the normal form, not the two
# separate strands.
#
# WHY IT IS NOT IN `CIRCULAR_RULES` — this is forced, not a style choice.
# `reduce_circular` (CircularDriver.jl) asserts after every step that `circular_weight =
# (#braid, #nodes, sum arm_count)` strictly FALLS — that is the termination
# proof. This rule turns two edges into a node, so it INCREASES components 2
# and 3 (the weight goes `(0,0,0) -> (0,1,4)`). It therefore runs — like
# `circular_fusion_step`, excluded from the registry for the same reason — in the
# driver `reduce_to_circular_leave`, at the end of a round and BEFORE the D4 check.
#
# ITS OWN TERMINATION MEASURE. Edge count does not work (it RISES,
# 2 -> 4, because two through-going strands become four stubs at the new
# node). What strictly falls is the number of edge-ENDS not yet attached to
# any node: the two merged edges lose their free ends to the new node, and no
# step creates new free ends. On the `1331` case: 2 leaf-leaf edges
# -> 0, fixed point after ONE step. Since the rule runs outside the
# `circular_weight` driver, this separate measure is fine — it only has to be
# well-founded on its own, and it is.
#
# ARM ORDER: the `[1,1,3,3]` node is unique up to rotation, and `CircularNode`
# equality holds exactly up to rotation (CircularGraph.jl). We fix the four ends in
# slot order (1a, 1b, 3a, 3b) and check the result against
# `euler`/`check_wiring`.

"""
    _circular_edge_is_leaf_to_dot(g::CircularGraph, e::Edge) -> Bool

`true` if `e` connects a leaf to a **dot** (1-armed node). Such edges are
excluded from [`circular_parallel_merge_step`](@ref) (the rule ignores all edges
except leaf-to-dot ones).
"""
function _circular_edge_is_leaf_to_dot(g::CircularGraph, e::Edge)
    is_dot(p) = p isa NodePort && arm_count(g.nodes[p.node]) == 1
    return (e.a isa Leaf && is_dot(e.b)) || (e.b isa Leaf && is_dot(e.a))
end

"""
    _circular_edges_touch(e1::Edge, e2::Edge) -> Bool

`true` if the two edges share an endpoint (same `Leaf`, same `NodePort`, or
the same node). Such pairs are NOT candidates — the rule needs two strands
that do not touch.
"""
function _circular_edges_touch(e1::Edge, e2::Edge)
    node_of(p) = p isa NodePort ? p.node : nothing
    ps1 = (e1.a, e1.b)
    ps2 = (e2.a, e2.b)
    for p in ps1, q in ps2
        p == q && return true
        n1, n2 = node_of(p), node_of(q)
        (n1 !== nothing && n1 == n2) && return true
    end
    return false
end

"""
    _circular_parallel_connection_clean(g, t, dart_edge, ei, ej, region) -> Bool

"Spindle" guard: walks BOTH boundary arcs of the shared `region` between `ei`
and `ej` (the face cycle `t.cycles[fi]` is a single cycle; region == cell ==
face) and requires that NEITHER arc carries a colour-1 or colour-3
edge — only colour-2 edges and plain boundary/leaf stretches (`arc[k]`, no
`dart_edge` entry) are allowed.

Calibrated against the spindle repro: the connection between the open
colour-1 stub at `Leaf(4)` and the next colour-3 edge runs right through the
1/3-arms of the `[1,1,3,3]` node created a step earlier — this guard rejects
that pair, P1 does not fire again on the same stub, and the recursion
terminates.

TERMINATION ARGUMENT (belongs with the file header above): every successful
application of P1 genuinely splits the shared region — the new node wall
separates it in two, because the connection between `ei` and `ej` is NOT made
of 1/3 material (otherwise this function returns `false`). If the connection
stays 1/3-coloured (the spindle), the guard refuses the pair and the recursion
cannot fire again on the same stub.

COLOURS AS A PARAMETER (for 2parallel, `Circular2Parallel.jl`): `forbidden` is the
set of colours that must not occur on the arcs. The default `Set([1, 3])` is
exactly the original, hard-wired P1 reading; 2parallel calls with `Set([s])`
for its one colour.
"""
function _circular_parallel_connection_clean(g::CircularGraph, t, dart_edge::Dict{Int,Int},
                                        ei::Int, ej::Int, region::Int;
                                        forbidden::Set{Int} = Set([1, 3]))
    fi = findfirst(==(region), t.cell_of_face)
    fi === nothing && return false
    cyc = t.cycles[fi]
    n = length(cyc)
    n == 0 && return false
    d_ei = findfirst(d -> get(dart_edge, d, 0) == ei, cyc)
    d_ej = findfirst(d -> get(dart_edge, d, 0) == ej, cyc)
    (d_ei === nothing || d_ej === nothing) && return false
    d_ei == d_ej && return false

    # The two arcs: once forward from `d_ei` to `d_ej`, once from `d_ej` to
    # `d_ei` — BOTH walked cyclically (with wraparound), independent of which
    # of `d_ei`/`d_ej` comes first in the cycle.
    forward(p, q) = (out = Int[]; k = mod1(p + 1, n);
                     while k != q; push!(out, k); k = mod1(k + 1, n); end; out)
    for rng in (forward(d_ei, d_ej), forward(d_ej, d_ei))
        for k in rng
            d = cyc[k]
            eidx = get(dart_edge, d, 0)
            eidx == 0 && continue                        # plain boundary/leaf stretch
            g.edges[eidx].colour in forbidden || continue # different colour — allowed
            return false                                  # forbidden material in the way
        end
    end
    return true
end

"""
    find_circular_parallel_merge(g::CircularGraph) -> Union{Nothing, NamedTuple}

Finds two edges to merge into a `[1,1,3,3]` node: one of colour **1**, one of
colour **3**, both at the **same region**, not touching, neither a
**leaf-to-dot** edge, AND whose connection along both boundary arcs of that
region is different-coloured (`_circular_parallel_connection_clean`, prevents the
"spindle", see there).

Returns `(e1 = edge index colour 1, e3 = edge index colour 3, region = R)` or
`nothing`. Details and derivation in the section header above.
"""
function find_circular_parallel_merge(g::CircularGraph)
    dart_region, t = _circular_dart_region(g)

    # Dart -> edge index, for the boundary walk in
    # `_circular_parallel_connection_clean`. All edges (not just colour 1/3), since
    # colour-2 edges on the arc must be recognised as allowed.
    dart_edge = Dict{Int, Int}()
    for (ei, e) in enumerate(g.edges)
        haskey(t.port_dart, e.a) || continue
        haskey(t.port_dart, e.b) || continue
        dart_edge[t.port_dart[e.a]] = ei
        dart_edge[t.port_dart[e.b]] = ei
    end

    # Regions per edge (both sides), only for candidate edges.
    regs_of = Dict{Int, Vector{Int}}()
    for (ei, e) in enumerate(g.edges)
        (e.colour == 1 || e.colour == 3) || continue
        _circular_edge_is_leaf_to_dot(g, e) && continue
        haskey(t.port_dart, e.a) || continue
        d  = t.port_dart[e.a]
        dr = t.darts[d].rev
        rs = Int[]
        for dd in (d, dr)
            r = get(dart_region, dd, nothing)
            r === nothing || push!(rs, r)
        end
        isempty(rs) || (regs_of[ei] = unique(rs))
    end

    # Deterministic: smallest edge index first, then smallest region.
    for ei in sort(collect(keys(regs_of)))
        g.edges[ei].colour == 1 || continue
        for ej in sort(collect(keys(regs_of)))
            g.edges[ej].colour == 3 || continue
            _circular_edges_touch(g.edges[ei], g.edges[ej]) && continue
            shared = sort(intersect(regs_of[ei], regs_of[ej]))
            for r in shared
                _circular_parallel_connection_clean(g, t, dart_edge, ei, ej, r) || continue
                return (e1 = ei, e3 = ej, region = r)
            end
        end
    end
    return nothing
end

"""
    circular_parallel_merge_step(fdm::CircularDecoratedMorphism) -> Union{Nothing, CircularComboR}

One step of rule **P1**: merges a 1-edge and a 3-edge that sit at the same
region and do not touch into ONE `[1,1,3,3]` node. `nothing` if no pair
matches.

The coefficient is **1** — the rule reshapes the same morphism, it is not a
relation with a prefactor: in the old, non-circular world the double-braid `13 ->
31 -> 13` is the identity.

Region labels are carried over via [`_circular_transfer_labels`](@ref) (when a
region wall disappears in the merge, "the non-trivial label wins"). See the
section header above for termination and where this rule fits — it does NOT
belong in `CIRCULAR_RULES`.
"""
function circular_parallel_merge_step(fdm::CircularDecoratedMorphism)
    g = fdm.m.graph
    m = find_circular_parallel_merge(g)
    m === nothing && return nothing

    e1 = g.edges[m.e1]
    e3 = g.edges[m.e3]

    # SLOT ORDER: the node is `[1,1,3,3]` (same-coloured arms ADJACENT), never
    # `[1,3,1,3]` (alternating) — the two strands run in parallel and their
    # ends must not cross at the node. That fixes everything up to rotation;
    # what remains open is which END of each edge lands on which of the two
    # same-coloured slots — and that decides whether the arms cross.
    #
    # The two admissible cases are therefore (a,b | c,d) and (a,b | d,c): the
    # colour blocks are fixed, only the 3-edge may be flipped. We take
    # whichever comes out planar and violation-free.
    nn = length(g.nodes) + 1

    base = Edge[]
    for (ei, e) in enumerate(g.edges)
        (ei == m.e1 || ei == m.e3) && continue
        push!(base, e)
    end

    # Arms always [1,1,3,3] — same-coloured slots adjacent.
    newnode = _unchecked_circularnode(:mixed, [1, 1, 3, 3])

    g2 = nothing
    for ends in ([(e1.a, 1), (e1.b, 1), (e3.a, 3), (e3.b, 3)],
                 [(e1.a, 1), (e1.b, 1), (e3.b, 3), (e3.a, 3)],
                 [(e1.b, 1), (e1.a, 1), (e3.a, 3), (e3.b, 3)],
                 [(e1.b, 1), (e1.a, 1), (e3.b, 3), (e3.a, 3)])
        keep = copy(base)
        for (slot, (port, col)) in enumerate(ends)
            push!(keep, Edge(col, port, NodePort(nn, slot)))
        end
        gc = CircularGraph(g.word, vcat(g.nodes, [newnode]), keep)
        if euler(gc) == 2 && isempty(check_wiring(gc))
            g2 = gc
            break
        end
    end
    # None of the four connection choices is clean: the rule does not fire
    # here (conservative — better not to merge than to produce a non-planar
    # diagram).
    g2 === nothing && return nothing

    labels2 = _circular_transfer_labels(g, g2, fdm.region_labels)
    # Carry the outer label along.
    return CircularComboR(CircularDecorated(g2, labels2, fdm.outer_label))
end
