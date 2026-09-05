# circular/rules/Circular2ParallelApply.jl — the pure merge surgery for two parallel
# same-colour edges and the CIRCULAR_RULES step that applies it (find + apply +
# label transfer into the two split regions).
#
# Included directly after Circular2Parallel.jl.

"""
    circular_2parallel_apply(g::CircularGraph, ei::Int, ej::Int)
        -> Union{Nothing, NamedTuple}

The pure SURGERY, without any condition — the two same-colour edges `ei`,
`ej` are merged into ONE 4-armed node of that colour. For manual use
(notebooks) and as a building block of [`circular_2parallel_step`](@ref).

Returns `(graph, node, ends, newregs)`; `newregs` are the **two regions the
old shared region splits into** — the sectors after slot 2 and after slot 4,
i.e. the two sectors flanked by arms of DIFFERENT original edges. Slot order
is `[end A of ei, end B of ei, end A of ej, end B of ej]`; since everything is
already the same colour, only the orientation is open, and it is not
guessed: all four ways of connecting are tried and only one is kept that
gives `euler == 2`, a clean `check_wiring`, and a genuine split
(`region_count + 1`, two distinct sectors). If none works, the result is
`nothing` — the rule then does not apply (conservative, like P1).
"""
function circular_2parallel_apply(g::CircularGraph, ei::Int, ej::Int)
    e1, e2 = g.edges[ei], g.edges[ej]
    c = e1.colour
    c == e2.colour || error(
        "circular_2parallel_apply: edge $ei has colour $c, edge $ej has " *
        "colour $(e2.colour) — the rule requires the SAME colour")
    nn   = length(g.nodes) + 1
    base = Edge[e for (k, e) in enumerate(g.edges) if k != ei && k != ej]
    newnode = circular_node([c, c, c, c])

    for ends in ([e1.a, e1.b, e2.a, e2.b], [e1.a, e1.b, e2.b, e2.a],
                 [e1.b, e1.a, e2.a, e2.b], [e1.b, e1.a, e2.b, e2.a])
        keep = copy(base)
        for (slot, port) in enumerate(ends)
            push!(keep, Edge(c, port, NodePort(nn, slot)))
        end
        g2 = CircularGraph(g.word, vcat(g.nodes, [newnode]), keep)
        (euler(g2) == 2 && isempty(check_wiring(g2))) || continue
        region_count(g2) == region_count(g) + 1 || continue
        sec = circular_node_sector_regions(g2)[nn]
        (sec[2] == 0 || sec[4] == 0 || sec[2] == sec[4]) && continue
        return (graph = g2, node = nn, ends = ends, newregs = (sec[2], sec[4]))
    end
    return nothing
end

"""
    circular_2parallel_step(fdm::CircularDecoratedMorphism; strict = true)
        -> Union{Nothing, CircularComboR}

One step of the **2parallel** rule: finds a candidate pair with
[`find_circular_2parallel`](@ref), merges it with [`circular_2parallel_apply`](@ref) and
returns the **two-term sum** — both terms the same graph, coefficient `1`,
with `α_s/2` once in one and once in the other of the two new regions.

`nothing` if no pair matches, the shared region carries a non-trivial label
(condition 5 in the file header), or no way of connecting comes out planar.

The remaining labels go through [`_circular_transfer_labels`](@ref); the two new
regions are set via `overrides` so the transfer leaves them alone.
"""
function circular_2parallel_step(fdm::CircularDecoratedMorphism; strict::Bool = true)
    g = fdm.m.graph
    hit = find_circular_2parallel(fdm.m; strict = strict)
    hit === nothing && return nothing
    # Condition 5: never decided over a region carrying a polynomial.
    isone(fdm.region_labels[hit.region]) || return nothing

    f = circular_2parallel_apply(g, hit.ei, hit.ej)
    f === nothing && return nothing

    half = (1//2) * alpha(hit.colour)
    out  = CircularComboR()
    for r in f.newregs
        labels = _circular_transfer_labels(g, f.graph, fdm.region_labels;
                                      overrides = Dict(r => half))
        out = out + CircularComboR(CircularDecorated(f.graph, labels))
    end
    return out
end
