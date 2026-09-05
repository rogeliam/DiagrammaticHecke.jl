# circular/rules/CircularRuleFixtures.jl — named example diagrams for the rules
# that no double leaf of a short word produces.
#
# WHY THIS FILE EXISTS. Most rules can be shown by running the driver on double
# leaves and catching the first diagram each rule matches. Seven cannot: their
# pattern needs a many-armed node, a specific double connection, or an unfolded
# preimage, and none of those turns up in the small pairs. The builders here
# produce one such diagram per rule, so the rule inventory in `examples/` and
# the tests can name a witness in one line instead of carrying their own copy.
#
# Every builder returns a diagram that is WIRED, planar (`euler == 2`) and free
# of `check_wiring` violations, or `nothing` if the requested shape does not
# embed. Nothing here is used by a rule; this is example data.

"""
    circular_braid_star(s, t, n) -> CircularGraph

A single braid-like node with `n` alternating `s`/`t` arms, every arm on its own
boundary leaf. The leaves run BACKWARDS in slot order, the house convention.

The seed for [`circular_glued_braids`](@ref): `n = 6` is an ordinary braid,
`n >= 8` a gbraid.
"""
function circular_braid_star(s::Int, t::Int, n::Int)
    a = [isodd(i) ? s : t for i in 1:n]
    lf(i) = i == 1 ? 1 : n + 2 - i
    word = [a[i] for i in 1:n][invperm([lf(i) for i in 1:n])]
    return CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
                         Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:n])
end

"""
    circular_glued_braids(nv, nb, nglue; s = 1, t = 2) -> Union{Nothing, CircularGraph}

Two braid-like nodes, `nv` and `nb` arms, joined by `nglue` edges at
consecutive slots on both sides; every remaining arm goes to a boundary leaf.

The planar embedding is SEARCHED, not guessed: all slot assignments are tried
and only a wired, planar, violation-free one is returned. `nothing` if none
exists.

The glue count picks the rule:

| `nglue` | `nb` | rule |
|---|---|---|
| 3 | 6 | C15 `braid_on_gen12` |
| 3 | >= 8 | C16 `gen12_on_gen12` |
| >= 4 | any | C22 `gen12_on_gen12_null` (the diagram is 0) |

Two edges match no rule directly — C18 wants one of the two nodes unfolded
first, see [`circular_two_edge_fixture`](@ref).
"""
function circular_glued_braids(nv::Int, nb::Int, nglue::Int; s::Int = 1, t::Int = 2)
    g0 = circular_braid_star(s, t, nv)
    lets = letters(g0.word)
    n = length(lets)
    for k in 1:n
        ks = [mod1(k + i, n) for i in 0:(nglue - 1)]
        lets[ks[1]] == lets[ks[2]] && continue
        for start in 1:nb, dir in (1, -1), newdir in (1, -1)
            barms = [isodd(i) ? lets[ks[1]] : lets[ks[2]] for i in 1:nb]
            sl = [mod1(start + dir * (i - 1), nb) for i in 1:nglue]
            all(barms[sl[i]] == lets[ks[i]] for i in 1:nglue) || continue
            rest = [q for q in 1:nb if !(q in sl)]
            ns = newdir == 1 ? rest : reverse(rest)
            rw = [lets[mod1(k + nglue - 1 + i, n)] for i in 1:(n - nglue)]
            word = vcat(rw, [barms[q] for q in ns])
            leafmap = Dict(mod1(k + nglue - 1 + i, n) => i for i in 1:(n - nglue))
            nodes = vcat(copy(g0.nodes), CircularNode[circular_node(barms)])
            bi = length(nodes)
            slot_of_old = Dict(ks[i] => sl[i] for i in 1:nglue)
            rp(p) = p isa Leaf ?
                    (haskey(slot_of_old, p.k) ? NodePort(bi, slot_of_old[p.k]) :
                                                Leaf(leafmap[p.k])) : p
            edges = Edge[Edge(e.colour, rp(e.a), rp(e.b)) for e in g0.edges]
            for (q, sq) in enumerate(ns)
                push!(edges, Edge(barms[sq], Leaf(length(word) - (nb - nglue) + q),
                                  NodePort(bi, sq)))
            end
            h = CircularGraph(CircularWord(word), nodes, edges)
            (is_wired(h) && euler(h) == 2 && isempty(check_wiring(h))) && return h
        end
    end
    return nothing
end

"""
    circular_two_edge_fixture(; s = 1, t = 2) -> Union{Nothing, CircularGraph}

The C18 `gen12_two_edges` pattern: take the two glued gbraids of
[`circular_glued_braids`](@ref) at `(8, 8, 3)` and unfold ONE of them
([`expand_gbraid`](@ref), C13 backwards). What is left is a gbraid and a braid
joined by two edges — the shape C18 matches. `nothing` if no rotation embeds.
"""
function circular_two_edge_fixture(; s::Int = 1, t::Int = 2)
    g = circular_glued_braids(8, 8, 3; s = s, t = t)
    g === nothing && return nothing
    for v in 1:2, r in 0:7
        h = expand_gbraid(g, v, 1; rot = r)
        h === nothing && continue
        (is_wired(h) && euler(h) == 2 && isempty(check_wiring(h))) && return h
    end
    return nothing
end

"""
    circular_bead_graph(c = 1) -> CircularGraph

A single 2-armed node of colour `c` with both arms on leaves — the "bead" C4
removes. `circular_node` rejects an arm count of 2, so the node is built with
the inner constructor; a bead is a transient state, not a valid end state
(`CircularGraph.jl`).
"""
circular_bead_graph(c::Int = 1) =
    CircularGraph(CircularWord([c, c]), CircularNode[CircularNode(:mixed, [c, c])],
                  Edge[Edge(c, Leaf(1), NodePort(1, 1)),
                       Edge(c, Leaf(2), NodePort(1, 2))])

"""
    circular_counter_joined(armsA, armsB) -> Union{Nothing, CircularGraph}

Two nodes joined by exactly TWO edges, at slots `(1, 2)` on the first and
`(2, 1)` on the second — cyclically adjacent on both and COUNTER-oriented. The
remaining arms go to leaves, in whichever of the four traversals embeds
planarly; `nothing` if none does.

That is the site C10 `commutation_merge` and C12 `two_adjacent_merge` share; C10
takes it when the two connecting edges carry the `{1,3}` pair, as in
`circular_counter_joined([1,3,1,3,1,3], [3,1,3,1,3,1])`.
"""
function circular_counter_joined(armsA::Vector{Int}, armsB::Vector{Int})
    A = circular_node(armsA)
    B = circular_node(armsB)
    # Which way round the free slots meet the boundary is not free to guess:
    # only one of the four traversals embeds without crossings. Try them and
    # keep the planar one, as `circular_glued_braids` does.
    for da in (1, -1), db in (1, -1)
        sa = da == 1 ? (3:length(armsA)) : reverse(3:length(armsA))
        sb = db == 1 ? (3:length(armsB)) : reverse(3:length(armsB))
        edges = Edge[Edge(1, NodePort(1, 1), NodePort(2, 2)),
                     Edge(3, NodePort(1, 2), NodePort(2, 1))]
        word = Int[]
        for s in sa
            push!(word, arm_colour(A, s))
            push!(edges, Edge(arm_colour(A, s), Leaf(length(word)), NodePort(1, s)))
        end
        for s in sb
            push!(word, arm_colour(B, s))
            push!(edges, Edge(arm_colour(B, s), Leaf(length(word)), NodePort(2, s)))
        end
        g = CircularGraph(CircularWord(word), [A, B], edges)
        (is_wired(g) && euler(g) == 2 && isempty(check_wiring(g))) && return g
    end
    return nothing
end

"""
    circular_mixed_at_braid(k, colours, armword) -> CircularGraph

A `2k`-armed braid node of `colours` with a mixed `{1,3}` node hanging on its
slot 1; `armword` is the arm sequence of that node and has to start on the braid
colour. Every other arm goes to a leaf.

The C23 `mixed_into_gen12` site: `circular_mixed_at_braid(3, (1, 2), [1, 3, 1, 1, 3])`
is a `121212` braid carrying a `13113` node, which C23 breaks open and merges.
"""
function circular_mixed_at_braid(k::Int, colours::Tuple{Int,Int}, armword::Vector{Int})
    n = 2k
    a = [isodd(s) ? colours[1] : colours[2] for s in 1:n]
    m = length(armword)
    armword[1] == a[1] ||
        throw(ArgumentError("circular_mixed_at_braid: armword must start on the braid " *
                            "colour $(a[1]), got $(armword[1])"))
    nodes = CircularNode[CircularNode(:braid, a)]
    edges = Edge[]
    word = Int[]
    for j in 1:n
        s = mod1(2 - j, n)
        if s == 1
            push!(nodes, circular_node(armword))
            v = length(nodes)
            push!(edges, Edge(a[1], NodePort(1, 1), NodePort(v, 1)))
            for t in m:-1:2
                push!(word, armword[t])
                push!(edges, Edge(armword[t], Leaf(length(word)), NodePort(v, t)))
            end
        else
            push!(word, a[s])
            push!(edges, Edge(a[s], Leaf(length(word)), NodePort(1, s)))
        end
    end
    return CircularGraph(CircularWord(word), nodes, edges)
end
