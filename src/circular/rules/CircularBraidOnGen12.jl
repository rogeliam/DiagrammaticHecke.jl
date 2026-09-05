# circular/rules/CircularBraidOnGen12.jl — C15 `braid_on_gen12`: a 12-braid that hangs off
# a gbraid node by THREE edges decomposes into an ordinary braid and two
# trivalents.
#
# THE EQUATION
#
#     gbraid (2k >= 8 arms) + `:braid` (6 arms), three shared edges
#         ↦  one braid-like node (2k − 2 arms) + two trivalents
#
# ONE term, coefficient 1. Both sides have the same boundary and `circular_degree = −2`.
#
# WHERE THE TRIVALENTS GO — this is the crux, and it's forced: braiding
# `121 → 212` at one spot of `12121212` gives, for instance, `12122122` — from
# which it is clear where the trivalents must sit.
#
#   * The patch boundary has `2k` connections (2k − 3 from the gbraid, 3
#     from the `:braid`). Their colours, read cyclically, contain EXACTLY TWO
#     adjacent same-coloured pairs — on the bare gbraid the sequence
#     alternates; the braid move `sts ↦ tst` creates the two doublings.
#   * A TRIVALENT of that colour goes at each such pair; it joins the two
#     connections into a single strand.
#   * `2k − 2` strands remain (the single connections plus two trivalent
#     stems), and their colours alternate again — that's the new node. For
#     `2k = 8` it is a `:braid`, above that a gbraid.
#
# So there's nothing to choose: the boundary colours fix the number, colour,
# and location of the trivalents.
#
# The four preimages
# of the gbraid (from `expand_gbraid`) reduce, at each of the eight candidate
# spots, in a 3 : 1 ratio — three land on `:braid` + 2 trivalents, one merges
# via C13 back to gbraid + `:braid` and stays there. Both sides are
# `CIRCULAR_RULES` fixed points. This rule closes exactly that gap; the result
# agrees bit-for-bit with the other three preimages (`circular_canonical_key`).
#
# TERMINATION. The structural identity: the single new node carries
# `nv + nb − 8` arms and the two trivalents are not braid-like. So the BRAID
# ARM SUM — component 1 of `circular_arm_weight` — goes `nv + nb ↦ nv + nb − 8`,
# a drop of exactly 8 at every pair of arm counts, not just at a measured one.
# Under the `sep` and `chain` readings that `CIRCULAR_WEIGHT_MODE = :auto`
# selects today, component 1 falls by 2 instead. Either way `circular_weight`
# falls and the assert in `CircularDriver.jl` holds. Measured at the arm-count
# pairs `(8,6)`, `(10,6)`, `(12,6)`, `(8,8)`, `(10,8)`, `(10,10)` and `(12,8)`:
# the rule fires, the result is planar and wiring-clean, `circular_degree` and
# the boundary word are unchanged, and the weight falls under all three
# readings.
#
# At `(8, 6)` the arm reading is `(14, 2, 14) ↦ (6, 3, 12)`.
#
# ARM COUNTS: the matcher takes any gbraid (`arm_count >= 8`) and any ordinary
# `:braid` (6 arms). Nothing in the surgery is tied to a particular count —
# `_circular_glue3_rewire_at` computes `p = nv + nb − 6` and checks the
# structural conditions itself, so a pair that does not fit is rejected there
# rather than excluded by the matcher.
#
# A gbraid is a `:braid` node with `2k >= 8` arms; there is no node kind of
# its own for it. The names `expand_gbraid`,
# `GBRAID_PREIMAGES`, `gbraid_tower_*` are its public API.

# The three glue edges between the gbraid v and `:braid` b, as
# (gbraid-slot, braid-slot) pairs. `nothing` if there aren't exactly three, or
# they don't connect cyclically on both sides.
function _circular_gen12_braid_glue(g::CircularGraph, v::Int, b::Int)
    glue_pairs = Tuple{Int,Int}[]
    for e in g.edges
        (e.a isa NodePort && e.b isa NodePort) || continue
        if e.a.node == v && e.b.node == b
            push!(glue_pairs, (e.a.slot, e.b.slot))
        elseif e.b.node == v && e.a.node == b
            push!(glue_pairs, (e.b.slot, e.a.slot))
        end
    end
    length(glue_pairs) == 3 || return nothing

    nv = arm_count(g.nodes[v]); nb = arm_count(g.nodes[b])
    # Sort by braid slot so they run 1,2,3 consecutively …
    sort!(glue_pairs; by = p -> p[2])
    b1 = glue_pairs[1][2]
    all(glue_pairs[i][2] == mod1(b1 + i - 1, nb) for i in 1:3) || begin
        # cyclic wraparound (e.g. slots 6,1,2): find the start point
        starts = [s for s in 1:nb if Set(mod1(s + i, nb) for i in 0:2) ==
                                     Set(p[2] for p in glue_pairs)]
        length(starts) == 1 || return nothing
        b1 = starts[1]
        sort!(glue_pairs; by = p -> mod(p[2] - b1, nb))
    end
    # … and then the gbraid slots must connect DESCENDING to match
    j1 = glue_pairs[1][1]
    all(glue_pairs[i][1] == mod1(j1 - (i - 1), nv) for i in 1:3) || return nothing
    return glue_pairs
end

"""
    _circular_glue3_rewire(g::CircularGraph, v::Int, b::Int) -> Union{Nothing, CircularGraph}

The **shared core of C15 and C16**: two braid-like nodes `v` (`nv` arms) and
`b` (`nb` arms), connected by exactly THREE edges, are replaced by ONE
braid-like node with `nv + nb − 8` arms and TWO trivalents.

The computation is the same in both cases and depends only on the boundary
colours:

  * The patch boundary has `p = nv + nb − 6` connections (`nv−3` from one,
    `nb−3` from the other node), read in boundary order.
  * Their colour sequence contains exactly TWO adjacent same-coloured pairs
    (on an untouched braid-like node it alternates; the three glue edges
    create the two doublings). A trivalent of that colour goes at each pair.
  * The remaining `p − 2` strands have colours that alternate again — that's
    the new node: `p − 2 = 6` is a `:braid`, anything larger a gbraid.

`nothing` if the pattern doesn't match (no three cyclically connected glue
edges, not exactly two disjoint same-coloured neighbour pairs, or the
remaining sequence doesn't alternate).
"""
function _circular_glue3_rewire(g::CircularGraph, v::Int, b::Int)
    glue_pairs = _circular_gen12_braid_glue(g, v, b)
    glue_pairs === nothing && return nothing
    return _circular_glue3_rewire_at(g, v, b, glue_pairs)
end

"""
    _circular_glue3_rewire_at(g, v, b, glue_pairs) -> Union{Nothing, CircularGraph}

The surgery of [`_circular_glue3_rewire`](@ref) with the three glue edges GIVEN
rather than searched. C22 (`circular/rules/CircularGen12Null.jl`) needs it: it
runs the same surgery on three of MORE than three shared edges.

The two callers see the same thing whenever there are exactly three shared
edges — "drop the three used ones" and "drop every `v`–`b` edge" are then the
same set. With more shared edges the surplus ones are rewired like any other
outer connection, and that is exactly what produces the needle C22 argues from.
"""
function _circular_glue3_rewire_at(g::CircularGraph, v::Int, b::Int,
                                   glue_pairs::Vector{Tuple{Int,Int}})
    length(glue_pairs) == 3 || return nothing
    ndv, ndb = g.nodes[v], g.nodes[b]
    nv, nb = arm_count(ndv), arm_count(ndb)

    p = nv + nb - 6                     # connections on the patch boundary
    p >= 8 || return nothing
    j3 = glue_pairs[3][1]; b1 = glue_pairs[1][2]

    # The patch's outer connections, in BOUNDARY ORDER: first the free slots
    # of `v` (descending from `j3 − 1`), then the free slots of `b`
    # (descending from `b1 − 1`). This order is the house convention "leaves
    # run backwards in slot order" (see CircularGen12Expand.jl, same point).
    ports = Tuple{Int,Int}[]                      # (node, slot)
    for i in 1:(nv - 3); push!(ports, (v, mod1(j3 - i, nv))); end
    for i in 1:(nb - 3); push!(ports, (b, mod1(b1 - i, nb))); end

    # the colour of each connection = colour of the edge hanging there
    outer = Dict{Tuple{Int,Int},Any}()
    for e in g.edges
        for (q1, q2) in ((e.a, e.b), (e.b, e.a))
            q1 isa NodePort || continue
            (q1.node, q1.slot) in ports || continue
            outer[(q1.node, q1.slot)] = (e.colour, q2)
        end
    end
    length(outer) == p || return nothing
    colours = [outer[q][1] for q in ports]

    # exactly two adjacent same-coloured pairs, disjoint, same colour
    pos = [i for i in 1:p if colours[i] == colours[mod1(i + 1, p)]]
    length(pos) == 2 || return nothing
    (mod1(pos[1] + 1, p) != pos[2] && mod1(pos[2] + 1, p) != pos[1]) || return nothing
    cp = colours[pos[1]]
    colours[pos[2]] == cp || return nothing

    # Rotate the boundary order so the first pair sits at position 1 —
    # otherwise one of the pairs would run across the seam (position `p` and
    # 1) and the strand loop below would count it twice. Rotating is harmless:
    # the order is cyclic.
    rot    = pos[1] - 1
    ports  = circshift(ports, -rot)
    colours = circshift(colours, -rot)
    pos    = [1, pos[2] - pos[1] + 1]

    # --- build the new wiring -----------------------------------------
    others = [u for u in 1:length(g.nodes) if u != v && u != b]
    idx   = Dict(u => t for (t, u) in enumerate(others))
    basis = length(others)
    NB = basis + 1                      # the new braid-like node
    T1 = basis + 2; T2 = basis + 3      # the two trivalents

    # Strands in boundary order: single connection or trivalent stem
    strands = Any[]; i = 1
    pair_of = Dict(pos[1] => T1, pos[2] => T2)
    while i <= p
        if haskey(pair_of, i)
            push!(strands, (:triv, pair_of[i], i, mod1(i + 1, p)))
            i += 2
        else
            push!(strands, (:port, ports[i], colours[i]))
            i += 1
        end
    end
    n_new = p - 2
    length(strands) == n_new || return nothing

    strand_colour(s) = s[1] === :triv ? cp : s[3]
    # Arms in slot order: strand `t` gets slot `n_new − t + 1`
    # (descending), so the boundary order runs backwards again.
    slot_of_strand = [mod1(n_new - t + 1, n_new) for t in 1:n_new]
    new_arms = Vector{Int}(undef, n_new)
    for t in 1:n_new; new_arms[slot_of_strand[t]] = strand_colour(strands[t]); end
    # alternating? (must hold, otherwise no braid shape fits here)
    all(new_arms[t] != new_arms[mod1(t + 1, n_new)] for t in 1:n_new) || return nothing

    nodes = vcat(CircularNode[g.nodes[u] for u in others],
                 CircularNode[circular_node(new_arms), circular_node([cp, cp, cp]),
                         circular_node([cp, cp, cp])])

    # target port for each old outer connection
    target = Dict{Tuple{Int,Int},NodePort}()
    for (t, s) in enumerate(strands)
        if s[1] === :port
            target[s[2]] = NodePort(NB, slot_of_strand[t])
        else
            (_, tn, i1, i2) = s
            target[ports[i1]] = NodePort(tn, 1)
            target[ports[i2]] = NodePort(tn, 3)
        end
    end

    shift(q) = q isa Leaf ? q :
               (haskey(target, (q.node, q.slot)) ? target[(q.node, q.slot)] :
                NodePort(idx[q.node], q.slot))

    edges = Edge[]
    used = Set(glue_pairs)
    for e in g.edges
        # the three USED glue edges drop out; a surplus `v`–`b` edge (C22)
        # survives and is rewired through `shift` like any other connection.
        if e.a isa NodePort && e.b isa NodePort &&
           (e.a.node in (v, b)) && (e.b.node in (v, b))
            pr = e.a.node == v ? (e.a.slot, e.b.slot) : (e.b.slot, e.a.slot)
            pr in used && continue
        end
        push!(edges, Edge(e.colour, shift(e.a), shift(e.b)))
    end
    # the two trivalent stems into the new node
    for (t, s) in enumerate(strands)
        s[1] === :triv || continue
        push!(edges, Edge(cp, NodePort(s[2], 2), NodePort(NB, slot_of_strand[t])))
    end

    return CircularGraph(g.word, nodes, edges)
end

"""
    _fr_braid_on_gen12(g::CircularGraph) -> Union{Nothing, CircularComboR}

**C15 `braid_on_gen12`.** A `:braid` (6 arms) that hangs off a gbraid
(`arm_count >= 8`) by **three** edges is replaced, together with it, by **one**
braid-like node and **two** trivalents — at the patch boundary's two
same-coloured neighbour pairs. One term, coefficient 1.

The surgery is the shared core [`_circular_glue3_rewire`](@ref); this function only
supplies the pattern (which two node kinds, which arm counts).

For the derivation and why the trivalent positions are forced: see the
file header.
"""
function _fr_braid_on_gen12(g::CircularGraph)
    for (v, ndv) in enumerate(g.nodes)
        (ndv.kind === :braid && arm_count(ndv) >= 8) || continue
        for (b, ndb) in enumerate(g.nodes)
            (ndb.kind === :braid && arm_count(ndb) == 6) || continue
            h = _circular_glue3_rewire(g, v, b)
            h === nothing && continue
            return CircularComboR(h)
        end
    end
    return nothing
end
