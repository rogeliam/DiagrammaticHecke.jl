# test/circularbraidongen12.jl — C15 `braid_on_gen12` (src/circular/rules/CircularBraidOnGen12.jl).
#
# The rule: `:gbraid`(8) + `:braid`(6) with three shared edges ↦ `:braid`(6) + two
# trivalents, at the two same-coloured neighbouring pairs of the patch boundary.
#
# The GROUND TRUTH check: the same diagram reduced over a preimage of
# the `:gbraid` (`expand_gbraid`, C13⁻¹) under the base rule set must give the same
# `circular_canonical_key`.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_braid_on_gen12, Edge, NodePort, Leaf, CIRCULAR_RULES

"The bare `n`-armed `:gbraid` (leaves run backwards in slot order)."
function _c15_bare(s::Int, t::Int, n::Int = 8)
    a = [isodd(i) ? s : t for i in 1:n]; lf(i) = i == 1 ? 1 : n + 2 - i
    word = [a[i] for i in 1:n][invperm([lf(i) for i in 1:n])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:n])
end

"""
Attaches an `m`-armed braid-like node to the leaves `k, k+1, k+2`. The planar embedding is not
guessed but searched for (all slot assignments, keeping only the convention-clean ones)
— there are always two, and both give the same
diagram.
"""
function _c15_with_braid(g::CircularGraph, k::Int, m::Int = 6)
    lets = letters(g.word); n = length(lets)
    ks = [mod1(k + i, n) for i in 0:2]
    c1, c2 = lets[ks[1]], lets[ks[2]]
    for start in 1:m, dir in (1, -1), newdir in (1, -1)
        barms = [isodd(i) ? c1 : c2 for i in 1:m]
        sl = [mod1(start + dir * (i - 1), m) for i in 1:3]
        all(barms[sl[i]] == lets[ks[i]] for i in 1:3) || continue
        rest = [s for s in 1:m if !(s in sl)]
        ns = newdir == 1 ? rest : reverse(rest)
        rw = [lets[mod1(k + 2 + i, n)] for i in 1:(n - 3)]
        neu = vcat(rw, [barms[s] for s in ns])
        leafmap = Dict(mod1(k + 2 + i, n) => i for i in 1:(n - 3))
        nodes = vcat(copy(g.nodes), CircularNode[circular_node(barms)]); bi = length(nodes)
        slot_of_old = Dict(ks[i] => sl[i] for i in 1:3)
        rp(p) = p isa Leaf ? (haskey(slot_of_old, p.k) ? NodePort(bi, slot_of_old[p.k]) :
                                                         Leaf(leafmap[p.k])) : p
        edges = Edge[Edge(e.colour, rp(e.a), rp(e.b)) for e in g.edges]
        for (j, s) in enumerate(ns)
            push!(edges, Edge(barms[s], Leaf(length(neu) - (m - 3) + j), NodePort(bi, s)))
        end
        h = CircularGraph(CircularWord(neu), nodes, edges)
        (is_wired(h) && euler(h) == 2 && isempty(check_wiring(h))) && return h
    end
    return nothing
end

@testset "C15 braid_on_gen12 (CircularBraidOnGen12.jl)" begin

    @testset "the equation, at all eight sites and in both colour pairs" begin
        hits = 0
        for (s, t) in ((1, 2), (2, 3)), k in 1:8
            g = _c15_with_braid(_c15_bare(s, t), k)
            res = _fr_braid_on_gen12(g)
            res === nothing && continue
            ts = collect(pairs_of(res))
            length(ts) == 1 && ts[1][2] == 1 || continue
            h = ts[1][1]
            sort([(nd.kind, arm_count(nd)) for nd in h.nodes]) ==
                sort([(:braid, 6), (h.nodes[2].kind, 3), (h.nodes[3].kind, 3)]) || continue
            (count(nd -> nd.kind === :braid, h.nodes) == 1 &&
             count(nd -> arm_count(nd) == 3, h.nodes) == 2) || continue
            (join(letters(h.word)) == join(letters(g.word)) &&
             circular_degree(h) == circular_degree(g) &&
             euler(h) == 2 && isempty(check_wiring(h)) &&
             circular_weight(h) < circular_weight(g)) || continue
            hits += 1
        end
        @test hits == 16
    end

    @testset "the weights" begin
        g = _c15_with_braid(_c15_bare(1, 2), 1)
        h = first(pairs_of(_fr_braid_on_gen12(g)))[1]
        @test circular_arm_weight(g) == (14, 2, 14)
        @test circular_arm_weight(h) == (6, 3, 12)
        @test circular_weight(h) < circular_weight(g)
        @test circular_degree(h) == -2
    end

    @testset "the two trivalents carry THE SAME colour (that of the braid move)" begin
        for k in (1, 2)
            g = _c15_with_braid(_c15_bare(1, 2), k)
            h = first(pairs_of(_fr_braid_on_gen12(g)))[1]
            tri = [nd for nd in h.nodes if arm_count(nd) == 3]
            @test length(tri) == 2
            @test Set(arms(tri[1])) == Set(arms(tri[2]))
            @test length(Set(arms(tri[1]))) == 1        # single-coloured
        end
    end

    # The next two blocks are STATEMENTS ABOUT THE DIRECTION: they require the
    # driver to land on the 3-node form and not on the `:gbraid`. Under
    # `CIRCULAR_TRIVALENT_MERGE = true` C19 rebuilds the `:gbraid`, so what is
    # checked here is the C15 surgery in the world
    # `CIRCULAR_TRIVALENT_MERGE = false`: C14 live, C19 dead, arm weight.
    _c15_direction = CIRCULAR_TRIVALENT_MERGE[]
    CIRCULAR_TRIVALENT_MERGE[] = false
    try

    @testset "GROUND TRUTH: the same as via a preimage without C15" begin
        without_c15 = [r for r in CIRCULAR_RULES if r.name !== :braid_on_gen12]
        g0 = _c15_bare(1, 2)
        k  = 1
        g  = _c15_with_braid(g0, k)
        h  = first(pairs_of(_fr_braid_on_gen12(g)))[1]
        target = nothing
        for v in 1:4
            gv = _c15_with_braid(expand_gbraid(g0, 1, v; rot = 0), k)
            c, _ = reduce_circular(circular_decorated(gv); rules = without_c15)
            for (tt, _) in pairs_of(c)
                any(nd.kind === :braid && arm_count(nd) >= 8 for nd in tt.graph.nodes) && continue
                target = tt.graph
            end
        end
        @test target !== nothing
        @test circular_canonical_key(h) == circular_canonical_key(target)
    end

    @testset "the driver unifies ALL four preimages" begin
        g0 = _c15_bare(1, 2); k = 1
        keys = Any[]
        for v in 1:4
            gv = _c15_with_braid(expand_gbraid(g0, 1, v; rot = 0), k)
            c, _ = reduce_circular(circular_decorated(gv))
            push!(keys, sort([circular_canonical_key(tt.graph) for (tt, _) in pairs_of(c)]))
        end
        @test length(unique(keys)) == 1
        # namely on the 3-node form, not on the :gbraid
        c, _ = reduce_circular(circular_decorated(_c15_with_braid(expand_gbraid(g0, 1, 1; rot = 0), k)))
        @test all(!any(nd.kind === :braid && arm_count(nd) >= 8 for nd in tt.graph.nodes)
                  for (tt, _) in pairs_of(c))
    end

    finally
        CIRCULAR_TRIVALENT_MERGE[] = _c15_direction
    end

    # ---- arm counts above 8 -------------------------------------------------
    #
    # The matcher takes ANY gbraid, not only the 8-armed one: nothing in
    # `_circular_glue3_rewire_at` is tied to a count. One site above the
    # measured one, `(10, 6)`, to keep that honest.
    @testset "C15 at a 10-armed gbraid" begin
        g = nothing
        for k in 1:10
            h = _c15_with_braid(_c15_bare(1, 2, 10), k, 6)
            h === nothing || (g = h; break)
        end
        @test g !== nothing
        @test sort([arm_count(nd) for nd in g.nodes]) == [6, 10]
        res = _fr_braid_on_gen12(g)
        @test res !== nothing
        h = first(pairs_of(res))[1]
        # one braid-like node with nv + nb - 8 arms, plus the two trivalents
        @test sort([arm_count(nd) for nd in h.nodes]) == [3, 3, 8]
        @test euler(h) == 2 && isempty(check_wiring(h)) && is_wired(h)
        @test h.word == g.word
        @test circular_degree(h) == circular_degree(g)
        # braid arm sum 16 -> 8: the drop of exactly 8 the file header argues
        @test circular_arm_weight(h)[1] == circular_arm_weight(g)[1] - 8
        @test circular_weight(h) < circular_weight(g)
    end

    @testset "no false alarm" begin
        @test _fr_braid_on_gen12(_c15_bare(1, 2)) === nothing     # no braid attached
        @test _fr_braid_on_gen12(circular(braid(1, 2; m = 3))) === nothing
    end

    # ---- the check of C6–C9 on the 8-arm case --------------------------------
    #
    # Mixed pair: 8-armed braid + 6-armed braid, three shared edges.
    # Exactly here C7 (`braid_back`, "two braids, three inner edges") or C8
    # (`braid_relation`) could formally apply, since both nodes have `kind === :braid`
    # and `kind` alone does not separate them. RESULT: none of the four touches it;
    # C15 is responsible, enforced by the `arm_count == 6` guard.
    @testset "C6–C9 on the mixed 8/6 pair" begin
        for k in 1:8
            g = _c15_with_braid(_c15_bare(1, 2), k)
            g === nothing && continue
            @test sort([arm_count(nd) for nd in g.nodes]) == [6, 8]
            @test DiagrammaticHecke._fr_dot_into_braid(g)     === nothing   # C6
            @test DiagrammaticHecke._fr_braid_back(g)         === nothing   # C7
            @test DiagrammaticHecke._fr_braid_back_at(g, 1, 2) === nothing  # C7 direct
            @test DiagrammaticHecke._fr_braid_relation(g)     === nothing   # C8
            @test _fr_braid_on_gen12(g) !== nothing                     # C15 it is
        end
    end
end

