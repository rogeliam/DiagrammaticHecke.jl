# test/circulargen12ongen12.jl — C16 `gen12_on_gen12` (src/circular/rules/CircularGen12OnGen12.jl).
#
# The rule: `:gbraid`(8) + `:gbraid`(8) with three shared edges ↦ `:gbraid`(8) + two
# trivalents, at the two same-coloured neighbouring pairs of the
# patch boundary.
#
# As for C15 the central test is the GROUND TRUTH: the same diagram reduced over
# the preimages of BOTH `:gbraid` (`expand_gbraid`, C13⁻¹) under the base rule set must
# give the same `circular_canonical_key`.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_gen12_on_gen12, Edge, NodePort, Leaf, CIRCULAR_RULES

"The bare 8-armed `:gbraid` (leaves run backwards in slot order)."
function _c16_bare(s::Int, t::Int, n::Int = 8)
    a = [isodd(i) ? s : t for i in 1:n]; lf(i) = i == 1 ? 1 : n + 2 - i
    word = [a[i] for i in 1:n][invperm([lf(i) for i in 1:n])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:n])
end

"""
Attaches a second braid-like node with `m` arms to the leaves `k, k+1, k+2`. The planar embedding
is not guessed but searched for (all slot assignments, keeping only the
convention-clean ones) — there are always two, and both give the same
diagram.
"""
function _c16_with_gbraid(g::CircularGraph, k::Int, m::Int = 8)
    lets = letters(g.word); n = length(lets)
    ks = [mod1(k + i, n) for i in 0:2]
    c1, c2 = lets[ks[1]], lets[ks[2]]
    (c1 != c2 && lets[ks[3]] == c1) || return nothing
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

"Unfold both `:gbraid`s (variants `v1`, `v2`); first planar rotation."
function _c16_both_preimages(h::CircularGraph, v1::Int, v2::Int)
    for r2 in 0:7
        h2 = expand_gbraid(h, 2, v2; rot = r2)
        h2 === nothing && continue
        for r1 in 0:7
            h1 = expand_gbraid(h2, 1, v1; rot = r1)
            h1 === nothing && continue
            (is_wired(h1) && euler(h1) == 2 && isempty(check_wiring(h1))) || continue
            return h1
        end
    end
    return nothing
end

@testset "C16 gen12_on_gen12 (CircularGen12OnGen12.jl)" begin

    @testset "the equation, at all eight sites and in both colour pairs" begin
        hits = 0
        for (s, t) in ((1, 2), (2, 3)), k in 1:8
            g = _c16_with_gbraid(_c16_bare(s, t), k)
            g === nothing && continue
            res = _fr_gen12_on_gen12(g)
            res === nothing && continue
            ts = collect(pairs_of(res))
            length(ts) == 1 && ts[1][2] == 1 || continue
            h = ts[1][1]
            (count(nd -> nd.kind === :braid && arm_count(nd) == 8, h.nodes) == 1 &&
             count(nd -> arm_count(nd) == 3, h.nodes) == 2 &&
             length(h.nodes) == 3) || continue
            (join(letters(h.word)) == join(letters(g.word)) &&
             circular_degree(h) == circular_degree(g) &&
             euler(h) == 2 && isempty(check_wiring(h)) &&
             circular_weight(h) < circular_weight(g)) || continue
            hits += 1
        end
        @test hits == 16
    end

    @testset "the weights" begin
        g = _c16_with_gbraid(_c16_bare(1, 2), 1)
        h = first(pairs_of(_fr_gen12_on_gen12(g)))[1]
        @test circular_arm_weight(g) == (16, 2, 16)
        @test circular_arm_weight(h) == (8, 3, 14)
        @test circular_weight(h) < circular_weight(g)
        @test circular_degree(g) == -4
        @test circular_degree(h) == -4
    end

    @testset "the two trivalents are monochrome and carry the same colour" begin
        for k in (1, 2)
            g = _c16_with_gbraid(_c16_bare(1, 2), k)
            h = first(pairs_of(_fr_gen12_on_gen12(g)))[1]
            tri = [nd for nd in h.nodes if arm_count(nd) == 3]
            @test length(tri) == 2
            @test Set(arms(tri[1])) == Set(arms(tri[2]))
            @test length(Set(arms(tri[1]))) == 1
        end
        # namely the colour of the doubled boundary letter:
        # k = 1 -> boundary word 1212212122 (colour 2), k = 2 -> 1121211212 (colour 1)
        h1 = first(pairs_of(_fr_gen12_on_gen12(_c16_with_gbraid(_c16_bare(1, 2), 1))))[1]
        h2 = first(pairs_of(_fr_gen12_on_gen12(_c16_with_gbraid(_c16_bare(1, 2), 2))))[1]
        @test only(unique(arms(first(nd for nd in h1.nodes if arm_count(nd) == 3)))) == 2
        @test only(unique(arms(first(nd for nd in h2.nodes if arm_count(nd) == 3)))) == 1
    end

    # The next two blocks pin the 3-node form with EXACTLY ONE `:gbraid` as the
    # right-hand side. Under `CIRCULAR_TRIVALENT_MERGE = true` C19 pulls the
    # trivalents back in and the driver lands elsewhere, so the C16 surgery and
    # the 16/16 run are checked in the world `CIRCULAR_TRIVALENT_MERGE = false`:
    # C14 live, C19 dead, arm weight.
    _c16_direction = CIRCULAR_TRIVALENT_MERGE[]
    CIRCULAR_TRIVALENT_MERGE[] = false
    try

    @testset "GROUND TRUTH: the same as via the preimages without C16" begin
        without_c16 = [r for r in CIRCULAR_RULES if r.name !== :gen12_on_gen12]
        g0 = _c16_bare(1, 2)
        for k in (1, 2)
            g = _c16_with_gbraid(g0, k)
            h = first(pairs_of(_fr_gen12_on_gen12(g)))[1]
            targets = Any[]
            for v1 in 1:4, v2 in 1:4
                e = _c16_both_preimages(g, v1, v2)
                e === nothing && continue
                c, _ = reduce_circular(circular_decorated(e); rules = without_c16)
                for (tt, _) in pairs_of(c)
                    gg = tt.graph
                    (length(gg.nodes) == 3 &&
                     count(nd -> nd.kind === :braid && arm_count(nd) >= 8, gg.nodes) == 1) || continue
                    push!(targets, circular_canonical_key(gg))
                end
            end
            @test length(unique(targets)) == 1          # the preimages agree
            @test circular_canonical_key(h) == first(unique(targets))
        end
    end

    @testset "with C16 ALL 16 preimage combinations land on the right-hand side" begin
        # Without C16 (`without_c16`), 6 of 16 preimage combinations land on the
        # right-hand side (four classes); with C16 alone, 14 of 16. The two remaining
        # outliers are the case `:gbraid`–`:braid` sharing TWO edges instead of three;
        # with C18 (`gen12_two_edges`) also active they fall onto the same right-hand
        # side too — the run is 16/16, only ONE class. Removing C18 from the registry
        # reads 14 again.
        g = _c16_with_gbraid(_c16_bare(1, 2), 1)
        h = first(pairs_of(_fr_gen12_on_gen12(g)))[1]
        zk = circular_canonical_key(h)
        hits = 0
        for v1 in 1:4, v2 in 1:4
            e = _c16_both_preimages(g, v1, v2)
            e === nothing && continue
            c, _ = reduce_circular(circular_decorated(e))
            all(circular_canonical_key(tt.graph) == zk for (tt, _) in pairs_of(c)) && (hits += 1)
        end
        @test hits == 16

        without18 = [r for r in CIRCULAR_RULES if r.name !== :gen12_two_edges]
        t14 = 0
        for v1 in 1:4, v2 in 1:4
            e = _c16_both_preimages(g, v1, v2)
            e === nothing && continue
            c, _ = reduce_circular(circular_decorated(e); rules = without18)
            all(circular_canonical_key(tt.graph) == zk for (tt, _) in pairs_of(c)) && (t14 += 1)
        end
        @test t14 == 14
    end

    finally
        CIRCULAR_TRIVALENT_MERGE[] = _c16_direction
    end

    @testset "the driver also takes C16 by itself" begin
        g = _c16_with_gbraid(_c16_bare(1, 2), 1)
        c, hist = reduce_circular(circular_decorated(g))
        @test :gen12_on_gen12 in hist
        @test length(collect(pairs_of(c))) == 1
    end

    # ---- arm counts above 8 -------------------------------------------------
    #
    # Two gbraids of DIFFERENT size, `(10, 8)`: the rule is not restricted to
    # the equal 8-armed pair, and the surgery derives everything from
    # `p = nv + nb - 6`.
    @testset "C16 at a 10/8 gbraid pair" begin
        g = nothing
        for k in 1:10
            h = _c16_with_gbraid(_c16_bare(1, 2, 10), k, 8)
            h === nothing || (g = h; break)
        end
        @test g !== nothing
        @test sort([arm_count(nd) for nd in g.nodes]) == [8, 10]
        res = _fr_gen12_on_gen12(g)
        @test res !== nothing
        h = first(pairs_of(res))[1]
        @test sort([arm_count(nd) for nd in h.nodes]) == [3, 3, 10]
        @test euler(h) == 2 && isempty(check_wiring(h)) && is_wired(h)
        @test h.word == g.word
        @test circular_degree(h) == circular_degree(g)
        @test circular_arm_weight(h)[1] == circular_arm_weight(g)[1] - 8
        @test circular_weight(h) < circular_weight(g)
    end

    @testset "no false alarm" begin
        @test _fr_gen12_on_gen12(_c16_bare(1, 2)) === nothing      # only ONE :gbraid
        @test _fr_gen12_on_gen12(circular(braid(1, 2; m = 3))) === nothing
    end

    # ---- C6–C9 do NOT touch the 8-arm pattern --------------------------------
    #
    # Two 8-armed braid-like nodes with three shared edges — EXACTLY the situation in
    # which C7 (`braid_back`, two braids with three inner edges) and C8
    # (`braid_relation`) could formally match. C16 is responsible. A "gbraid" is
    # simply a `:braid` node with `arm_count >= 8`, so `kind` alone cannot separate
    # the two cases — the explicit `arm_count == 6` guard on C7/C8 does. This block
    # pins the behaviour down, without reference to the kind.
    @testset "C6–C9 do not fire at the 8-arm pair" begin
        for k in 1:8
            g = _c16_with_gbraid(_c16_bare(1, 2), k)
            g === nothing && continue
            @test arm_count(g.nodes[1]) == 8 && arm_count(g.nodes[2]) == 8
            @test DiagrammaticHecke._fr_braid_back(g) === nothing
            @test DiagrammaticHecke._fr_braid_back_at(g, 1, 2) === nothing
            @test DiagrammaticHecke._fr_braid_relation(g) === nothing
            @test DiagrammaticHecke._fr_dot_into_braid(g) === nothing
            # the other half of the check: C16 is responsible.
            @test _fr_gen12_on_gen12(g) !== nothing
        end
    end
end

