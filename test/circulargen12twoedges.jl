# test/circulargen12twoedges.jl — C18 `gen12_two_edges` (src/circular/rules/CircularGen12TwoEdges.jl).
#
# The rule: an 8-armed `:gbraid` and a `:braid` with EXACTLY TWO shared edges. It has
# no surgery of its own — the `:gbraid` is unfolded to a preimage via `expand_gbraid`
# (C13⁻¹), the existing rules compute there, and only what strictly lowers `circular_weight`
# is taken.
#
# The fixture is the outlier of the C16 run: the pattern "two `:gbraid`,
# three edges" unfolded over both preimages and reduced without C18 — two of the 16
# combinations remain at `:braid` + `:gbraid` + two trivalents. Those are exactly what
# C18 catches.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_gen12_two_edges, _circular_gen12_braid_two_edges,
                     _circular_gen12_partner_pairs, _circular_gen12_channel_pairs,
                     _circular_gen12_expansions, _circular_braid_channels, _circular_braid_pairs,
                     _fr_braid_relation, _fr_braid_relation_gen,
                     _fr_gen12_on_gen12, Edge, NodePort, Leaf, CIRCULAR_RULES

"The bare 8-armed `:gbraid` (as in test/circulargen12ongen12.jl)."
function _c18_bare(s::Int, t::Int)
    a = [isodd(i) ? s : t for i in 1:8]; lf(i) = i == 1 ? 1 : 10 - i
    word = [a[i] for i in 1:8][invperm([lf(i) for i in 1:8])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:8])
end

"A second 8-armed `:gbraid` at the leaves `k, k+1, k+2` (searched for planarity)."
function _c18_with_gbraid(g::CircularGraph, k::Int)
    lets = letters(g.word); n = length(lets)
    ks = [mod1(k + i, n) for i in 0:2]
    c1, c2 = lets[ks[1]], lets[ks[2]]
    (c1 != c2 && lets[ks[3]] == c1) || return nothing
    for start in 1:8, dir in (1, -1), newdir in (1, -1)
        barms = [isodd(i) ? c1 : c2 for i in 1:8]
        sl = [mod1(start + dir * (i - 1), 8) for i in 1:3]
        all(barms[sl[i]] == lets[ks[i]] for i in 1:3) || continue
        rest = [s for s in 1:8 if !(s in sl)]
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
            push!(edges, Edge(barms[s], Leaf(length(neu) - 5 + j), NodePort(bi, s)))
        end
        h = CircularGraph(CircularWord(neu), nodes, edges)
        (is_wired(h) && euler(h) == 2 && isempty(check_wiring(h))) && return h
    end
    return nothing
end

"Unfold both `:gbraid`s (variants `v1`, `v2`); first planar rotation."
function _c18_both_preimages(h::CircularGraph, v1::Int, v2::Int)
    for r2 in 0:7
        h2 = expand_gbraid(h, 2, v2; rot = r2); h2 === nothing && continue
        for r1 in 0:7
            h1 = expand_gbraid(h2, 1, v1; rot = r1); h1 === nothing && continue
            (is_wired(h1) && euler(h1) == 2 && isempty(check_wiring(h1))) || continue
            return h1
        end
    end
    return nothing
end

"""
Attach a partner node with arms `barms` to the `j` consecutive leaves `k, …, k+j-1`
(a generalisation of `_c18_with_gbraid`; the fixtures of the widened matcher).
The first planar wiring, otherwise `nothing`.
"""
function _c18_with_partner(g::CircularGraph, k::Int, j::Int, barms::Vector{Int})
    lets = letters(g.word); n = length(lets); m = length(barms)
    ks = [mod1(k + i - 1, n) for i in 1:j]
    for start in 1:m, dir in (1, -1), newdir in (1, -1)
        sl = [mod1(start + dir * (i - 1), m) for i in 1:j]
        all(barms[sl[i]] == lets[ks[i]] for i in 1:j) || continue
        rest = [s for s in 1:m if !(s in sl)]
        ns = newdir == 1 ? rest : reverse(rest)
        rw = [lets[mod1(ks[end] + i, n)] for i in 1:(n - j)]
        neu = vcat(rw, [barms[s] for s in ns])
        leafmap = Dict(mod1(ks[end] + i, n) => i for i in 1:(n - j))
        nodes = vcat(copy(g.nodes), CircularNode[circular_node(barms)]); bi = length(nodes)
        slot_of_old = Dict(ks[i] => sl[i] for i in 1:j)
        rp(p) = p isa Leaf ? (haskey(slot_of_old, p.k) ? NodePort(bi, slot_of_old[p.k]) :
                                                         Leaf(leafmap[p.k])) : p
        edges = Edge[Edge(e.colour, rp(e.a), rp(e.b)) for e in g.edges]
        for (jj, s) in enumerate(ns)
            push!(edges, Edge(barms[s], Leaf(n - j + jj), NodePort(bi, s)))
        end
        h = CircularGraph(CircularWord(neu), nodes, edges)
        (is_wired(h) && euler(h) == 2 && isempty(check_wiring(h))) && return h
    end
    return nothing
end

"""
The outlier for the attachment site `k`: the C16 pattern unfolded over the preimages
`(v1,v2)`, reduced without C18, and of that the 4-node term.
"""
function _c18_outlier(k::Int, v1::Int, v2::Int)
    without_c18 = [r for r in CIRCULAR_RULES if r.name !== :gen12_two_edges]
    g = _c18_with_gbraid(_c18_bare(1, 2), k)
    e = _c18_both_preimages(g, v1, v2)
    e === nothing && return nothing
    c, _ = reduce_circular(circular_decorated(e); rules = without_c18)
    for (t, _) in pairs_of(c)
        length(t.graph.nodes) == 4 && return t.graph
    end
    return nothing
end

# The fixture `A` is built through the DRIVER (`_c18_outlier`), which contains
# C19 `trivalent_into_gen12`. Under `CIRCULAR_TRIVALENT_MERGE = true` C19 pulls
# the two trivalents into the 8-armed node while building, giving a 12-armed node
# with two dots, and the two-edge site C18 looks for no longer exists. The outlier
# is by construction an object of the world `CIRCULAR_TRIVALENT_MERGE = false`:
# C14 live, C19 dead, arm weight — so it is built and checked there. C18 itself is
# unaffected and keeps firing under `true`.
_c18_direction = CIRCULAR_TRIVALENT_MERGE[]
CIRCULAR_TRIVALENT_MERGE[] = false
try

@testset "C18 gen12_two_edges (CircularGen12TwoEdges.jl)" begin

    A = _c18_outlier(1, 1, 1)          # k = 1, preimages (1,1)

    @testset "the fixture really is the two-edge case" begin
        @test A !== nothing
        @test sort([(nd.kind, arm_count(nd)) for nd in A.nodes], by = string) ==
              [(:braid, 6), (:braid, 8), (:mixed, 3), (:mono, 3)]
        @test circular_weight(A) == (14, 4, 20)
        @test length(_circular_gen12_braid_two_edges(A)) == 1        # exactly ONE pair
    end

    @testset "the rule fires and lowers the weight" begin
        res = _fr_gen12_two_edges(A)
        @test res !== nothing
        ts = collect(pairs_of(res))
        @test length(ts) == 1
        @test ts[1][2] == 1
        h = ts[1][1]
        @test sort([(nd.kind, arm_count(nd)) for nd in h.nodes], by = string) ==
              [(:braid, 8), (:mono, 3), (:mono, 3)]
        @test circular_weight(h) == (8, 3, 14)
        @test circular_weight(h) < circular_weight(A)
        @test join(letters(h.word)) == join(letters(A.word))     # boundary unchanged
        @test circular_degree(h) == circular_degree(A)
        @test euler(h) == 2 && isempty(check_wiring(h))
    end

    @testset "GROUND TRUTH: this is C16's right-hand side" begin
        g = _c18_with_gbraid(_c18_bare(1, 2), 1)
        zk = circular_canonical_key(first(pairs_of(_fr_gen12_on_gen12(g)))[1])
        @test circular_canonical_key(first(pairs_of(_fr_gen12_two_edges(A)))[1]) == zk
    end

    @testset "WELL DEFINED: all weight-lowering preimages agree" begin
        # 4 variants × 8 rotations; 32 of them planar, 16 weight-lowering, and those
        # give ONE key.
        without_c18 = [r for r in CIRCULAR_RULES if r.name !== :gen12_two_edges]
        w0 = circular_weight(A)
        keys_ = Any[]; planar = 0
        v = only(x[1] for x in _circular_gen12_braid_two_edges(A))
        for variant in 1:4, rot in 0:7
            e = expand_gbraid(A, v, variant; rot = rot)
            e === nothing && continue
            (is_wired(e) && euler(e) == 2 && isempty(check_wiring(e))) || continue
            planar += 1
            c, _ = reduce_circular(circular_decorated(e); rules = without_c18)
            all(circular_weight(t.graph) < w0 for (t, _) in pairs_of(c)) || continue
            push!(keys_, sort([circular_canonical_key(t.graph) for (t, _) in pairs_of(c)], by = string))
        end
        @test planar == 32
        @test length(keys_) == 16
        @test length(unique(keys_)) == 1
    end

    @testset "the driver takes C18 by itself" begin
        c, hist = reduce_circular(circular_decorated(A))
        @test :gen12_two_edges in hist
        @test length(collect(pairs_of(c))) == 1
    end

    @testset "widening the matcher: 1/4 edges and :gbraid partners" begin
        # The synthetic minimal patterns are CIRCULAR_RULES fixed points, and NONE
        # of the 32 planar unfoldings lowers the weight. The widened matcher SEES
        # them, the weight guard refuses — the rule is a no-op there.
        #
        # `h4` is the exception, and not because of C18: four shared edges are the
        # site of C22 `gen12_on_gen12_null`, which stands before C18 and makes the
        # term 0. That is a hit of C22's FIRST branch, so it is a proof by
        # construction — the C15/C16 surgery on either three-edge window really is
        # carried out and leaves a wired self-loop (checked below).
        g0 = _c18_bare(1, 2)
        h1 = _c18_with_partner(g0, 1, 1, [isodd(i) ? 1 : 2 for i in 1:6])  # one edge
        h4 = _c18_with_partner(g0, 1, 4, [isodd(i) ? 1 : 2 for i in 1:6])  # four edges
        hg = _c18_with_partner(g0, 1, 2, [isodd(i) ? 1 : 2 for i in 1:8])  # gbraid–gbraid, two
        for (h, np, narms) in ((h1, 1, 6), (h4, 4, 6), (hg, 2, 8))
            @test h !== nothing
            trip = _circular_gen12_partner_pairs(h)
            @test (1, 2, np) in trip                       # the matcher sees the pair
            narms == 8 && @test (2, 1, np) in trip         # … in both orderings
            @test _fr_gen12_two_edges(h) === nothing       # guard: nothing lighter
            c, hist = reduce_circular(circular_decorated(h))
            if np == 4
                @test hist == [:gen12_on_gen12_null]       # C22, not C18
                @test isempty(collect(pairs_of(c)))        # the term is 0
            else
                @test isempty(hist)                        # driver fixed point
                @test length(collect(pairs_of(c))) == 1
            end
        end

        # C22 on `h4` in detail: 8 + 6 arms, four shared edges, two three-edge
        # windows per node ordering. Every window's surgery is planar, wired and
        # carries the self-loop that C1 turns into 0.
        for (v, w) in ((1, 2), (2, 1))
            @test length(DiagrammaticHecke._circular_shared_slot_pairs(h4, v, w)) == 4
            wins = DiagrammaticHecke._circular_glue3_windows(h4, v, w)
            @test length(wins) == 2
            for win in wins
                hh = DiagrammaticHecke._circular_glue3_rewire_at(h4, v, w, win)
                @test hh !== nothing
                @test euler(hh) == 2 && isempty(check_wiring(hh))
                @test DiagrammaticHecke._circular_has_self_loop(hh)
            end
            @test DiagrammaticHecke._circular_has_glue4_window(h4, v, w)
        end
        # the two-edge matcher is a filter of the wider triples
        @test _circular_gen12_braid_two_edges(A) == [(v, b) for (v, b, n) in
              _circular_gen12_partner_pairs(A) if n == 2]
    end

    @testset "channel sites + targeted preimage choice" begin
        # At the fixture there sits an F2-3 channel (`nc = 3`)
        # between `:gbraid` and `:braid`, which
        # `_circular_braid_pairs` never looks at, because it requires `:braid` at
        # both ends. The widened site search sees it.
        v = only(x[1] for x in _circular_gen12_braid_two_edges(A))
        b = only(x[2] for x in _circular_gen12_braid_two_edges(A))
        chs = _circular_braid_channels(A, v, b)
        @test !isempty(chs)                       # channel detection is generic
        @test any(ch -> ch.nc == 3, chs)          # t = nc - 2 = 1  ⇒  F2-3
        @test isempty(_circular_braid_pairs(A))        # … that filter sees nothing
        @test (v, b) in _circular_gen12_channel_pairs(A)

        # The choice of preimage really is a choice: of the 32 planar
        # unfoldings, exactly 8 allow F2-3 DIRECTLY — and those come first.
        exps = _circular_gen12_expansions(A, v)
        @test length(exps) == 32
        f23 = [_fr_braid_relation(e) !== nothing || _fr_braid_relation_gen(e) !== nothing
               for e in exps]
        @test count(f23) == 8
        @test all(f23[1:8]) && !any(f23[9:end])   # the F2-3-capable ones first

        # … and all eight compute to the SAME right-hand side (the general rule).
        without_c18 = [r for r in CIRCULAR_RULES if r.name !== :gen12_two_edges]
        ks = [sort([circular_canonical_key(t.graph)
                    for (t, _) in pairs_of(reduce_circular(circular_decorated(e); rules = without_c18)[1])],
                   by = string) for e in exps[1:8]]
        @test length(unique(ks)) == 1
    end

    @testset "no false alarm" begin
        # only ONE :gbraid, no :braid beside it
        @test _fr_gen12_two_edges(_c18_bare(1, 2)) === nothing
        # the C16 pattern: two :gbraid, THREE edges — not this pattern
        @test _fr_gen12_two_edges(_c18_with_gbraid(_c18_bare(1, 2), 1)) === nothing
        # an ordinary braid on its own
        @test _fr_gen12_two_edges(circular(braid(1, 2; m = 3))) === nothing
    end
end

finally
    CIRCULAR_TRIVALENT_MERGE[] = _c18_direction
end
