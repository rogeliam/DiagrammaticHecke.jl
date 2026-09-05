# test/circulargen12expand.jl — C13⁻¹ `expand_gbraid` (src/circular/rules/CircularGen12Expand.jl).
#
# The one test that matters is the ROUND TRIP: unfolding and merging again must give
# the same `:gbraid` — over all four variants and all eight rotations, in BOTH colour
# pairs (`{1,2}` and `{2,3}`). That also shows the frozen literals are genuine preimages
# and not merely planar fantasy constructs.
#
# The four literals fall into four classes of three, not twelve.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_gen12_merge, Edge, NodePort, Leaf

"""
The bare 8-armed `:gbraid` with colour pair `(s,t)`, all eight arms on the boundary.

⚠️ The leaves run **backwards** in slot order (slot 1 ↦ leaf 1, slot `i` ↦ leaf `10-i`)
— that is the house convention `check_wiring` tests as `:leaf_ring`. With the naive
assignment `slot i ↦ leaf i` even this base graph is non-planar (`euler = −4`, three
violations).
"""
function _ge_bare(s::Int, t::Int)
    a = [isodd(i) ? s : t for i in 1:8]
    lf(i) = i == 1 ? 1 : 10 - i
    word = [a[i] for i in 1:8][invperm([lf(i) for i in 1:8])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:8])
end

@testset "C13^-1 expand_gbraid (CircularGen12Expand.jl)" begin

    @testset "the base graph itself keeps the convention" begin
        for (s, t) in ((1, 2), (2, 3))
            g = _ge_bare(s, t)
            @test is_wired(g) && euler(g) == 2 && isempty(check_wiring(g))
            @test g.nodes[1].kind === :braid && arm_count(g.nodes[1]) == 8
            @test circular_degree(g) == -2
        end
    end

    @testset "four frozen preimages" begin
        @test length(GBRAID_PREIMAGES) == 4
        for (nds, inner, outer) in GBRAID_PREIMAGES
            @test length(nds) == 4                      # 2 braids + 2 trivalents
            @test count(a -> length(a) == 6, nds) == 2
            @test count(a -> length(a) == 3, nds) == 2
            @test length(inner) == 5                    # five inner edges
            @test length(outer) == 8                    # eight outer ports
        end
    end

    @testset "round trip {1,2}: unfold, merge, the same :gbraid again" begin
        g = _ge_bare(1, 2)
        k = circular_canonical_key(g)
        hits = 0
        for variant in 1:4, rot in 0:7
            h = expand_gbraid(g, 1, variant; rot = rot)
            @test h !== nothing
            r = _fr_gen12_merge(h)
            r === nothing && continue
            circular_canonical_key(first(pairs_of(r))[1]) == k && (hits += 1)
        end
        @test hits == 32                             # 4 variants x 8 rotations
    end

    @testset "round trip {2,3}: the same literals, a different colour pair" begin
        g = _ge_bare(2, 3)
        k = circular_canonical_key(g)
        hits = 0
        for variant in 1:4, rot in 0:7
            h = expand_gbraid(g, 1, variant; rot = rot)
            r = _fr_gen12_merge(h)
            r === nothing && continue
            circular_canonical_key(first(pairs_of(r))[1]) == k && (hits += 1)
        end
        @test hits == 32
    end

    @testset "the invariants of every unfolding" begin
        g = _ge_bare(1, 2)
        for variant in 1:4, rot in 0:7
            h = expand_gbraid(g, 1, variant; rot = rot)
            @test is_wired(h)
            @test euler(h) == 2
            @test isempty(check_wiring(h))
            @test join(letters(h.word)) == join(letters(g.word))   # boundary word
            @test circular_degree(h) == circular_degree(g)                   # degree
            @test length(h.nodes) == 4
            @test circular_weight(h) > circular_weight(g)                    # unfolds UP
        end
    end

    @testset "the four are pairwise distinct — and exhaustive" begin
        g = _ge_bare(1, 2)
        ks = [circular_canonical_key(expand_gbraid(g, 1, v; rot = 0)) for v in 1:4]
        @test length(unique(ks)) == 4

        # And there are no more than these four: all 4 x 8 unfoldings (every
        # variant at every rotation) together produce EXACTLY these same four
        # diagrams. On the bare, rotation-symmetric `:gbraid`, rotations only
        # permute the four variants — `rot` produces nothing new. (At
        # an embedded `:gbraid` this differs: there the surroundings fix
        # which rotation fits.)
        every = [circular_canonical_key(expand_gbraid(g, 1, v; rot = r)) for v in 1:4, r in 0:7]
        @test length(unique(every)) == 4
        @test Set(unique(every)) == Set(ks)
    end

    @testset "Randfaelle" begin
        # no :gbraid present
        g6 = circular(braid(1, 2; m = 3))
        @test expand_gbraid(g6, 1, 1) === nothing
        # variant outside 1:4
        g = _ge_bare(1, 2)
        @test_throws ArgumentError expand_gbraid(g, 1, 5)
        # the short form finds the :gbraid itself
        @test expand_gbraid(g, 2) !== nothing
    end
end

# ---- the TOWER preimage for an arbitrary arm count ---------------------------
#
# From ten arms on, nothing is searched for any more but constructed: the position of
# the fusion cluster is forced, so the preimage can be written down. The
# test is again the ROUND TRIP — unfold, merge with C13, and get back
# the same node.

"The bare `2m`-armed `:gbraid` with colour pair `(s,t)` (leaves backwards)."
function _ge_bare_n(n::Int, s::Int = 1, t::Int = 2)
    a = [isodd(i) ? s : t for i in 1:n]
    lf(i) = mod1(2 - i, n)
    word = a[invperm([lf(i) for i in 1:n])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:n])
end

@testset "C13^-1 tower preimages for 2k (CircularGen12Expand.jl)" begin

    @testset "gbraid_tower_decompositions" begin
        @test gbraid_tower_decompositions(4) == [(3, 3)]
        @test gbraid_tower_decompositions(5) == [(3, 4), (4, 3)]
        @test gbraid_tower_decompositions(6) == [(3, 5), (4, 4), (5, 3)]
        @test all(k + l - 2 == 7 for (k, l) in gbraid_tower_decompositions(7))
    end

    @testset "the cluster is a diagram — and it contracts correctly" begin
        for (k, l) in [(3, 3), (4, 3), (3, 4), (4, 4), (5, 3), (3, 5)]
            m = k + l - 2
            g = gbraid_tower_cluster(k, l)
            @test is_wired(g)
            @test euler(g) == 2
            @test isempty(check_wiring(g))
            @test length(g.nodes) == 4
            @test sort([arm_count(nd) for nd in g.nodes]) == sort([2k, 2l, 3, 3])
            # both trivalents of the same colour, that of the outer slots
            triv = [nd for nd in g.nodes if arm_count(nd) == 3]
            @test length(unique(arms(nd)[1] for nd in triv)) == 1
            r = contract_circular_cluster(g, [1, 2, 3, 4])
            @test r !== nothing
            (h, v) = r
            @test h.nodes[v].kind === :braid && arm_count(h.nodes[v]) >= 8
            @test arm_count(h.nodes[v]) == 2m
            @test circular_canonical_key(h) == circular_canonical_key(_ge_bare_n(2m))
            @test circular_degree(h) == circular_degree(g)          # degree preserved (Q6!)
        end
    end

    @testset "round trip: unfold, merge, the same node again" begin
        for (s, t) in ((1, 2), (2, 3))
            for m in 4:6
                g = _ge_bare_n(2m, s, t)
                key = circular_canonical_key(g)
                for (k, l) in gbraid_tower_decompositions(m), r in 0:(2m - 1)
                    h = expand_gbraid_tower(g, 1, k, l; rot = r)
                    @test h !== nothing
                    @test is_wired(h)
                    @test euler(h) == 2
                    @test isempty(check_wiring(h))
                    @test join(letters(h.word)) == join(letters(g.word))
                    @test circular_degree(h) == circular_degree(g)
                    @test circular_weight(h) > circular_weight(g)          # unfolds UP
                    res = _fr_gen12_merge(h)
                    @test res !== nothing
                    @test circular_canonical_key(first(pairs_of(res))[1]) == key
                end
            end
        end
    end

    @testset "expand_gbraid_variants — all unfoldings, deduplicated" begin
        # 8 arms: the four frozen preimages AND the tower cluster (3,3). The latter is
        # one of the four (trivalents of equal colour), so 4 remain.
        g8 = _ge_bare_n(8)
        vs8 = expand_gbraid_variants(g8, 1)
        @test length(vs8) == 4
        @test circular_canonical_key(gbraid_tower_cluster(3, 3)) in
              Set(circular_canonical_key(h) for h in vs8)

        # 10 arms: (3,4) and (4,3) — at the ROTATION-SYMMETRIC bare node these collapse
        # into ONE diagram (the two differ only by a reflection/rotation of the
        # boundary).
        vs10 = expand_gbraid_variants(_ge_bare_n(10), 1)
        @test !isempty(vs10)
        for h in vs10
            @test euler(h) == 2 && isempty(check_wiring(h))
            @test circular_degree(h) == -4
        end

        vs12 = expand_gbraid_variants(_ge_bare_n(12), 1)
        @test !isempty(vs12)
        for h in vs12
            @test circular_degree(h) == -6
        end

        # no :gbraid ⇒ empty
        @test isempty(expand_gbraid_variants(circular(braid(1, 2; m = 3)), 1))
    end

    @testset "Randfaelle" begin
        @test_throws ArgumentError gbraid_tower_literal(2, 3)
        # wrong arm count ⇒ nothing
        @test expand_gbraid_tower(_ge_bare_n(8), 1, 4, 3) === nothing
    end
end

