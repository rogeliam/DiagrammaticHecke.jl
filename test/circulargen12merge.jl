# test/circulargen12merge.jl — C13 `gen12_merge` (src/circular/rules/CircularGen12Merge.jl).
#
# The fixtures are built with the planar constructor and frozen here as literals —
# each is wiring-clean with `euler == 2`.
#
# They are checked against the two invariants EVERY rule must preserve (the boundary
# word and `circular_degree`), plus the termination evidence (`circular_weight` falls) — without
# which the rule could not stand in `CIRCULAR_RULES`.
#
# The NEGATIVE test matters as much here as the positive ones: `nc ≥ 4` in ONE
# channel must NOT become a `:gbraid`, but must land at C8.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_gen12_merge, _circular_gen12_clusters, _circular_first_rule_match,
      _fr_braid_relation_at, _circular_is_pure_trivalent, NodePort, Edge, Leaf, CircularNode, circular_node,
      CircularWord

# Classic case C: both outer channels a PURE trivalent [3,3,3].
_g12_caseC_classic() = CircularGraph(CircularWord([2, 3, 2, 3, 2, 3, 2, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]),
            circular_node([3, 3, 3]), circular_node([3, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 3), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(4, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(4, 2), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(4, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(1, 2)),
    ])

# Both channels a general-13 node, side-by-side variant [1,1,3,3,3].
_g12_both_beside() = CircularGraph(CircularWord([1, 1, 3, 2, 3, 2, 1, 1, 3, 2, 3, 2]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]),
            circular_node([1, 1, 3, 3, 3]), circular_node([1, 1, 3, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(3, 4)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 3), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(4, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(4, 4), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(4, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(4, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(4, 5)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(3, 5)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(10), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(11), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(12), DiagrammaticHecke.NodePort(1, 1)),
    ])

# Both channels a general-13 node, CROSSING variant [1,3,1,3,3].
_g12_both_crossing() = CircularGraph(CircularWord([1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 1, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]),
            circular_node([1, 3, 1, 3, 3]), circular_node([1, 3, 1, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(3, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(4, 4)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(4, 5), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(4, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(3, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(10), DiagrammaticHecke.NodePort(1, 1)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(11), DiagrammaticHecke.NodePort(4, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(12), DiagrammaticHecke.NodePort(4, 2)),
    ])

# nc = 4 in ONE channel — the case that must NOT become a :gbraid.
_g12_nc4_one_channel() = CircularGraph(CircularWord([1, 1, 3, 3, 2, 3, 2, 2, 3, 2]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]),
            circular_node([1, 1, 3, 3, 3, 3])],
    DiagrammaticHecke.Edge[
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 4), DiagrammaticHecke.NodePort(2, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.NodePort(1, 5), DiagrammaticHecke.NodePort(2, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(1, 6), DiagrammaticHecke.NodePort(3, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.NodePort(3, 4), DiagrammaticHecke.NodePort(2, 6)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(1), DiagrammaticHecke.NodePort(3, 2)),
        DiagrammaticHecke.Edge(1, DiagrammaticHecke.Leaf(2), DiagrammaticHecke.NodePort(3, 1)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(3), DiagrammaticHecke.NodePort(3, 6)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(4), DiagrammaticHecke.NodePort(3, 5)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(5), DiagrammaticHecke.NodePort(2, 5)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(6), DiagrammaticHecke.NodePort(2, 4)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(7), DiagrammaticHecke.NodePort(2, 3)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(8), DiagrammaticHecke.NodePort(1, 3)),
        DiagrammaticHecke.Edge(3, DiagrammaticHecke.Leaf(9), DiagrammaticHecke.NodePort(1, 2)),
        DiagrammaticHecke.Edge(2, DiagrammaticHecke.Leaf(10), DiagrammaticHecke.NodePort(1, 1)),
    ])

"The single term of a rule application, checked against the invariants."
function _g12_einziger_term(g::CircularGraph, res)
    @test res !== nothing
    ts = collect(pairs_of(res))
    @test length(ts) == 1
    h = ts[1][1]
    @test ts[1][2] == 1
    @test join(letters(h.word)) == join(letters(g.word))     # boundary word preserved
    @test circular_degree(h) == circular_degree(g)                     # degree preserved
    @test euler(h) == 2
    @test isempty(check_wiring(h))
    @test circular_weight(h) < circular_weight(g)                      # termination
    return h
end

"The :gbraid in the result."
function _g12_gbraid(h::CircularGraph)
    gb = [nd for nd in h.nodes if nd.kind === :braid && arm_count(nd) >= 8]
    @test length(gb) == 1
    @test arm_count(gb[1]) == 8
    @test circular_degree(gb[1]) == -2
    # alternating 2,3,2,3,…  — the defining property
    a = collect(arms(gb[1]))
    @test all(a[i] != a[mod1(i + 1, 8)] for i in eachindex(a))
    @test Set(a) == Set([2, 3])
    return gb[1]
end

@testset "C13 gen12_merge (CircularGen12Merge.jl)" begin

    @testset "the fixtures are wiring-clean" begin
        for f in (_g12_caseC_classic, _g12_both_beside,
                  _g12_both_crossing, _g12_nc4_one_channel)
            g = f()
            @test is_wired(g)
            @test euler(g) == 2
            @test isempty(check_wiring(g))
            @test circular_degree(g) == -2          # t = 2 in each of the four
        end
    end

    @testset "case C classic: 2 braids + 2 pure trivalents -> :gbraid" begin
        g = _g12_caseC_classic()
        # The coarse filter finds exactly ONE cluster.
        @test length(_circular_gen12_clusters(g)) == 1
        h = _g12_einziger_term(g, _fr_gen12_merge(g))
        @test length(h.nodes) == 1             # only the :gbraid left
        _g12_gbraid(h)
        @test join(letters(h.word)) == "23232323"
        @test circular_arm_weight(g) == (12, 4, 18)
        @test circular_arm_weight(h) == (8, 1, 8)
        @test circular_weight(h) < circular_weight(g)
    end

    @testset "both channels a general-13 node (channel normalisation)" begin
        # Without normalisation there is NO cluster — the coarse filter sees no pure
        # trivalents. Only splitting off in BOTH channels makes one.
        for f in (_g12_both_beside, _g12_both_crossing)
            g = f()
            @test isempty(_circular_gen12_clusters(g))
            h = _g12_einziger_term(g, _fr_gen12_merge(g))
            _g12_gbraid(h)
            # the :gbraid plus the two trimmed 13-nodes
            @test length(h.nodes) == 3
            @test count(nd -> arm_count(nd) == 4 && Set(arms(nd)) == Set([1, 3]),
                        h.nodes) == 2
            @test circular_arm_weight(g) == (12, 4, 22)
            @test circular_arm_weight(h) == (8, 3, 16)
            @test circular_weight(h) < circular_weight(g)
        end
    end

    @testset "these cases were FIXED POINTS and now fire :gen12_merge" begin
        for f in (_g12_caseC_classic, _g12_both_beside, _g12_both_crossing)
            g = f()
            hit = _circular_first_rule_match(circular_decorated(g), CIRCULAR_RULES)
            @test hit !== nothing
            @test hit[1] === :gen12_merge
        end
    end

    @testset "NEGATIVE: nc >= 4 in ONE channel does not become a :gbraid" begin
        g = _g12_nc4_one_channel()
        # the contraction runs but gives eight NON-alternating arms — C13
        # must not match here.
        @test _fr_gen12_merge(g) === nothing
        # C8 remains responsible, and it fires (CircularBraidChannels.jl).
        hit = _circular_first_rule_match(circular_decorated(g), CIRCULAR_RULES)
        @test hit !== nothing
        @test hit[1] === :braid_relation
        @test !any(nd.kind === :braid && arm_count(nd) >= 8
                   for (t, _) in pairs_of(hit[2]) for nd in t.graph.nodes)
    end

    @testset "contract_circular_cluster: building blocks" begin
        g = _g12_caseC_classic()
        r = contract_circular_cluster(g, [1, 2, 3, 4])
        @test r !== nothing
        (h, v) = r
        @test h.nodes[v].kind === :braid && arm_count(h.nodes[v]) >= 8
        @test arm_count(h.nodes[v]) == 8
        # A cluster that is not a disc gives `nothing`.
        @test contract_circular_cluster(g, [1, 3]) === nothing
    end

    @testset "the :gbraid is a fixed point afterwards (nothing eats it further)" begin
        g = _g12_caseC_classic()
        h = collect(pairs_of(_fr_gen12_merge(g)))[1][1]
        @test _circular_first_rule_match(circular_decorated(h), CIRCULAR_RULES) === nothing
        # and reduce_circular_full runs through instead of throwing in merge_at_edge
        c = reduce_circular_full(_g12_both_crossing())
        @test length(collect(pairs_of(c))) == 1
        @test any(nd.kind === :braid && arm_count(nd) >= 8
                  for (t, _) in pairs_of(c) for nd in t.graph.nodes)
    end

    # ---- C13 for arbitrary arm counts --------------------------------------
    #
    # `2k` and `2l`, triply connected with two trivalents outside, give `2(k+l-2)`.
    # The fixture is the tower cluster itself (`gbraid_tower_cluster`), the only
    # planar situation of that shape.
    @testset "C13 generalised: 2k + 2l -> 2(k+l-2)" begin
        for (k, l) in [(3, 3), (4, 3), (3, 4), (4, 4), (5, 3), (3, 5)]
            m = k + l - 2
            g = gbraid_tower_cluster(k, l)
            @test is_wired(g)
            @test euler(g) == 2
            @test isempty(check_wiring(g))
            @test circular_degree(g) == 6 - 2m               # (6-2k)+(6-2l)-1-1
            @test !isempty(_circular_gen12_clusters(g))

            h = _g12_einziger_term(g, _fr_gen12_merge(g))
            @test length(h.nodes) == 1
            @test h.nodes[1].kind === :braid && arm_count(h.nodes[1]) >= 8
            @test arm_count(h.nodes[1]) == 2m
            @test circular_degree(h.nodes[1]) == 6 - 2m
            # termination: component 1 of the ARM reading is the arm sum over the
            # braid-like nodes, and it ALWAYS falls by exactly 4. The active
            # reading falls too, which is what the driver checks.
            @test circular_arm_weight(g)[1] == 2k + 2l
            @test circular_arm_weight(h)[1] == 2m
            @test circular_weight(h) < circular_weight(g)
        end
    end

    # The 8-armed case is the base case: the tower cluster (3,3) in the colour pair
    # {2,3} is exactly the classical fixture up to `circular_canonical_key`, and both
    # yield the same `:gbraid`.
    @testset "(3,3) matches the classical fixture" begin
        gt = gbraid_tower_cluster(3, 3; s = 2, t = 3)
        @test join(letters(gt.word)) == join(letters(_g12_caseC_classic().word))
        @test sort([arm_count(nd) for nd in gt.nodes]) ==
              sort([arm_count(nd) for nd in _g12_caseC_classic().nodes])
        h  = collect(pairs_of(_fr_gen12_merge(gt)))[1][1]
        h0 = collect(pairs_of(_fr_gen12_merge(_g12_caseC_classic())))[1][1]
        @test circular_canonical_key(h) == circular_canonical_key(h0)
        @test arm_count(h.nodes[1]) == 8
        @test circular_degree(h.nodes[1]) == -2
    end
end

# ---- the double-tower induction ----------------------------------------------
#
# With the generalised C13, ALL FIVE rotations of the double tower `T` reduce via
# `reduce_to_circular_leave` to ONE node (the bare 10-armed `:gbraid`), coefficient
# `+1` throughout: `ρ²(π₅) = π₅` exactly and `ρσ(π₅) = +π₅`, for 10 as well as for 8
# arms.
@testset "double tower: all five rotations collapse onto ONE node" begin
    function _tower_merge_at(v, k)
        ps = MorphismGraph[identity_strand(v[j]) for j in 1:k-1]
        push!(ps, merge_morphism(v[k]))
        append!(ps, [identity_strand(v[j]) for j in k+2:length(v)])
        foldl(tensor, ps)
    end
    _tower_step((v, art, i, j)) = art === :braid ? braid_move_morphism(v, i, j) : _tower_merge_at(v, i)
    TOWER = [([1,2,1,2,1,2,1], :braid, 3, 5), ([1,2,2,1,2,2,1], :merge, 2, 0),
            ([1,2,1,2,2,1],   :merge, 4, 0), ([1,2,1,2,1],     :braid, 2, 4),
            ([1,1,2,1,1],     :merge, 1, 0), ([1,2,1,1],       :merge, 3, 0),
            ([1,2,1],         :braid, 1, 3)]
    TM = circular_morphism(foldl(compose, map(_tower_step, TOWER))); T = TM.graph
    @test length(T.nodes) == 7
    @test circular_degree(T) == -4

    sw(c) = c == 1 ? 2 : c == 2 ? 1 : c
    function _tower_rotate(g::CircularGraph, r::Int; swap::Bool = false)
        n = length(letters(g.word)); cols = collect(letters(g.word))
        neu = [(swap ? sw : identity)(cols[mod1(k + r, n)]) for k in 1:n]
        neu == cols || error("the boundary word changes")
        rp(p) = p isa Leaf ? Leaf(mod1(p.k - r, n)) : p
        nodes = CircularNode[swap ? circular_node(sw.(nd.arms)) : nd for nd in g.nodes]
        CircularGraph(g.word, nodes, Edge[Edge(swap ? sw(e.colour) : e.colour, rp(e.a), rp(e.b)) for e in g.edges])
    end

    results = Dict{String,Tuple{Any,Any}}()
    for (nm, g) in [("T", T), ("rho2T", _tower_rotate(T, 2)), ("rho4T", _tower_rotate(T, 4)),
                    ("rhosT", _tower_rotate(T, 1; swap = true)), ("rho3sT", _tower_rotate(T, 3; swap = true))]
        @test euler(g) == 2
        @test isempty(check_wiring(g))
        fdm = CircularDecoratedMorphism(CircularMorphismGraph(g, TM.cut1, TM.cut2),
                                   fill(one(SoergelPoly), region_count(g)))
        c = reduce_to_circular_leave(fdm)
        terms = collect(pairs_of(c))
        @test length(terms) == 1
        (d, coef) = terms[1]
        h = hasproperty(d, :graph) ? d.graph : d
        @test h.nodes[1].kind === :braid && arm_count(h.nodes[1]) >= 8
        @test arm_count(h.nodes[1]) == 10
        @test circular_degree(h) == -4
        results[nm] = (coef, circular_canonical_key(h))
    end
    @test length(Set(v[2] for v in values(results))) == 1     # one node
    @test all(v[1] == 1 for v in values(results))              # coefficient +1 everywhere
end

# ---- both collapse routes of the double tower are confluent ------------------
#
# `T` (7 nodes) has EXACTLY two four-node clusters contracting to an 8-armed
# `:gbraid`: {1,4,2,3} ("lower") and {4,7,5,6} ("upper"). Fully contracting either
# remainder must lead to the same bare 10-armed `:gbraid` (associativity of the
# fusion), and every genuine 8-armed preimage
# (`GBRAID_PREIMAGES`, `expand_gbraid`) in place of the collapsed intermediate node must
# lead back to it as well.
@testset "double tower: lower-first and upper-first collapse confluently" begin
    function _tower_merge_at2(v, k)
        ps = MorphismGraph[identity_strand(v[j]) for j in 1:k-1]
        push!(ps, merge_morphism(v[k]))
        append!(ps, [identity_strand(v[j]) for j in k+2:length(v)])
        foldl(tensor, ps)
    end
    _tower_step2((v, art, i, j)) = art === :braid ? braid_move_morphism(v, i, j) : _tower_merge_at2(v, i)
    TOWER = [([1,2,1,2,1,2,1], :braid, 3, 5), ([1,2,2,1,2,2,1], :merge, 2, 0),
            ([1,2,1,2,2,1],   :merge, 4, 0), ([1,2,1,2,1],     :braid, 2, 4),
            ([1,1,2,1,1],     :merge, 1, 0), ([1,2,1,1],       :merge, 3, 0),
            ([1,2,1],         :braid, 1, 3)]
    T = circular_morphism(foldl(compose, map(_tower_step2, TOWER))).graph

    (gA, vA) = contract_circular_cluster(T, [1, 4, 2, 3])
    @test gA.nodes[vA].kind === :braid && arm_count(gA.nodes[vA]) == 8
    (gB, vB) = contract_circular_cluster(T, [4, 7, 5, 6])
    @test gB.nodes[vB].kind === :braid && arm_count(gB.nodes[vB]) == 8

    (g1A, v1A) = contract_circular_cluster(gA, collect(1:length(gA.nodes)))
    (g1B, v1B) = contract_circular_cluster(gB, collect(1:length(gB.nodes)))
    @test g1A.nodes[v1A].kind === :braid && arm_count(g1A.nodes[v1A]) == 10
    @test g1B.nodes[v1B].kind === :braid && arm_count(g1B.nodes[v1B]) == 10
    @test circular_canonical_key(g1A) == circular_canonical_key(g1B)

    # and reinserting EVERY genuine 8-armed preimage in place of vA leads
    # back to the same 10-armed node (32 = 4 variants x 8 rotations).
    hits = 0
    for variant in 1:4, rot in 0:7
        h = expand_gbraid(gA, vA, variant; rot = rot)
        h === nothing && continue
        r = contract_circular_cluster(h, collect(1:length(h.nodes)))
        r === nothing && continue
        (hh, vh) = r
        (euler(hh) == 2 && isempty(check_wiring(hh)) &&
         circular_canonical_key(hh) == circular_canonical_key(g1A)) && (hits += 1)
    end
    @test hits == 32
end

# ---- preimages via RELATIONS --------
#
# Rotatability is a claim about RELATIONS. Shown here: the four 8-armed preimages
# (`GBRAID_PREIMAGES`) are connected to each other by F2-3^-1 (unfold, ONE site)
# followed by a DIFFERENT F2-3 (fold, a DIFFERENT site).
@testset "8-armed preimages linked by F2-3^-1/F2-3" begin
    "The bare 8-armed :gbraid (leaves in reverse slot order, house convention)."
    function _bare8(s::Int, t::Int)
        a = [isodd(i) ? s : t for i in 1:8]
        lf(i) = mod1(2 - i, 8)
        word = a[invperm([lf(i) for i in 1:8])]
        CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
                 Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:8])
    end
    N8 = _bare8(1, 2)
    PRE = [expand_gbraid(N8, 1, v; rot = 0) for v in 1:4]
    keys4 = [circular_canonical_key(g) for g in PRE]
    @test length(Set(keys4)) == 4     # pairwise distinct

    f23_sites(g) = [(t, b) for t in 1:length(g.nodes), b in 1:length(g.nodes)
                    if _circular_is_pure_trivalent(g.nodes[t]) && g.nodes[b].kind === :braid &&
                       count(e -> e.a isa NodePort && e.b isa NodePort &&
                                  ((e.a.node == t && e.b.node == b) || (e.a.node == b && e.b.node == t)),
                             g.edges) == 1]
    function _c8_sited(g)
        out = Tuple{Int,Int,Int,CircularGraph}[]
        brs = [i for (i, nd) in enumerate(g.nodes) if nd.kind === :braid]
        trs = [i for (i, nd) in enumerate(g.nodes) if _circular_is_pure_trivalent(nd)]
        for a in brs, b in brs, t in trs
            a < b || continue
            r = _fr_braid_relation_at(g, a, b, t); r === nothing && continue
            for (h, c) in pairs_of(r); push!(out, (a, b, t, hasproperty(h, :graph) ? h.graph : h)) end
        end
        return out
    end
    "ONE move: F2-3^-1 at one site, F2-3 at ANOTHER site."
    function ein_zug(g0::CircularGraph, zielkey)
        for (tv, bi) in f23_sites(g0)
            for h1 in expand_circular_braid_relation_variants(g0, tv, bi)
                h1g = hasproperty(h1, :graph) ? h1.graph : h1
                for (a, b, t, h2) in _c8_sited(h1g)
                    circular_canonical_key(h2) == zielkey && return true
                end
            end
        end
        return false
    end

    # the four edges of a 4-ring in ONE move
    @test ein_zug(PRE[1], keys4[2])
    @test ein_zug(PRE[1], keys4[3])
    @test ein_zug(PRE[2], keys4[4])
    @test ein_zug(PRE[3], keys4[4])
    # the diagonal 1-4 not in one move, but via preimage 2 in two
    @test !ein_zug(PRE[1], keys4[4])
    @test ein_zug(PRE[1], keys4[2]) && ein_zug(PRE[2], keys4[4])   # 1 -> 2 -> 4
end

