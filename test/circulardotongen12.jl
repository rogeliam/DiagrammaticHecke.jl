# test/circulardotongen12.jl — C14 `dot_on_gen12` (src/circular/rules/CircularDotOnGen12.jl).
#
# The `:gbraid` is not frozen as a literal but produced via C13
# from the fixture `_g12_caseC_classic` (test/circulargen12merge.jl) — so the two
# rules stay tied together.
#
# It is checked against the GROUND TRUTH: "dot on the preimage first, then reduce
# fully" must give the same `circular_canonical_key` as "C13 first, then dot, then C14". Plus
# the two invariants of every rule (boundary word, `circular_degree`), planarity/wiring and
# the termination evidence.
#
# Why `check_wiring` must not be dropped here: the other cyclic order of the three
# trivalent arms gives the SAME node types but is not planar. A test on
# node kinds alone would let that through — exactly the mistake R8 failed on.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_dot_on_gen12, _fr_gen12_merge, Edge, NodePort, Leaf

# Fixture as in test/circulargen12merge.jl (called `_g12_caseC_classic` there).
_c14_preimage() = CircularGraph(CircularWord([2, 3, 2, 3, 2, 3, 2, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]),
            circular_node([3, 3, 3]), circular_node([3, 3, 3])],
    Edge[
        Edge(3, NodePort(1, 4), NodePort(3, 1)),
        Edge(3, NodePort(3, 3), NodePort(2, 2)),
        Edge(2, NodePort(1, 5), NodePort(2, 1)),
        Edge(3, NodePort(1, 6), NodePort(4, 1)),
        Edge(3, NodePort(4, 2), NodePort(2, 6)),
        Edge(2, Leaf(1), NodePort(1, 1)),
        Edge(3, Leaf(2), NodePort(4, 3)),
        Edge(2, Leaf(3), NodePort(2, 5)),
        Edge(3, Leaf(4), NodePort(2, 4)),
        Edge(2, Leaf(5), NodePort(2, 3)),
        Edge(3, Leaf(6), NodePort(3, 2)),
        Edge(2, Leaf(7), NodePort(1, 3)),
        Edge(3, Leaf(8), NodePort(1, 2)),
    ])

"Dot on boundary leaf `k`."
function _c14_attach_dot(g::CircularGraph, k::Int)
    lets = letters(g.word); c = lets[k]
    nodes = vcat(copy(g.nodes), CircularNode[circular_node([c])]); dotidx = length(nodes)
    shift(p) = p isa Leaf ? (p.k == k ? NodePort(dotidx, 1) : Leaf(p.k > k ? p.k - 1 : p.k)) : p
    CircularGraph(CircularWord(vcat(lets[1:k-1], lets[k+1:end])), nodes,
             Edge[Edge(e.colour, shift(e.a), shift(e.b)) for e in g.edges])
end

_c14_graph(t) = t isa CircularGraph ? t : t.graph
_c14_key(c) = sort([(string(co), circular_canonical_key(_c14_graph(t)))
                    for (t, co) in pairs_of(c)], by = string)

# ⚠ the DOT POLICY (`circular_dots_reducible`) blocks C14 and the
# generic preimage stage on a `2k >= 8`-node with ONE dot — there the dot
# stays put. This block tests the SURGERY, not the gate, and hence
# runs with the switch turned off. The gate itself has its own
# tests in `test/circulardotmerge.jl`.
#
# ⚠ the DIRECTION has to be turned as well. `CIRCULAR_TRIVALENT_MERGE` is the
# one switch for three things (circular/rules/CircularWeight.jl): `false` means
# C14 is live, C19 is dead, and `circular_weight` is the ARM weight. The
# statements of this block are formulated in exactly that world — under `true`
# the C14 step RISES in the separating-edge reading, and the driver no longer
# takes it.
_old_policy = CIRCULAR_DOT_POLICY[]
_old_direction = CIRCULAR_TRIVALENT_MERGE[]
CIRCULAR_DOT_POLICY[] = false
CIRCULAR_TRIVALENT_MERGE[] = false
try

@testset "C14 dot_on_gen12 (CircularDotOnGen12.jl)" begin

    G8 = first(pairs_of(_fr_gen12_merge(_c14_preimage())))[1]

    @testset "the :gbraid that is computed on" begin
        @test length(G8.nodes) == 1
        @test G8.nodes[1].kind === :braid
        @test arm_count(G8.nodes[1]) == 8
        @test isempty(check_wiring(G8)) && euler(G8) == 2
    end

    @testset "colour 2 (leaf 1): :braid + trivalent of colour 3" begin
        gd = _c14_attach_dot(G8, 1)
        res = _fr_dot_on_gen12(gd)
        @test res !== nothing
        ts = collect(pairs_of(res))
        @test length(ts) == 1 && ts[1][2] == 1
        h = ts[1][1]
        @test sort([(nd.kind, arm_count(nd)) for nd in h.nodes]) ==
              [(:braid, 6), (:mixed, 3)]
        @test all(arm_colour(nd, s) == 3 for nd in h.nodes if nd.kind === :mixed
                  for s in 1:arm_count(nd))
        @test join(letters(h.word)) == join(letters(gd.word))   # boundary word preserved
        @test circular_degree(h) == circular_degree(gd)                   # degree preserved
        @test isempty(check_wiring(h)) && euler(h) == 2
        @test circular_weight(h) < circular_weight(gd)                    # termination
    end

    @testset "colour 3 (leaf 2): :braid + trivalent of colour 2" begin
        gd = _c14_attach_dot(G8, 2)
        h = first(pairs_of(_fr_dot_on_gen12(gd)))[1]
        @test sort([(nd.kind, arm_count(nd)) for nd in h.nodes]) ==
              [(:braid, 6), (:mono, 3)]
        @test circular_degree(h) == circular_degree(gd)
        @test isempty(check_wiring(h)) && euler(h) == 2
    end

    @testset "all eight leaves == ground truth (dot on the preimage, fully reduced)" begin
        matches = 0
        for k in 1:8
            gd = _c14_attach_dot(G8, k)
            out = _fr_dot_on_gen12(gd)
            gt  = reduce_circular_full(_c14_attach_dot(_c14_preimage(), k))
            _c14_key(out) == _c14_key(gt) && (matches += 1)
        end
        @test matches == 8
    end

    @testset "C14 is registered — the driver takes it by itself" begin
        # C14 is registered in `CIRCULAR_RULES`, behind the dot policy, so the driver
        # reaches it directly. Within this block the policy is OFF (top of the file),
        # otherwise the dots here would stay put.
        gd = _c14_attach_dot(G8, 1)
        (c, hist) = reduce_circular(gd)
        @test :dot_on_gen12 in hist
        @test length(collect(pairs_of(c))) == 1
        # the same canonical key as the direct C14 call.
        out14 = _fr_dot_on_gen12(gd)
        @test _c14_key(c) == _c14_key(out14)
    end

    @testset "8/8 boundary leaves: driver == direct C14 call" begin
        # Assurance test: for EACH of the eight boundary leaves, reduce_circular (via the
        # generic preimage stage) gives the same circular_canonical_key as the direct
        # _fr_dot_on_gen12 call.
        matches = 0
        for k in 1:8
            gd = _c14_attach_dot(G8, k)
            (c, hist) = reduce_circular(gd)
            out14 = _fr_dot_on_gen12(gd)
            (:dot_on_gen12 in hist) && _c14_key(c) == _c14_key(out14) && (matches += 1)
        end
        @test matches == 8
    end

    @testset "the termination weight (the reason for the circular_weight change)" begin
        gd = _c14_attach_dot(G8, 1)
        h  = first(pairs_of(_fr_dot_on_gen12(gd)))[1]
        @test circular_weight(gd) == (8, 2, 9)     # braid(8) + dot(1)
        @test circular_weight(h)  == (6, 2, 9)     # braid(6) + trivalent(3)
        # The number of `:braid` nodes is not a good measure here: it stays flat at 1 on
        # both sides (`gd` and `h` each have exactly one `:braid` node) — not strictly
        # falling. What falls strictly is the ARM SUM over the braid-like nodes, 8 -> 6.
        @test count(nd -> nd.kind === :braid, gd.nodes) ==
              count(nd -> nd.kind === :braid, h.nodes) == 1
        @test circular_weight(h)[1] < circular_weight(gd)[1]
    end

    @testset "no false alarm" begin
        # without a dot: nothing to do
        @test _fr_dot_on_gen12(G8) === nothing
        # a dot at an ordinary :braid is C6, not C14
        @test _fr_dot_on_gen12(_c14_attach_dot(_c14_preimage(), 1)) === nothing
    end
end

finally
    CIRCULAR_DOT_POLICY[] = _old_policy
    CIRCULAR_TRIVALENT_MERGE[] = _old_direction
end

