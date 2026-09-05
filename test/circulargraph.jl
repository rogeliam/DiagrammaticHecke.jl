# test/circulargraph.jl — CircularNode/CircularGraph basics, canonical keys, merge
# machinery, and the phase-2 rotation/render system (render/CircularTutte.jl).

@testset "CircularNode — kind derivation and constructor invariants (step 1)" begin
    # circular_node derives the kind from the arm colours.
    @test circular_node([2, 2, 2, 2]).kind === :mono
    # [1,3,1]: two colours, but colour 3 only once — violates the "every colour at
    # least twice" invariant (see below), so it throws instead of giving :mixed.
    @test_throws ArgumentError circular_node([1, 3, 1])
    @test circular_node([1]).kind === :mixed
    @test circular_node([1, 2, 1, 2, 1, 2]).kind === :braid
    @test circular_node([2, 3, 2, 3, 2, 3]).kind === :braid

    # every ArgumentError path, direct constructor AND circular_node.
    @test_throws ArgumentError CircularNode(:mono, Int[])
    @test_throws ArgumentError CircularNode(:mono, [2, 1])            # :mono is colour 2 only
    @test_throws ArgumentError CircularNode(:mixed, [1, 2])           # :mixed is 1/3 only
    @test_throws ArgumentError CircularNode(:braid, [1, 2, 1])        # :braid needs exactly 6 arms
    @test_throws ArgumentError CircularNode(:braid, [1, 1, 1, 1, 1, 1])  # colour pair must be {1,2}/{2,3}
    @test_throws ArgumentError CircularNode(:braid, [1, 2, 1, 2, 2, 1])  # not alternating
    @test_throws ArgumentError CircularNode(:oops, [1])               # unknown kind
    @test_throws ArgumentError circular_node(Int[])
    @test_throws ArgumentError circular_node([1, 2])                  # fits no kind

    # rotate_arms(nd, r) == nd for every r (equality holds up to rotation only).
    nd = circular_node([1, 1, 3, 3, 1])
    for r in 0:4
        @test rotate_arms(nd, r) == nd
    end
    @test arms(rotate_arms(nd, 2)) == [3, 3, 1, 1, 1]

    # circular_key/hash consistency: same key ⇒ same hash, different rotations give the
    # same key. [1,1,3,3] carries each colour twice, as the invariant demands.
    a = circular_node([1, 1, 3, 3]); b = rotate_arms(a, 1)
    @test circular_key(a) == circular_key(b)
    @test hash(a) == hash(b)
    @test a == b
    c = circular_node([1, 1, 1, 3, 3])
    @test circular_key(a) != circular_key(c)
    @test a != c

    # circular_opposite_slot: braid-like, mod1(slot + k, 2k). On the 6-armed braid this
    # is mod1(slot+3, 6).
    braidnd = circular_node([1, 2, 1, 2, 1, 2])
    @test circular_opposite_slot(braidnd, 1) == 4
    @test circular_opposite_slot(braidnd, 6) == 3
    @test_throws ErrorException circular_opposite_slot(circular_node([1, 1, 1]), 1)
    # 8 and 10 arms: k = 4 resp. 5.
    g8nd = circular_node([1, 2, 1, 2, 1, 2, 1, 2])
    @test circular_opposite_slot(g8nd, 1) == 5
    @test circular_opposite_slot(g8nd, 8) == 4
    g10nd = circular_node([2, 3, 2, 3, 2, 3, 2, 3, 2, 3])
    @test circular_opposite_slot(g10nd, 1) == 6
    @test circular_opposite_slot(g10nd, 10) == 5
    # involution and colour fidelity (the arms alternate, so odd k swaps colours and
    # even k preserves them — only the involution is checked here).
    for nd in (braidnd, g8nd, g10nd), s in 1:arm_count(nd)
        @test circular_opposite_slot(nd, circular_opposite_slot(nd, s)) == s
    end
end

@testset "is_wired (step 1)" begin
    # positive: braid(1,3;m=2) → CircularGraph, every slot occupied.
    fg = circular(braid(1, 3; m = 2))
    @test is_wired(fg)

    # `is_wired` reads only `g.nodes`/`g.edges`, not `g.cells` — the negative
    # fixtures below violate the wired invariant (doubly occupied slot / colour
    # error / missing slot), which the generic face tracer canNOT take (it requires
    # full wiring as a precondition). They are therefore built via the 4-argument
    # constructor with a placeholder `Cells`; the 3-argument one calls the tracer.
    dummy_cells = DiagrammaticHecke.Cells(0, Int[], Vector{Int}[], Dict{Port,NTuple{2,Int}}())
    nd = circular_node([1, 1, 1])

    # negative: one slot occupied twice.
    bad = CircularGraph(CircularWord([1, 1]), [nd],
        Edge[Edge(1, Leaf(1), NodePort(1, 1)), Edge(1, Leaf(2), NodePort(1, 1))], dummy_cells)
    @test !is_wired(bad)

    # negative: edge colour does not match the arm colour.
    bad2 = CircularGraph(CircularWord([1, 1, 1]), [nd],
        Edge[Edge(1, Leaf(1), NodePort(1, 1)), Edge(1, Leaf(2), NodePort(1, 2)),
             Edge(3, Leaf(3), NodePort(1, 3))], dummy_cells)
    @test !is_wired(bad2)

    # negative: one slot not occupied at all.
    bad3 = CircularGraph(CircularWord([1, 1]), [nd],
        Edge[Edge(1, Leaf(1), NodePort(1, 1)), Edge(1, Leaf(2), NodePort(1, 2))], dummy_cells)
    @test !is_wired(bad3)
end

@testset "circular(g) — round trip WordGraph → CircularGraph (step 1)" begin
    for g in (dot(1), trivalent(2), braid(1, 2; m = 3), braid(2, 3; m = 3))
        fg = circular(g)
        @test face_count(fg) == face_count(g)
        @test fg.cells.gap_cell == g.cells.gap_cell
        @test fg.cells.sector_cell == g.cells.sector_cell
        @test fg.word == g.word
    end

    # m=2 braid(1,3) → 4-armed :mixed, NOT :braid (the core point of step 1).
    fg13 = circular(braid(1, 3; m = 2))
    @test length(fg13.nodes) == 1
    @test fg13.nodes[1].kind === :mixed
    @test arms(fg13.nodes[1]) == [1, 3, 1, 3]

    # m=3 braids stay :braid with 6 arms.
    fg12 = circular(braid(1, 2; m = 3))
    @test fg12.nodes[1].kind === :braid
    @test arms(fg12.nodes[1]) == [1, 2, 1, 2, 1, 2]

    # an invalid m=2 combination (not {1,3}) is rejected loudly.
    N = DiagrammaticHecke.Node
    @test_throws ArgumentError circular(WordGraph(CircularWord([1, 2, 1, 2]),
        [N(:braid, [1, 2], 2)],
        Edge[Edge(1, Leaf(1), NodePort(1, 1)), Edge(2, Leaf(2), NodePort(1, 2)),
             Edge(1, Leaf(3), NodePort(1, 3)), Edge(2, Leaf(4), NodePort(1, 4))]))
end

@testset "merge_nodes — the asymmetric splice assertion" begin
    # [1,1,3] slot 3 (colour 3) + [3,1,3,3] slot 1 (colour 3) → [1,1,1,3,3],
    # compared SLOT-EXACTLY as a vector. The `111` block is rotation-symmetric, so an
    # off-by-one in the splice looks there like a mere rotation (`[1,1,1] ==
    # rotate_arms([1,1,1], r)` for every r) and does NOT show up in a `CircularNode`
    # `==` comparison, so the check compares `arms(...)` as a vector.
    #
    # [1,1,3]/[3,1,3,3] themselves violate the "every colour at least twice"
    # invariant (only one 3 resp. one 1) — invalid as STANDALONE nodes, but exactly
    # the intermediate state merge_nodes is about to splice (see the other
    # `_unchecked_circularnode` places in this file). The RESULT, by contrast, must be
    # valid: `circular_node` throws otherwise. The fourth 3-arm on `nd_j` is what
    # makes the result valid (`[1,1,1,3]` is again a two-coloured node with only
    # one 3), without touching the `111` point.
    NP = DiagrammaticHecke.NodePort
    nd_i = DiagrammaticHecke._unchecked_circularnode(:mixed, [1, 1, 3])      # slots 1,2,3 = 1,1,3
    nd_j = DiagrammaticHecke._unchecked_circularnode(:mixed, [3, 1, 3, 3])   # slots 1..4 = 3,1,3,3
    g = CircularGraph(CircularWord([1, 1, 1, 1, 3, 3, 3]),
        [nd_i, nd_j],
        Edge[
            Edge(3, NP(1, 3), NP(2, 1)),                 # the edge to be contracted
            Edge(1, Leaf(1), NP(1, 1)),
            Edge(1, Leaf(2), NP(1, 2)),
            Edge(1, Leaf(3), NP(2, 2)),
            Edge(3, Leaf(4), NP(2, 3)),
            Edge(3, Leaf(5), NP(2, 4)),
        ])
    merged = merge_nodes(g, 1, 2)
    @test length(merged.nodes) == 1
    @test arms(merged.nodes[1]) == [1, 1, 1, 3, 3]        # SLOT-EXACT, not merely up to rotation
    @test merged.nodes[1].kind === :mixed
end

@testset "what slot rotation is NOT" begin
    # For [1,1,1] the canonical rotation is the identity — NOTHING is normalised
    # away. Which LEAF hangs on which slot is genuine information: two CircularGraphs with
    # the same node kind [1,1,1] but different leaf assignments to the slots are
    # DIFFERENT diagrams.
    nd = circular_node([1, 1, 1])
    g_a = CircularGraph(CircularWord([1, 1, 1]), [nd],
        Edge[Edge(1, Leaf(1), NodePort(1, 1)), Edge(1, Leaf(2), NodePort(1, 2)),
             Edge(1, Leaf(3), NodePort(1, 3))])
    # same node kind, but leaf 1 and leaf 2 on swapped slots.
    g_b = CircularGraph(CircularWord([1, 1, 1]), [nd],
        Edge[Edge(1, Leaf(2), NodePort(1, 1)), Edge(1, Leaf(1), NodePort(1, 2)),
             Edge(1, Leaf(3), NodePort(1, 3))])
    # rotate_arms alone sees no difference (EXACTLY the mistake this test pins
    # down): the node kind is the same in both graphs, and for [1,1,1] the rotation
    # canonicalisation normalises nothing away.
    @test nd == rotate_arms(nd, 1)          # node kind: equal up to rotation
    @test arms(nd) == arms(rotate_arms(nd, 1))  # for [1,1,1] even slot-exactly equal
    # BUT: the DIAGRAMS differ — which leaf hangs on which slot is genuine
    # information the canonical key MUST see.
    @test circular_canonical_key(g_a) != circular_canonical_key(g_b)
    @test g_a != g_b
end

@testset "circular_canonical_key is index-independent" begin
    fg = circular(braid(1, 3; m = 2))
    @test circular_canonical_key(fg) == circular_canonical_key(fg)

    # two trivalents joined by an edge — nodes permuted, the key must stay equal
    # (the permutation is a pure renaming of array indices, not a structural
    # change).
    c = 1
    nd1 = circular_node([c, c, c]); nd2 = circular_node([c, c, c])
    w = CircularWord([c, c, c, c])
    edges_orig = Edge[
        Edge(c, Leaf(1), NodePort(1, 1)), Edge(c, Leaf(2), NodePort(1, 2)),
        Edge(c, Leaf(3), NodePort(2, 1)), Edge(c, Leaf(4), NodePort(2, 2)),
        Edge(c, NodePort(1, 3), NodePort(2, 3)),
    ]
    g_orig = CircularGraph(w, [nd1, nd2], edges_orig)
    swap(p::NodePort) = NodePort(p.node == 1 ? 2 : p.node == 2 ? 1 : p.node, p.slot)
    swap(p::Leaf) = p
    edges_swapped = [Edge(e.colour, swap(e.a), swap(e.b)) for e in edges_orig]
    g_swapped = CircularGraph(w, [nd2, nd1], edges_swapped)   # node order swapped
    @test circular_canonical_key(g_orig) == circular_canonical_key(g_swapped)
    @test g_orig == g_swapped
end

@testset "merge_all — fixed point (step 1)" begin
    # two trivalents of the same colour joined by an edge — merge_all must fuse them
    # into ONE [1,1,1,1] node (the associative 111→1 situation).
    c = 1
    nd1 = circular_node([c, c, c]); nd2 = circular_node([c, c, c])
    w = CircularWord([c, c, c, c])
    edges = Edge[
        Edge(c, Leaf(1), NodePort(1, 1)), Edge(c, Leaf(2), NodePort(1, 2)),
        Edge(c, Leaf(3), NodePort(2, 1)), Edge(c, Leaf(4), NodePort(2, 2)),
        Edge(c, NodePort(1, 3), NodePort(2, 3)),
    ]
    g = CircularGraph(w, [nd1, nd2], edges)
    m = merge_all(g)
    @test length(m.nodes) == 1
    @test arm_count(m.nodes[1]) == 4
    @test all(==(1), arms(m.nodes[1]))

    # already merged ⇒ fixed point (merge_all(merge_all(g)) == merge_all(g)).
    @test merge_all(m) == m

    # :braid nodes are NEVER merged, not even with a connected mergeable neighbour
    # (merge_nodes rejects :braid explicitly).
    fgb = circular(braid(1, 2; m = 3))
    @test merge_all(fgb) == fgb   # nothing to do, the only node is :braid
end

# ---------------------------------------------------------------------------
# Step 2 — render/CircularTutte.jl: the phase-2 rotation system (departure angles).
# ---------------------------------------------------------------------------

# fixtures, shared by several blocks.
_circular_fixtures() = (
    circular(dot(1)),
    circular(trivalent(1)),
    circular(braid(1, 2; m = 3)),
    circular(braid(1, 3; m = 2)),
    merge_all(let
        c = 1
        nd1 = circular_node([c, c, c]); nd2 = circular_node([c, c, c])
        w = CircularWord([c, c, c, c])
        edges = Edge[
            Edge(c, Leaf(1), NodePort(1, 1)), Edge(c, Leaf(2), NodePort(1, 2)),
            Edge(c, Leaf(3), NodePort(2, 1)), Edge(c, Leaf(4), NodePort(2, 2)),
            Edge(c, NodePort(1, 3), NodePort(2, 3)),
        ]
        CircularGraph(w, [nd1, nd2], edges)
    end),
)

# the two connected m=3 braids from test/diagrams.jl — full angle degeneration
# (three legs of one node all pointing at the same target), converted by `circular`.
function _circular_g23()
    N = DiagrammaticHecke.Node; E = DiagrammaticHecke.Edge; LP = DiagrammaticHecke.Leaf; NP = DiagrammaticHecke.NodePort
    sc(nd, s) = DiagrammaticHecke._slot_colour(nd, s)
    nodes23 = [N(:braid, [2, 3], 3), N(:braid, [2, 3], 3)]
    edges23 = E[]
    for (s, lf) in [(1, 1), (2, 2), (3, 3)]
        push!(edges23, E(sc(nodes23[1], s), LP(lf), NP(1, s)))
    end
    for s in [4, 5, 6]
        push!(edges23, E(sc(nodes23[1], s), NP(1, s), NP(2, s)))
    end
    for (s, lf) in [(1, 4), (2, 5), (3, 6)]
        push!(edges23, E(sc(nodes23[2], s), LP(lf), NP(2, s)))
    end
    g23 = WordGraph(CircularWord([2, 3, 2, 2, 3, 2]), nodes23, edges23)
    return circular(g23)
end

@testset "monotonicity of the departure angles (step 2, central)" begin
    # for every inner node of all fixtures: cyclically through the slots, each
    # forward gap > 0, sum = 2π ± 1e-9. This rules out several legs sharing an
    # identical raw angle.
    for fg in (_circular_fixtures()..., _circular_g23())
        for v in 1:length(fg.nodes)
            d = arm_count(fg.nodes[v])
            leafxy, nodexy = circular_tutte_positions(fg)
            θ = DiagrammaticHecke._circular_leg_angles(fg, leafxy, nodexy)[v]
            gaps = Float64[]
            for s in 1:d
                s2 = mod1(s + 1, d)
                gap = θ[s2] - θ[s]
                s == d && (gap += 2pi)
                push!(gaps, gap)
            end
            @test all(g -> g > 0, gaps)
            @test isapprox(sum(gaps), 2pi; atol = 1e-9)
        end
    end
end

@testset "full angle degeneration (step 2)" begin
    # a node whose legs ALL point at a single point (three legs of a braid node all
    # land at a distant node that itself sits close enough to ONE point) still gets
    # d distinct angles.
    fg = circular(braid(1, 2; m = 3))
    nodexy = Dict(1 => (0.0, 0.0))
    leafxy = Dict(k => (1.0, 1.0) for k in 1:6)     # all leaves at the SAME point
    θ = DiagrammaticHecke._circular_leg_angles(fg, leafxy, nodexy)[1]
    @test length(unique(round.(θ, digits = 8))) == 6
end

# `_circular_g23()` joins slots 4→4, 5→5, 6→6 of the two braid nodes, so the strands
# cross: check_wiring reports 6 violations, euler(g23) = -4 and the embedding is not
# planar. Planarity is checked instead by the testsets "the embedding is planar:
# V − E + F = 2" and "check_wiring — the convention computed directly".

@testset "SVG surface (step 2)" begin
    fg = circular(braid(1, 2; m = 3))
    svg = circular_tutte_svg(fg)
    @test occursin("<svg", svg)

    # number of legs in the SVG (<path ... C ...>, no <circle> edges) = number of
    # non-circle edges.
    npaths = length(collect(eachmatch(r"<path d=\"M [^\"]*\" stroke=\"#\w+\" stroke-width=\"3.4\" fill=\"none\"/>", svg)))
    nondrawable = count(e -> e.a isa Circle || e.b isa Circle, fg.edges)
    @test npaths == length(fg.edges) - nondrawable

    # free circles show up.
    fgloop = CircularGraph(CircularWord(Int[]), CircularNode[], Edge[Edge(1, Circle(1), Circle(1))])
    svgloop = circular_tutte_svg(fgloop)
    @test occursin("<circle", svgloop)

    # an empty CircularGraph gives the circle-only output (early return).
    empty_fg = CircularGraph(CircularWord(Int[]), CircularNode[], Edge[])
    svgempty = circular_tutte_svg(empty_fg)
    @test occursin("<svg", svgempty)
    @test !occursin("<path", svgempty)
    @test !occursin("<text", svgempty)
end

@testset "display contract (step 2)" begin
    fg = circular(braid(1, 3; m = 2))
    @test showable("image/svg+xml", fg)
    @test showable("image/png", fg)
    @test occursin("<svg", repr("image/svg+xml", fg))
    @test length(repr("image/png", fg)) > 100
end

@testset "colour rule is total (step 2)" begin
    # every constructible arm set gives a `#` string without error.
    @test startswith(DiagrammaticHecke._circular_node_fill(circular_node([2, 2, 2, 2])), "#")   # :mono {2}
    @test startswith(DiagrammaticHecke._circular_node_fill(circular_node([1, 1, 1])), "#")      # :mixed {1}
    @test startswith(DiagrammaticHecke._circular_node_fill(circular_node([3, 3, 3])), "#")      # :mixed {3}
    @test startswith(DiagrammaticHecke._circular_node_fill(circular_node([1, 3, 1, 3])), "#")   # :mixed {1,3}
    @test startswith(DiagrammaticHecke._circular_node_fill(circular_node([1, 2, 1, 2, 1, 2])), "#")  # :braid {1,2}
    @test startswith(DiagrammaticHecke._circular_node_fill(circular_node([2, 3, 2, 3, 2, 3])), "#")  # :braid {2,3}
end

# ---------------------------------------------------------------------------
# Step 3 — the rule files (circular/CircularGraph.jl arm_count/circular_degree, circular/CircularCombo.jl,
# circular/CircularDecorated.jl, circular/CircularRules.jl).
# ---------------------------------------------------------------------------

@testset "arm_count / circular_degree — the degree formula" begin
    # arm_count is the ARM COUNT, circular_degree is the DEGREE.
    @test arm_count(circular_node([1, 1, 1])) == 3
    @test arm_count(circular_node([1])) == 1

    # degree formula: 2·(#colours) − #arms, :braid fixed at 0. Must reproduce
    # `degree(::WordGraph)` (dot_count − trivalent_count) on EVERY existing node
    # type — exactly the test that backs test/pairing.jl (Soergel grading against the
    # Hecke ground truth).
    @test circular_degree(circular_node([1])) == 1                    # Dot: 2·1−1 = +1
    @test circular_degree(circular_node([2])) == 1
    @test circular_degree(circular_node([1, 1, 1])) == -1              # Trivalent: 2·1−3 = −1
    @test circular_degree(circular_node([2, 2, 2])) == -1
    @test circular_degree(circular_node([1, 3, 1, 3])) == 0            # m=2 crossing: 2·2−4 = 0
    @test circular_degree(circular_node([1, 2, 1, 2, 1, 2])) == 0      # m=3 braid: fixed 0
    @test circular_degree(circular_node([2, 3, 2, 3, 2, 3])) == 0
    @test circular_degree(circular_node([1, 1, 1, 1])) == -2           # merged 111→1: 2·1−4
    @test circular_degree(circular_node([1, 3, 3, 1, 3, 3])) == -2     # 2·2−6

    # against degree(::WordGraph) on every fixture node type.
    for (wg_node, wg) in ((dot(1), dot(1)), (trivalent(2), trivalent(2)),
                          (braid(1, 3; m = 2), braid(1, 3; m = 2)),
                          (braid(1, 2; m = 3), braid(1, 2; m = 3)))
        fg = circular(wg)
        @test circular_degree(fg) == degree(wg)
    end

    # circular_degree(g::CircularGraph) sums over all nodes.
    nd1 = circular_node([1]); nd2 = circular_node([1, 1, 1])
    NP = DiagrammaticHecke.NodePort
    g = CircularGraph(CircularWord([1, 1]), [nd1, nd2],
        Edge[Edge(1, NP(1, 1), NP(2, 1)), Edge(1, Leaf(1), NP(2, 2)), Edge(1, Leaf(2), NP(2, 3))])
    @test circular_degree(g) == circular_degree(nd1) + circular_degree(nd2)   # 1 + (−1) = 0
end

@testset "circular_degree of a 2k-armed braid node — the formula 6 − #arms" begin
    # The tower family `π_k` has degree `6 − #arms = −(2k−6)`.
    g8 = circular_node([1, 2, 1, 2, 1, 2, 1, 2])
    @test g8.kind === :braid && arm_count(g8) == 8
    @test circular_degree(g8) == -2      # both formulas: 6−8 == 2−8/2 == −2

    # `6 − #arms` holds for EVERY braid-like node, the 6-armed braid included
    # (`6 − 6 == 0`).
    b6 = circular_node([1, 2, 1, 2, 1, 2])
    @test circular_degree(b6) == 0 == 6 - arm_count(b6)
    @test circular_degree(circular_node([2, 3, 2, 3, 2, 3])) == 0
    for n in (6, 8, 10, 12, 14)
        @test circular_degree(circular_node([isodd(i) ? 2 : 3 for i in 1:n])) == 6 - n
    end

    g10 = circular_node([1, 2, 1, 2, 1, 2, 1, 2, 1, 2])
    @test g10.kind === :braid && arm_count(g10) == 10
    @test circular_degree(g10) == -4     # double tower π₅ (2−10/2 = −3 would be odd)

    g12 = circular_node([1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2])
    @test g12.kind === :braid && arm_count(g12) == 12
    @test circular_degree(g12) == -6     # π₆

    # The parity bound on which `2 − #arms/2` fails: every generator has
    # `arms ≡ degree (mod 2)`, and every inner edge consumes two arms.
    for n in (6, 8, 10, 12, 14)
        @test iseven(circular_degree(circular_node([isodd(i) ? 1 : 2 for i in 1:n])) - n)
    end

    # C13's fusion balance works out with this formula: (6−2k)+(6−2l)−1−1
    # == 6 − 2(k+l−2).
    for k in 3:6, l in 3:6
        gk = circular_node([isodd(i) ? 1 : 2 for i in 1:2k])
        gl = circular_node([isodd(i) ? 1 : 2 for i in 1:2l])
        gm = circular_node([isodd(i) ? 1 : 2 for i in 1:(2 * (k + l - 2))])
        @test circular_degree(gk) + circular_degree(gl) - 2 == circular_degree(gm)
    end
end

@testset "merge_at_edge — bigon contraction, self-loop" begin
    NP = DiagrammaticHecke.NodePort
    c = 1
    nd1 = circular_node([c, c, c]); nd2 = circular_node([c, c, c])
    edges_bigon = Edge[
        Edge(c, Leaf(1), NP(1, 1)),
        Edge(c, Leaf(2), NP(2, 1)),
        Edge(c, NP(1, 2), NP(2, 2)),
        Edge(c, NP(1, 3), NP(2, 3)),
    ]
    gb = CircularGraph(CircularWord([c, c]), [nd1, nd2], edges_bigon)
    @test is_wired(gb)
    m = merge_at_edge(gb, 3)
    @test length(m.nodes) == 1
    @test arm_count(m.nodes[1]) == 4
    @test all(==(c), arms(m.nodes[1]))
    @test is_wired(m)
    # the surviving edge is a self-loop; the two boundary-leaf
    # edges stay ordinary leaf edges.
    e_surv = only(e for e in m.edges if e.a isa NodePort && e.b isa NodePort)
    @test e_surv.a.node == e_surv.b.node

    # merge_nodes is a thin wrapper: on multiple edges (the bigon itself) it throws,
    # and one uses the edge index directly via merge_at_edge.
    @test_throws ArgumentError merge_nodes(gb, 1, 2)

    # merge_at_edge on a non-NodePort/NodePort edge throws.
    ndx = circular_node([c])
    gx = CircularGraph(CircularWord([c]), [ndx], Edge[Edge(c, Leaf(1), NP(1, 1))])
    @test_throws ArgumentError merge_at_edge(gx, 1)

    # merge_at_edge on a :braid endpoint throws (braid(1,2;m=3), slot 1 (colour 1),
    # joined to a trivalent(1) at one of its external legs — the braid keeps 5
    # external slots, the trivalent 2, giving 7 leaves).
    fgb2 = circular(braid(1, 2; m = 3))
    braid_plus_triv = CircularGraph(CircularWord([1, 2, 1, 2, 1, 1, 1]),
        [fgb2.nodes[1], circular_node([1, 1, 1])],
        Edge[
            Edge(1, NP(1, 1), NP(2, 1)),          # braid slot 1 -- trivalent leg 1
            Edge(2, Leaf(1), NP(1, 2)),
            Edge(1, Leaf(2), NP(1, 3)),
            Edge(2, Leaf(3), NP(1, 4)),
            Edge(1, Leaf(4), NP(1, 5)),
            Edge(2, Leaf(5), NP(1, 6)),
            Edge(1, Leaf(6), NP(2, 2)),
            Edge(1, Leaf(7), NP(2, 3)),
        ])
    @test is_wired(braid_plus_triv)
    @test_throws ArgumentError merge_at_edge(braid_plus_triv, 1)
end

@testset "merge_at_edges — EXACTLY TWO neighbouring connections" begin
    # Builds A (arms armsA, n arms) and B (arms armsB, m arms), joined at TWO
    # cyclically neighbouring slot pairs: A slot k <-> B slot l+1, A slot k+1 <-> B
    # slot l (COUNTER-ORIENTED — the only planar orientation, see the CircularGraph.jl
    # header above merge_at_edges). All remaining arms go to their own leaves so it
    # can be checked slot by slot which old neighbour moves where.
    function build_two_adjacent(armsA, k, armsB, l)
        n, m = length(armsA), length(armsB)
        ndA = DiagrammaticHecke._unchecked_circularnode(:mixed, armsA)
        ndB = DiagrammaticHecke._unchecked_circularnode(:mixed, armsB)
        k2 = mod1(k + 1, n); l2 = mod1(l + 1, m)
        edges = Edge[
            Edge(armsA[k], NodePort(1, k), NodePort(2, l2)),
            Edge(armsA[k2], NodePort(1, k2), NodePort(2, l)),
        ]
        leaf = 0
        leaf_slot = Dict{Int,Tuple{Symbol,Int}}()
        for s in 1:n
            (s == k || s == k2) && continue
            leaf += 1
            push!(edges, Edge(armsA[s], Leaf(leaf), NodePort(1, s)))
            leaf_slot[leaf] = (:A, s)
        end
        for s in 1:m
            (s == l || s == l2) && continue
            leaf += 1
            push!(edges, Edge(armsB[s], Leaf(leaf), NodePort(2, s)))
            leaf_slot[leaf] = (:B, s)
        end
        g = CircularGraph(CircularWord(fill(1, max(leaf, 1))), [ndA, ndB], edges)
        return g, leaf_slot
    end

    # ---- main case: the same colour (1,1) at both connections (the case C10 does
    # NOT cover — Set((1,3)) fails) — n=m=5, result ONE node [3,3,3,3,3,3], no
    # self-loop, no bead left behind.
    armsA = [1, 1, 3, 3, 3]; armsB = [1, 1, 3, 3, 3]
    g, leaf_slot = build_two_adjacent(armsA, 1, armsB, 1)
    @test is_wired(g)
    m1 = DiagrammaticHecke.merge_at_edges(g, 1, 2)
    @test length(m1.nodes) == 1
    @test arm_count(m1.nodes[1]) == 6
    @test arms(m1.nodes[1]) == [3, 3, 3, 3, 3, 3]      # SLOT-EXACT (rotation-symmetric, but trivial here)
    @test is_wired(m1)
    @test isempty([e for e in m1.edges if e.a isa NodePort && e.b isa NodePort])  # no self-connection

    # Wiring: A from slot 3 forwards (slots 3,4,5), then B from slot 3 forwards
    # (slots 3,4,5) — each new slot must carry the old leaf that hung at exactly that
    # old (node, slot).
    expected_src = [(:A, 3), (:A, 4), (:A, 5), (:B, 3), (:B, 4), (:B, 5)]
    for (slot, src) in enumerate(expected_src)
        expected_leaf = only(k for (k, v) in leaf_slot if v == src)
        far = only(q for e in m1.edges for (p, q) in ((e.a, e.b), (e.b, e.a))
                     if p isa NodePort && p.node == 1 && p.slot == slot)
        @test far isa Leaf && far.k == expected_leaf
    end

    # ---- different colours at the two connections (overlaps C10's remit but gives
    # the same result — only merge_at_edges is checked here, not the rule driver).
    armsA2 = [1, 3, 1, 1, 3]; armsB2 = [3, 1, 1, 3, 1]
    g2, _ = build_two_adjacent(armsA2, 1, armsB2, 1)
    @test is_wired(g2)
    m2 = DiagrammaticHecke.merge_at_edges(g2, 1, 2)
    @test length(m2.nodes) == 1
    @test arms(m2.nodes[1]) == [1, 1, 3, 1, 3, 1]
    @test is_wired(m2)

    # ---- four-armed nodes (n=m=4): result ONE node [3,3,3,3].
    armsA3 = [1, 1, 3, 3]; armsB3 = [1, 1, 3, 3]
    g3, _ = build_two_adjacent(armsA3, 1, armsB3, 1)
    m3 = DiagrammaticHecke.merge_at_edges(g3, 1, 2)
    @test length(m3.nodes) == 1
    @test arms(m3.nodes[1]) == [3, 3, 3, 3]
    @test is_wired(m3)

    # ---- negative case: a CO-ORIENTED assignment (k<->l, k+1<->l+1 instead of
    # counter-oriented) throws — that is the topologically inconsistent orientation
    # (see the CircularGraph.jl header: the rest would fall into two non-neighbouring arcs
    # instead of fusing into one node).
    function build_two_same(armsA, k, armsB, l)
        n, m = length(armsA), length(armsB)
        ndA = DiagrammaticHecke._unchecked_circularnode(:mixed, armsA)
        ndB = DiagrammaticHecke._unchecked_circularnode(:mixed, armsB)
        k2 = mod1(k + 1, n); l2 = mod1(l + 1, m)
        edges = Edge[
            Edge(armsA[k], NodePort(1, k), NodePort(2, l)),
            Edge(armsA[k2], NodePort(1, k2), NodePort(2, l2)),
        ]
        leaf = 0
        for s in 1:n
            (s == k || s == k2) && continue
            leaf += 1
            push!(edges, Edge(armsA[s], Leaf(leaf), NodePort(1, s)))
        end
        for s in 1:m
            (s == l || s == l2) && continue
            leaf += 1
            push!(edges, Edge(armsB[s], Leaf(leaf), NodePort(2, s)))
        end
        return CircularGraph(CircularWord(fill(1, max(leaf, 1))), [ndA, ndB], edges)
    end
    g_same = build_two_same(copy(armsA), 1, copy(armsB), 1)
    @test is_wired(g_same)
    @test_throws ArgumentError DiagrammaticHecke.merge_at_edges(g_same, 1, 2)

    # ---- negative case: the two slots on A are NOT neighbours (slots 1 and 3 at
    # n=5) — throws.
    ndA_far = circular_node([1, 3, 1, 3, 3])
    ndB_far = circular_node([1, 3, 1, 3, 3])
    g_far = CircularGraph(CircularWord([3, 3, 3, 3, 3, 3]), [ndA_far, ndB_far],
        Edge[Edge(1, NodePort(1, 1), NodePort(2, 3)), Edge(1, NodePort(1, 3), NodePort(2, 1)),
             Edge(3, Leaf(1), NodePort(1, 2)), Edge(3, Leaf(2), NodePort(1, 4)), Edge(3, Leaf(3), NodePort(1, 5)),
             Edge(3, Leaf(4), NodePort(2, 2)), Edge(3, Leaf(5), NodePort(2, 4)), Edge(3, Leaf(6), NodePort(2, 5))])
    @test_throws ArgumentError DiagrammaticHecke.merge_at_edges(g_far, 1, 2)

    # ---- negative case: identical edge indices.
    @test_throws ArgumentError DiagrammaticHecke.merge_at_edges(g, 1, 1)

    # ---- negative case: :braid nodes are not mergeable. The :braid check in
    # merge_at_edges comes BEFORE the colour-fidelity check, so the fixture need not
    # be is_wired-correct itself; it only needs two structurally neighbouring
    # NodePort edges between the two nodes (word length = leaf count, as the generic
    # cell tracing in the 3-arg constructor demands).
    fgb = circular(braid(1, 2; m = 3))
    braid_pair = CircularGraph(CircularWord(fill(1, 6)),
        [fgb.nodes[1], DiagrammaticHecke._unchecked_circularnode(:mixed, [1, 1, 3, 3])],
        Edge[Edge(1, NodePort(1, 1), NodePort(2, 1)), Edge(2, NodePort(1, 2), NodePort(2, 2)),
             Edge(2, Leaf(1), NodePort(1, 3)), Edge(1, Leaf(2), NodePort(1, 4)),
             Edge(2, Leaf(3), NodePort(1, 5)), Edge(1, Leaf(4), NodePort(1, 6)),
             Edge(3, Leaf(5), NodePort(2, 3)), Edge(3, Leaf(6), NodePort(2, 4))])
    @test_throws ArgumentError DiagrammaticHecke.merge_at_edges(braid_pair, 1, 2)
end

@testset "CircularCombo{T} — arithmetic, term merging, zero coefficients (step 3, B2)" begin
    g1 = circular(dot(1))
    g2 = circular(dot(2))
    c1 = CircularComboR(g1)
    c2 = CircularComboR(g2)
    @test length(c1) == 1
    @test !isempty(c1)
    @test isempty(CircularComboR())

    s = c1 + c2
    @test length(s) == 2
    @test coefficient(s, g1) == one(SoergelPoly)
    @test coefficient(s, g2) == one(SoergelPoly)

    # term merging via circular_canonical_key: adding the same graph twice doubles the
    # coefficient instead of creating a second term.
    doubled = c1 + c1
    @test length(doubled) == 1
    @test coefficient(doubled, g1) == 2 * one(SoergelPoly)

    # zero coefficients are dropped.
    zeroed = c1 - c1
    @test isempty(zeroed)
    @test coefficient(zeroed, g1) == zero(SoergelPoly)

    # α_i arithmetic as a coefficient.
    scaled = alpha(1) * c1
    @test coefficient(scaled, g1) == alpha(1)

    # Int coefficients are lifted to SoergelPoly.
    @test coefficient(3 * c1, g1) == 3 * one(SoergelPoly)

    # == compares terms AND coefficients.
    @test c1 + c2 == c2 + c1
    @test c1 != c2
end

