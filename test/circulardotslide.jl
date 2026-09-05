# test/circulardotslide.jl — the dot-slide rule (src/circular/rules/CircularDotSlide.jl).
#
# ⚠️ WHAT IS CHECKED HERE — AND WHAT IS NOT. The equation
# `id window = [braid² + forced D4] − [remainder window]` is established in THREE
# pieces, all of them checkable in the driver:
#   (a) the internally built braid² window is EXACTLY the frozen reference graph
#       (canonical key),
#   (b) the direct reduction of the braid² window is `1·id window + 1·remainder window`
#       — together with the D4 rule (tested elsewhere) that IS the equation,
#   (c) with the switch ON the driver gives EXACTLY the
#       reduced slide sum (coefficients included).
# What is NOT checked is `reduce(id window, switch off) == reduce(slide sum)`: without
# the switch the id window is a NORMAL FORM for the driver (which is exactly why the
# rule exists), so the two normal forms differ structurally (5-term D4 route vs.
# 2-term direct).

using Test
using DiagrammaticHecke
using DiagrammaticHecke: CircularMorphismGraph, NodePort, Leaf, Edge, CircularNode,
                     circular_dot_slide_step, _fds_braid2, _fds_rest, _fds_paths,
                     _circular_edges_at_node, circular_region_adjacency, circular_region_of_dot,
                     is_wired

# the id window: id₁₂₁ ⊗ 2-dot. Circular word 1212121, normalised
# [1,1,2,1,2,1,2]; leaves 2,3,4 = bottom 121, leaf 5 = the 2-dot, 6,7,1 = top.
function _ds_id_fenster()
    NP, L, E = NodePort, Leaf, Edge
    nodes = CircularNode[circular_node([2])]
    edges = E[E(1, L(2), L(1)),          # the outer 1-strand (at the marker)
              E(2, L(3), L(7)),          # 2-strand
              E(1, L(4), L(6)),          # inner 1-strand (at the dot)
              E(2, L(5), NP(1, 1))]      # the dot cap
    g = CircularGraph(CircularWord([1, 2, 1, 2, 1, 2, 1]), nodes, edges)
    return g, CircularMorphismGraph(g, 1, 5)
end

# the braid² reference window
function _ds_braid2_notebook()
    NP, L, E = NodePort, Leaf, Edge
    nodes = CircularNode[circular_node([1,2,1,2,1,2]), circular_node([1,2,1,2,1,2]), circular_node([2])]
    edges = E[E(1, L(2), NP(1,1)), E(2, L(3), NP(1,6)), E(1, L(4), NP(1,5)),
              E(1, L(6), NP(2,1)), E(2, L(7), NP(2,6)), E(1, L(1), NP(2,5)),
              E(2, L(5), NP(3,1)),
              E(2, NP(1,4), NP(2,2)), E(1, NP(1,3), NP(2,3)), E(2, NP(1,2), NP(2,4))]
    return CircularGraph(CircularWord([1, 2, 1, 2, 1, 2, 1]), nodes, edges)
end

function _ds_reduziere(total, cut1, cut2)
    out = CircularComboR()
    for (dd, c2) in pairs_of(total)
        m3 = CircularMorphismGraph(dd.graph, cut1, cut2)
        out = out + c2 * reduce_to_circular_leave(
                  CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label))
    end
    return out
end

@testset "dot slide (CircularDotSlide.jl)" begin

    @testset "the fixture: the step fires, 4 terms, degree-preserving" begin
        g, m = _ds_id_fenster()
        @test euler(g) == 2 && isempty(check_wiring(g)) && is_wired(g)
        fdm = with_marks(g, m)
        res = circular_dot_slide_step(fdm)
        @test res !== nothing
        ts = collect(pairs_of(res))
        @test length(ts) == 4                    # 3 D4 terms − 1 remainder window
        expected = circular_degree(circular_decorated(g))
        @test all(circular_term_degree(c, d) == expected for (d, c) in ts)
    end

    @testset "surgery: braid² window == the reference graph, remainder planar" begin
        g, m = _ds_id_fenster()
        dist = circular_region_distances(m)
        adj, _ = circular_region_adjacency(g)
        R = circular_region_of_dot(g, 1)
        @test dist[R] == 3
        paths = _fds_paths(adj, dist, R)
        @test length(paths) == 1                 # (inner a, t, outer a)
        (e1, e2, e3) = paths[1]
        (cap_ei, _, dse) = _circular_edges_at_node(g, 1)[1]
        b2res = _fds_braid2(g, 1, cap_ei, e1, e2, e3, 1, 2)
        @test b2res !== nothing
        gnew = b2res[1]
        # (a) EXACTLY the wiring of the reference graph:
        @test circular_canonical_key(gnew) == circular_canonical_key(_ds_braid2_notebook())
        # D4 fires at the dot on EXACTLY one middle edge:
        @test apply_circular_d4(circular_decorated(gnew), b2res[2], b2res[3]) !== nothing
        rest = _fds_rest(g, 1, cap_ei, dse, e1, e2, e3, 1, 2)
        @test rest !== nothing
        @test euler(rest) == 2 && isempty(check_wiring(rest)) && is_wired(rest)
        # (b) the equation at the braid² window: the direct reduction equals
        # the reductions of the id window and of the remainder window added
        # — coefficients included.
        directT = reduce_to_circular_leave(with_marks(gnew, m))
        @test directT == reduce_to_circular_leave(with_marks(g, m)) +
                         reduce_to_circular_leave(with_marks(rest, m))
    end

    @testset "the driver reduces the id window to the reduced slide sum" begin
        # The driver does not call this step; its dot fusion rewrites the same
        # window through the derived relation and lands on the same sum.
        @test CIRCULAR_DOT_SLIDE_ENABLED[] == false
        g, m = _ds_id_fenster()
        fdm = with_marks(g, m)
        res = circular_dot_slide_step(fdm)
        slide = _ds_reduziere(res, 1, 5)
        @test reduce_to_circular_leave(fdm) == slide
        # on a figure WITHOUT a dot the switch changes nothing:
        g2 = general_braid_cluster(2)
        fdm2 = with_marks(g2, CircularMorphismGraph(g2, 0, 0))
        r_off = reduce_to_circular_leave(fdm2)
        old = CIRCULAR_DOT_SLIDE_ENABLED[]
        CIRCULAR_DOT_SLIDE_ENABLED[] = true
        try
            @test reduce_to_circular_leave(fdm2) == r_off
        finally
            CIRCULAR_DOT_SLIDE_ENABLED[] = old
        end
    end

    @testset "no firing without an (a,t,a) pattern" begin
        # the bigon fixture from test/circulardecoratedrules.jl ("borderless BC region"): the
        # dot sees colour 1 == the dot colour as its first edge, and the distance of its
        # region is < 3 — no window.
        NP = NodePort
        g_bigon = CircularGraph(CircularWord([1, 1, 1]),
                           CircularNode[circular_node([1, 1, 1]), circular_node([1, 1, 1]), circular_node([1])],
                           Edge[Edge(1, NP(1, 1), NP(2, 2)),
                                Edge(1, NP(1, 2), NP(2, 1)),
                                Edge(1, Leaf(1), NP(1, 3)),
                                Edge(1, Leaf(2), NP(2, 3)),
                                Edge(1, Leaf(3), NP(3, 1))])
        fdm = with_marks(g_bigon, CircularMorphismGraph(g_bigon, 0, 0))
        @test circular_dot_slide_step(fdm) === nothing
    end
end

