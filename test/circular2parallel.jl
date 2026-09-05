# 2parallel — two parallel edges of the SAME colour (src/circular/rules/Circular2Parallel.jl).
#
# CONDITION 4 IS THE REGION WORD: the pair must be the (tree parent edge,
# transition edge) of the region carrying the shortest unreduced transition —
# the pair rex fusion fuses at `case == :direct`, so the two rules agree on the
# site. What that means for the fixtures here:
#   * the minimal figure `1111` (cut `0|2`) is a POSITIVE fixture: the region
#     word `11` is unreduced and the transition names the pair `(1, 2)` at
#     region 2;
#   * the same figure with the marker in the MIDDLE (cut `1|3`) is negative:
#     the region word is reduced everywhere, so there is no transition;
#   * `crossing_figure_13` fires at edges 2/7, region 4, colour 1;
#   * `ex22_figure` fires at edges 17/18, region 5, colour 2;
#   * on `ex22_figure` the step tracer reaches a FIXED POINT at `maxsteps = 120`.
# A marking must be set; without one there is no transition and nothing fires.

using DiagrammaticHecke: Leaf, Edge, CircularMorphismGraph, CircularDecorated,
                     CircularDecoratedMorphism, _circular_transfer_labels

"""
    ex22_figure() -> (gU, mU)

The CLEANED-UP figure: the double leaf over `21211212` without the leaves
1, 8, 9, 16 — boundary word `121121|121121`, 12 leaves. It is the only known place
where 2parallel fires at an INNER region, and hence the positive fixture of this
file.
"""
function ex22_figure()
    w  = [2, 1, 2, 1, 1, 2, 1, 2]
    D  = double_leaf(w, ones(Int, 8), w, ones(Int, 8))
    fD = circular_morphism(D)
    gD = fD.graph
    keep_slots = vcat(2:7, 10:15)
    neu    = Dict(k => i for (i, k) in enumerate(keep_slots))
    umbau(p) = p isa Leaf ? Leaf(neu[p.k]) : p
    keep = Edge[Edge(ed.colour, umbau(ed.a), umbau(ed.b))
                for (i, ed) in enumerate(gD.edges) if !(i in (1, 10, 11, 20))]
    push!(keep, Edge(2, NodePort(1, 3), NodePort(3, 3)))
    push!(keep, Edge(2, NodePort(2, 4), NodePort(4, 2)))
    gU = CircularGraph(CircularWord(letters(gD.word)[keep_slots]), gD.nodes, keep)
    return gU, CircularMorphismGraph(gU, 0, 6)
end

"""
    crossing_figure_13() -> (g, m)

The 1/3 variant: 1 and 3 lie in ONE component.

A 3-strand crosses two nested 1-strands — boundary word `[3,1,1,3,1,1]`, two
`[1,3,1,3]` crossings, nothing else. Through the two crossings and the 3-edge between
them EVERYTHING hangs in ONE `circular_connected_components` component.
"""
function crossing_figure_13()
    g = CircularGraph(CircularWord([3, 1, 1, 3, 1, 1]),
                 CircularNode[circular_node([1, 3, 1, 3]), circular_node([1, 3, 1, 3])],
                 Edge[Edge(1, Leaf(2), NodePort(1, 1)),
                      Edge(1, NodePort(1, 3), Leaf(6)),
                      Edge(3, Leaf(1), NodePort(1, 2)),
                      Edge(3, NodePort(1, 4), NodePort(2, 2)),
                      Edge(3, NodePort(2, 4), Leaf(4)),
                      Edge(1, Leaf(3), NodePort(2, 1)),
                      Edge(1, NodePort(2, 3), Leaf(5))])
    return g, CircularMorphismGraph(g, 0, 3)
end

@testset "2parallel — parallel edges of equal colour" begin
    # ---- THE MINIMAL FIGURE: boundary word 1111, two parallel 1-strands ------
    # As a morphism (cut 0|2) this is the IDENTITY on `11`.
    par = CircularGraph(CircularWord([1, 1, 1, 1]), CircularNode[],
                   Edge[Edge(1, Leaf(1), Leaf(4)), Edge(1, Leaf(2), Leaf(3))])
    m   = CircularMorphismGraph(par, 0, 2)

    # The three regions sit in the distance pattern (i−1, i, i+1) = (0, 1, 2) — the
    # shape of the figure, not a condition.
    @test circular_region_distances(m) == [0, 1, 2]

    # ---- THE IDENTITY ON 1,1 FIRES -------------------------------------------
    # It fires although the two edges do not sit in one component: the two
    # leaf–leaf strands are two separate strands, and `:d4` asks nothing about
    # where the edges live — only that one side is closer to the marking.
    hm = find_circular_2parallel(m)
    @test hm !== nothing
    @test (hm.ei, hm.ej, hm.region, hm.colour) == (1, 2, 2, 1)
    cm = circular_2parallel_step(CircularDecoratedMorphism(
             m, fill(one(SoergelPoly), region_count(par))))
    @test cm !== nothing && length(pairs_of(cm)) == 2
    @test all(circular_degree(d) == circular_degree(par) for (d, _) in pairs_of(cm))
    @test sort([string(x) for (d, _) in pairs_of(cm)
                for x in d.region_labels if !isone(x)]) ==
          [string((1//2) * alpha(1)), string((1//2) * alpha(1))]
    # The dropped conditions really are gone: the middle region is NOT interior,
    # and both edges end at leaves — neither blocks the hit.
    bw = circular_region_boundary_word(par, 2)
    @test length(bw) == 4 && !is_interior(bw)

    # The surgery: ONE 4-armed 1-node, planar and violation-free, and the shared
    # region splits into TWO. (`circular_2parallel_apply` checks NO condition, so the
    # figure stays usable as a surgery fixture.)
    f = circular_2parallel_apply(par, 1, 2)
    @test f !== nothing
    @test f.graph.nodes[f.node].arms == [1, 1, 1, 1]
    @test euler(f.graph) == 2
    @test isempty(check_wiring(f.graph))
    @test region_count(f.graph) == region_count(par) + 1
    @test f.newregs[1] != f.newregs[2]

    # THE MARKER IS A GENUINE INPUT TO THE RULE. The same figure with the marker in
    # the MIDDLE (cut 1|3) is not the identity as a morphism but cup–cap; the boundary
    # circle is cut open DIFFERENTLY there, and the two leaf sets sit as SEPARATE
    # blocks (`A A B B`) instead of nested, and the distances are (1, 0, 1). So: no
    # hit.
    m_middle = CircularMorphismGraph(par, 1, 3)
    @test circular_region_distances(m_middle) == [1, 0, 1]
    @test find_circular_2parallel(m_middle) === nothing

    # WITHOUT A MARKER (cut 0|0). `find_circular_2parallel` carries on even when all
    # distances are `-1`, and `cut1 = 0` is the canonical gap between leaf n and 1 —
    # the same situation as above, hence the same hit.
    @test find_circular_2parallel(CircularMorphismGraph(par, 0, 0)) !== nothing

    # ---- the 1/3 variant -----------------------------------------------------
    # 1 and 3 lie in ONE component: a 3-strand crosses two nested 1-strands.
    gK, mK = crossing_figure_13()
    @test euler(gK) == 2 && isempty(check_wiring(gK))
    # The figure FIRES: the transition names region 4 and the pair (2, 7). The
    # distances are recorded because the returned `far_lo`/`far_hi` are ordered
    # by them, not because any condition reads them.
    @test circular_region_distances(mK) == [1, 2, 0, 1, 2, 3]
    hK = find_circular_2parallel(mK)
    @test hK !== nothing
    @test (hK.ei, hK.ej, hK.region, hK.colour) == (2, 7, 4, 1)
    # ---- the FIGURE: the hit ----------------------------------------
    gU, mU = ex22_figure()
    @test bottom(mU) == [1, 2, 1, 1, 2, 1] && top(mU) == [1, 2, 1, 1, 2, 1]
    @test euler(gU) == 2 && isempty(check_wiring(gU))
    @test isempty(reduce_circular(gU)[2])            # still the CIRCULAR_RULES normal form

    hU = find_circular_2parallel(mU)
    @test hU !== nothing

    # WHICH PAIR IS HIT: the one the region word names. The transition sits at
    # region 5 and gives the pair (17, 18) in colour 2 — the same pair rex
    # fusion fuses there. This is the only known place where 2parallel fires at
    # an INNER region, which is what makes this figure the positive fixture.
    T = circular_unreduced_transition(mU)
    @test T !== nothing && T.case === :direct && T.kind === :edge
    @test T.region == 5 && minmax(T.prev_edge, T.edge) == (17, 18)
    @test hU.region == 5
    @test (hU.ei, hU.ej) == (17, 18)
    @test gU.edges[hU.ei].colour == 2 && gU.edges[hU.ej].colour == 2
    bwU = circular_region_boundary_word(gU, hU.region)
    @test is_interior(bwU)

    cU = circular_2parallel_step(CircularDecoratedMorphism(mU,
             fill(one(SoergelPoly), region_count(gU))))
    @test cU !== nothing
    @test length(pairs_of(cU)) == 2
    @test all(circular_degree(d) == circular_degree(gU) for (d, _) in pairs_of(cU))
    @test all(euler(d.graph) == 2 && isempty(check_wiring(d.graph))
              for (d, _) in pairs_of(cU))
    # one 4-armed MONOCHROME node and one α/2 per term, in the colour of the
    # fused pair; the degree check is unchanged.
    @test all(count(n -> arm_count(n) == 4, d.graph.nodes) == 1
              for (d, _) in pairs_of(cU))
    @test all(count(x -> !isone(x), d.region_labels) == 1 for (d, _) in pairs_of(cU))
    @test sort([string(x) for (d, _) in pairs_of(cU)
                for x in d.region_labels if !isone(x)]) ==
          [string((1//2) * alpha(2)), string((1//2) * alpha(2))]

    # condition 5: the rule does not act across a region carrying a polynomial.
    labsU = fill(one(SoergelPoly), region_count(gU))
    labsU[hU.region] = alpha(1)
    @test circular_2parallel_step(CircularDecoratedMorphism(mU, labsU)) === nothing

    # ---- THE TERMINATION MEASURE (condition 6) ------------------------------
    # The rule's header claims: the number of `s`-edges hanging on NO monochrome
    # `s`-node with ≥ 3 arms falls by 2 per application. The rule CHAINS here, which
    # is what makes the claim checkable: here it is one application, the measure
    # falls 18 → 16, and then the figure is a fixed point. The proof only requires
    # the strict fall in steps of two — how many steps a given figure allows is a
    # property of the figure, not of the proof.
    _2p_mass(g) = count(e -> !DiagrammaticHecke._circular_edge_at_general_2parallel(
                                 g, e, e.colour), g.edges)
    masse = Int[_2p_mass(gU)]
    gk = gU
    for _ in 1:20
        h = find_circular_2parallel(CircularMorphismGraph(gk, mU.cut1, mU.cut2))
        h === nothing && break
        fk = circular_2parallel_apply(gk, h.ei, h.ej)
        fk === nothing && break
        gk = fk.graph
        push!(masse, _2p_mass(gk))
    end
    @test masse == [18, 16]
    @test all(masse[k] > masse[k + 1] for k in 1:length(masse) - 1)
    @test find_circular_2parallel(CircularMorphismGraph(gk, mU.cut1, mU.cut2)) === nothing

    # ---- THE EDGE CASE FIRES REGULARLY -------------------------
    # An edge ending at a leaf is no obstacle: the rule fires without a bypass.
    @test find_circular_2parallel(mU) !== nothing
end

# ---- The step tracer (src/render/CircularSteps.jl) -----------------------------
# The picture equation `D -> c1*D1 + c2*D2 + ...` as general machinery. Checked
# here: `CIRCULAR_STEPS` finds 2parallel, and the tracer runs into a fixed point
# without violations.
#
# The positive fixture is the figure (`ex22_figure`). The minimal figure `1111` is
# a positive fixture too, so the step tracer finds a step on it as well; the
# negative claim holds only for the BASIC RULES (`circular_steps(:circular_rules)`).
@testset "step tracer CIRCULAR_STEPS / circular_trace" begin
    gU, mU = ex22_figure()

    st = circular_step_once(circular_decorated(gU), mU)
    @test st !== nothing
    @test st.name === Symbol("2parallel")
    @test length(pairs_of(st.after)) == 2

    par = CircularGraph(CircularWord([1, 1, 1, 1]), CircularNode[],
                   Edge[Edge(1, Leaf(1), Leaf(4)), Edge(1, Leaf(2), Leaf(3))])
    m   = CircularMorphismGraph(par, 0, 2)
    # NO BASIC RULE applies to the minimal figure — it consists of just two
    # leaf–leaf edges; there is nothing to merge, to dot, or to turn.
    @test circular_step_once(circular_decorated(par), m; steps = circular_steps(:circular_rules)) === nothing
    # The full step set does find 2parallel, though: two terms, each with `α₁/2` in
    # one of the two new regions.
    st_min = circular_step_once(circular_decorated(par), m)
    @test st_min !== nothing
    @test st_min.name === Symbol("2parallel")
    @test length(pairs_of(st_min.after)) == 2
    @test sort([string(x) for (d, _) in pairs_of(st_min.after)
                for x in d.region_labels if !isone(x)]) ==
          [string((1//2) * alpha(1)), string((1//2) * alpha(1))]
    @test first.(circular_steps(:d4, :fusion)) == [:fusion, :d4]   # the order of CIRCULAR_STEPS
    @test_throws ErrorException circular_steps(:gibtsnicht)

    # The full tracer: every intermediate state is planar and violation-free (which
    # `circular_trace` checks itself and counts).
    #
    # Under the `:d4` default the tracer reaches a FIXED POINT: it stops after 77
    # steps at 8 terms (node counts 3…6), well within `maxsteps = 120`. The step
    # sequence contains EXACTLY ONE 2parallel step, and after it only basic rules,
    # fusion and scalars.
    endc = circular_trace(gU, mU; draw = false, maxsteps = 120)
    @test !isempty(pairs_of(endc))
    @test all(euler(d.graph) == 2 && isempty(check_wiring(d.graph))
              for (d, _) in pairs_of(endc))
    @test all(circular_step_once(d, mU) === nothing for (d, _) in pairs_of(endc))
end

# ---- condition 4: the region word ------------------------------------
# The default stays `:d4`.
@testset "2parallel — condition 4 is the region word" begin

    # (a) case-1 fixture: the minimal figure. Transition (region 2, word [1],
    #     colour 1, :edge, :direct); :regionword hits the SAME pair as :d4.
    par = CircularGraph(CircularWord([1, 1, 1, 1]), CircularNode[],
                   Edge[Edge(1, Leaf(1), Leaf(4)), Edge(1, Leaf(2), Leaf(3))])
    m = CircularMorphismGraph(par, 0, 2)
    T = circular_unreduced_transition(m)
    @test T !== nothing
    @test (T.region, T.word, T.colour, T.kind, T.case) == (2, [1], 1, :edge, :direct)
    hm = find_circular_2parallel(m)
    @test hm !== nothing && (hm.ei, hm.ej, hm.region, hm.colour) == (1, 2, 2, 1)

    # marker in the middle: no transition (the middle region is the closest), no hit
    # — as under every other variant.
    m_middle = CircularMorphismGraph(par, 1, 3)
    @test circular_unreduced_transition(m_middle) === nothing
    @test find_circular_2parallel(m_middle) === nothing

    # (b) the two other fixtures: :regionword picks DIFFERENT places than :d4.
    gK, mK = crossing_figure_13()
    hK = find_circular_2parallel(mK)
    @test hK !== nothing && (hK.ei, hK.ej, hK.region, hK.colour) == (2, 7, 4, 1)
    gU, mU = ex22_figure()
    hU = find_circular_2parallel(mU)
    @test hU !== nothing && (hU.ei, hU.ej, hU.region, hU.colour) == (17, 18, 5, 2)


    # (c) the id window (id₁₂₁ ⊗ 2-dot, marker 1|5): the case-2 fixture.
    #     Transition = dot transition (word 121, colour 2, :dihedral) — and the move is
    #     the DOT-SLIDE rule (the dot region carries NO 2-edge, so the "forced
    #     D4" is constructively empty there; see the docstring of
    #     circular_regionword_dihedral_step).
    NP, L, E = NodePort, Leaf, Edge
    g30 = CircularGraph(CircularWord([1, 2, 1, 2, 1, 2, 1]), CircularNode[circular_node([2])],
                   E[E(1, L(2), L(1)), E(2, L(3), L(7)), E(1, L(4), L(6)),
                     E(2, L(5), NP(1, 1))])
    m30 = CircularMorphismGraph(g30, 1, 5)
    T30 = circular_unreduced_transition(m30)
    @test T30 !== nothing
    @test (T30.region, T30.word, T30.colour, T30.kind, T30.case) ==
          (4, [1, 2, 1], 2, :dot, :dihedral)
    fdm30 = CircularDecoratedMorphism(m30, fill(one(SoergelPoly), region_count(g30)))
    r30 = circular_regionword_dihedral_step(fdm30)
    @test r30 !== nothing && r30.dot == 1
    @test r30.combo == circular_dot_slide_step(fdm30)          # the move IS the slide
    @test length(pairs_of(r30.combo)) == 4
    # under :regionword 2parallel does NOT fire there (no :direct transition):
    @test find_circular_2parallel(m30) === nothing

    # (d) L4..L7 (marker 0|5 — the marker is part of the input).
    #     The SATURATED reading (`circular_region_words_2k`, `CIRCULAR_REGIONWORD_2K`,
    #     default on): around one node the word never grows past that node's own
    #     longest element — 1 letter at a single-colour node, 2 at `{1,3}`, 3 at
    #     `{1,2}`/`{2,3}`.
    #     Under that reading L4 reads `2,2,2,ε,2` (the target value, see
    #     `test/circularmvalentwords.jl`) and has NO transition; L5 keeps its
    #     `:direct` transition; L6 and L7 are quiet. So exactly ONE of the four
    #     fires. The strict guard still rejects every pair.
    ldir = joinpath(@__DIR__, "data", "plan36")
    if isfile(joinpath(ldir, "L5.cd"))
        for (n, expected) in ((4, nothing),
                          (5, (3, [2], 2, :edge, :direct)),
                          (6, nothing),
                          (7, nothing))
            fd = load_circulardecorated(joinpath(ldir, "L$n.cd"))
            mm = CircularMorphismGraph(fd.graph, 0, 5)
            Tn = circular_unreduced_transition(mm)
            if expected === nothing
                @test Tn === nothing
            else
                @test Tn !== nothing
                @test (Tn.region, Tn.word, Tn.colour, Tn.kind, Tn.case) == expected
            end
            @test find_circular_2parallel(mm) === nothing
        end
    end
end

# ---- the region word AS A LIST OF EDGES ----------------------
# Path and region word are the SAME object — the colours
# of the tree path from the mark to R are exactly `circular_region_words(m)[1][R]`.
@testset "circular_path_edges — path and region word are the same object" begin
    # (a) the minimal figure (two nested 1-edges, mark 0|2): region 2 has word [1].
    par = CircularGraph(CircularWord([1, 1, 1, 1]), CircularNode[],
                   Edge[Edge(1, Leaf(1), Leaf(4)), Edge(1, Leaf(2), Leaf(3))])
    m = CircularMorphismGraph(par, 0, 2)
    words, parent = circular_region_words(m)
    p2 = circular_path_edges(m, 2)
    @test p2 !== nothing
    @test p2.word == words[2] == [1]
    @test p2.edges == [parent[2]]
    @test p2.regions[end] == 2

    # the marked region itself: the empty path (and its word is empty).
    start = findfirst(w -> w !== nothing && isempty(w), words)
    @test start !== nothing
    p0 = circular_path_edges(m, start)
    @test p0 !== nothing && isempty(p0.edges) && isempty(p0.word) && p0.regions == [start]

    # (b) the probe over EVERY reachable region of the bigger fixtures: colours of
    #     the path == region word, path starts at the mark, ends at R, and each
    #     step is one level further out.
    for (g, mm) in (crossing_figure_13(), ex22_figure())
        ws, _ = circular_region_words(mm)
        st = findfirst(w -> w !== nothing && isempty(w), ws)
        for R in 1:region_count(mm.graph)
            if ws[R] === nothing
                @test circular_path_edges(mm, R) === nothing
                continue
            end
            p = circular_path_edges(mm, R)
            @test p !== nothing
            @test p.word == ws[R]
            @test length(p.edges) == length(ws[R]) == length(p.regions) - 1
            @test p.regions[1] == st && p.regions[end] == R
            @test allunique(p.edges)
            for (j, Rj) in enumerate(p.regions)
                @test length(ws[Rj]::Vector{Int}) == j - 1
            end
        end
    end

    # (c) out of range / no mark: nothing, no exception.
    @test circular_path_edges(m, 0) === nothing
    @test circular_path_edges(m, region_count(par) + 1) === nothing
end

@testset "circular_path_edges — the route as an edge list" begin
    # THE CHECK OF THE STEP: path and region word are the SAME object —
    # `colours(edges) == word` and `word == circular_region_words(m)[1][R]`, for
    # EVERY reachable region. Fixtures: the figure (the only known inner-region
    # case), the id window, and L4..L7 (marker 0|5).
    function check_paths(m::CircularMorphismGraph)
        g = m.graph
        words, _ = circular_region_words(m)
        for R in 1:region_count(g)
            p = circular_path_edges(m, R)
            if words[R] === nothing
                @test p === nothing
            else
                @test p !== nothing
                @test [g.edges[e].colour for e in p.edges] == p.word
                @test p.word == words[R]
                @test length(p.regions) == length(p.edges) + 1
                @test p.regions[end] == R
                # the path starts at the marked region (its word is empty)
                @test words[p.regions[1]] == Int[]
            end
        end
    end

    gU, mU = ex22_figure()
    check_paths(mU)
    # the marked region itself: the empty path
    words_U, _ = circular_region_words(mU)
    start = findfirst(w -> w == Int[], words_U)
    pU = circular_path_edges(mU, start)
    @test pU !== nothing && isempty(pU.edges) && pU.regions == [start]
    # out of range
    @test circular_path_edges(mU, 0) === nothing
    @test circular_path_edges(mU, region_count(gU) + 1) === nothing

    NP, L, E = NodePort, Leaf, Edge
    g30 = CircularGraph(CircularWord([1, 2, 1, 2, 1, 2, 1]), CircularNode[circular_node([2])],
                   E[E(1, L(2), L(1)), E(2, L(3), L(7)), E(1, L(4), L(6)),
                     E(2, L(5), NP(1, 1))])
    check_paths(CircularMorphismGraph(g30, 1, 5))

    ldir = joinpath(@__DIR__, "data", "plan36")
    if isfile(joinpath(ldir, "L5.cd"))
        for n in 4:7
            fd = load_circulardecorated(joinpath(ldir, "L$n.cd"))
            check_paths(CircularMorphismGraph(fd.graph, 0, 5))
        end
    end
end

