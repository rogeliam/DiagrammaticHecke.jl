# test/circular2parallelapply.jl — circular_splice, circular_rex_relation, and the
# fusion steps (circular_rex_fusion_step, circular_dot_fusion_step).

@testset "circular_splice — splicing into the cut curve" begin
    # The self-test of the splice: cutting a tree path and splicing the IDENTITY
    # back in must reproduce the graph exactly (canonical key) — for every
    # reachable region of the fixtures.
    function identity_splice_ok(m::CircularMorphismGraph)
        g = m.graph
        for R in 1:region_count(g)
            p = circular_path_edges(m, R)
            (p === nothing || isempty(p.edges)) && continue
            r = circular_splice(g, p.edges, circular_splice_identity(p.word))
            r === nothing && return false
            circular_canonical_key(r.graph) == circular_canonical_key(g) || return false
        end
        return true
    end

    NP, L, E = NodePort, Leaf, Edge
    g30 = CircularGraph(CircularWord([1, 2, 1, 2, 1, 2, 1]), CircularNode[circular_node([2])],
                   E[E(1, L(2), L(1)), E(2, L(3), L(7)), E(1, L(4), L(6)),
                     E(2, L(5), NP(1, 1))])
    @test identity_splice_ok(CircularMorphismGraph(g30, 1, 5))

    gU, mU = ex22_figure()
    @test identity_splice_ok(mU)

    ldir = joinpath(@__DIR__, "data", "plan36")
    if isfile(joinpath(ldir, "L5.cd"))
        for n in 4:7
            fd = load_circulardecorated(joinpath(ldir, "L$n.cd"))
            @test identity_splice_ok(CircularMorphismGraph(fd.graph, 0, 5))
        end
    end

    # malformed input throws
    mpar = CircularMorphismGraph(g30, 1, 5)
    p = circular_path_edges(mpar, findfirst(R -> begin
            q = circular_path_edges(mpar, R); q !== nothing && !isempty(q.edges)
        end, 1:region_count(g30)))
    @test_throws ArgumentError circular_splice(g30, Int[], circular_splice_identity(Int[]))
    @test_throws ArgumentError circular_splice(g30, p.edges,
                                               circular_splice_identity(vcat(p.word, 1)))

    # PINNED assignment (`sides`): pinning the one the search found reproduces
    # it; a wrong length is a caller error.
    hit = circular_splice(g30, p.edges, circular_splice_identity(p.word))
    @test hit !== nothing
    pinned = circular_splice(g30, p.edges, circular_splice_identity(p.word);
                             sides = hit.sides)
    @test pinned !== nothing && pinned.sides == hit.sides
    @test circular_canonical_key(pinned.graph) == circular_canonical_key(hit.graph)
    @test_throws ArgumentError circular_splice(g30, p.edges,
                                               circular_splice_identity(p.word);
                                               sides = fill(true, length(p.edges) + 1))
end

"""
    rex_cycle_fdm(w, i) -> CircularDecoratedMorphism

The rex-cycle figure of `(w, i)`: `w → w' → w` tensored with the `i`-strand —
the diagram the relation was derived on, and the only known figure where the
region-word criterion fires with an `:edge` transition AND a fusion partner.
"""
function rex_cycle_fdm(w::Vector{Int}, i::Int)
    wp, _ = braid_to_end_with(w, i)
    f = braid_top_to(light_leaf_up(w), wp)
    fm = circular_morphism(tensor(compose(f, flip(f)), identity_strand(i)))
    return CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
end

@testset "circular_rex_relation — derived and cached" begin
    # Naming and the two caller errors; the derivation itself is only checked
    # against the CACHE, as for `.cdb`.
    @test endswith(circular_rex_cache_path([1, 2, 1], 2), "rex-121-i2.rexrel")
    @test circular_rex_cache_path([1, 2, 1], 2; dir = "X") == joinpath("X", "rex-121-i2.rexrel")
    @test_throws ArgumentError circular_rex_relation([1, 1], 1)      # w not reduced
    @test_throws ArgumentError circular_rex_relation([1, 2], 1)      # w·i still reduced

    # The carrying i-pair is (9,10) for one braid move and (22,23) for three — the
    # corrected partner, NOT `prev_edge`. And `rhs` contains no identity term.
    for (w, i, pair, nrhs) in (([1, 2, 1], 2, (9, 10), 4),
                               ([2, 1, 2, 3, 2, 1], 2, (22, 23), 5))
        isfile(circular_rex_cache_path(w, i)) || continue
        r = circular_rex_relation(w, i)
        @test r.pair == pair
        @test length(pairs_of(r.id)) == 1
        @test length(pairs_of(r.rhs)) == nrhs
        idk = Set(circular_canonical_key(d) for (d, _) in pairs_of(r.id))
        @test isempty(intersect(Set(circular_canonical_key(d) for (d, _) in pairs_of(r.rhs)), idk))
        # version 1 is label-free — the polynomials sit in the COEFFICIENTS
        @test all(all(isone, d.region_labels) && isone(d.outer_label)
                  for (d, _) in pairs_of(r.rhs))
    end
end

@testset "circular_rex_fusion_step — the relation locally" begin
    @test CIRCULAR_REX_FUSION_ENABLED[] == true           # own switch, default ON

    # THE POSITIVE FIXTURE: the rex-cycle figure of (121, 2). The criterion hits
    # region 4 / word `121` / colour 2 / edge 10, the fusion partner is edge 9
    # (the OTHER 2-edge of that region, not `prev_edge` = 3), and the cut curve
    # is the tree path plus the transition edge.
    fdm = rex_cycle_fdm([1, 2, 1], 2)
    T = circular_unreduced_transition(fdm.m)
    @test (T.region, T.word, T.colour, T.kind, T.edge, T.case) ==
          (4, [1, 2, 1], 2, :edge, 10, :dihedral)

    if isfile(circular_rex_cache_path([1, 2, 1], 2))
        r = circular_rex_fusion_step(fdm)
        @test r !== nothing
        @test r.route === :relation
        @test r.curve == [1, 2, 3, 10] && r.partner == 9 && r.pair == (9, 10)
        @test length(pairs_of(r.combo)) == 4
        # the real safeguard of an INSERT operation ("check"):
        for (d, _) in pairs_of(r.combo)
            @test euler(d.graph) == 2
            @test isempty(check_wiring(d.graph))
            @test all(isone, d.region_labels) && isone(d.outer_label)
        end
        # version 1 is label-free: one non-trivial label and the step declines
        lab = fill(one(SoergelPoly), length(fdm.region_labels))
        lab[1] = alpha(1)
        @test circular_rex_fusion_step(
                  CircularDecoratedMorphism(fdm.m, lab)) === nothing
    end

    # THE PARTNER IS NOT A GATE. On the IDENTITY diagram `id_1212` the criterion fires
    # (`:dihedral`, region word `121` + colour 2) but there is NO second 2-edge at
    # that region — the braid move of the rex cycle `1212 → 2122 → 1212` is what
    # creates it, inside the relation. The step must fire here, with `partner`
    # reported as `nothing`. Same picture one level up at `id_1213213` (`:long`).
    for (w, i, nterms) in (([1, 2, 1, 2], 2, 4), ([1, 2, 1, 3, 2, 1, 3], 3, 5))
        fm = circular_morphism(light_leaf_up(w))
        mid = CircularMorphismGraph(fm.graph, fm.cut1, fm.cut2)
        fid = CircularDecoratedMorphism(mid, fill(one(SoergelPoly),
                                                 region_count(mid.graph)))
        Tid = circular_unreduced_transition(mid)
        @test Tid !== nothing && Tid.colour == i && Tid.kind === :edge
        adjid, _ = circular_region_adjacency(mid.graph)
        @test isempty([e for (_, e) in adjid[Tid.region]
                       if e != Tid.edge && mid.graph.edges[e].colour == i])
        if isfile(circular_rex_cache_path(w[1:end-1], i))
            rid = circular_rex_fusion_step(fid)
            @test rid !== nothing
            @test rid.route === :relation
            @test rid.partner === nothing
            @test length(pairs_of(rid.combo)) == nterms
            @test isempty(circular_degree_violations(circular_degree(mid.graph),
                                                     rid.combo))
            for (d, _) in pairs_of(rid.combo)
                @test euler(d.graph) == 2
                @test isempty(check_wiring(d.graph))
            end
        end
    end

    # ROUTE `:fusion` — the trivial rex cycle (if the cycle is trivial, i.e. wi ends
    # with w on i, then plain fusion suffices). At `id_11` the two 1-edges already lie at the
    # transition region, so the step does the plain 2parallel surgery — and must
    # deliver EXACTLY what `circular_2parallel_step` delivers there.
    m11 = let fm = circular_morphism(light_leaf_up([1, 1]))
        CircularMorphismGraph(fm.graph, fm.cut1, fm.cut2)
    end
    f11 = CircularDecoratedMorphism(m11, fill(one(SoergelPoly), region_count(m11.graph)))
    r11 = circular_rex_fusion_step(f11)
    @test r11 !== nothing
    @test r11.route === :fusion
    @test r11.transition.case === :direct
    @test r11.pair == (1, 2) && r11.sides === nothing
    @test r11.combo == circular_2parallel_step(f11)      # same surgery, same terms
    @test isempty(circular_degree_violations(circular_degree(m11.graph), r11.combo))

    # NEGATIVE: the window is a `:dot` transition — that is
    # CircularDotSlide's case, and the step declines, it does not raise.
    NP, L, E = NodePort, Leaf, Edge
    g30 = CircularGraph(CircularWord([1, 2, 1, 2, 1, 2, 1]), CircularNode[circular_node([2])],
                   E[E(1, L(2), L(1)), E(2, L(3), L(7)), E(1, L(4), L(6)),
                     E(2, L(5), NP(1, 1))])
    @test circular_rex_fusion_step(
              CircularDecoratedMorphism(CircularMorphismGraph(g30, 1, 5),
                                        fill(one(SoergelPoly), region_count(g30)))) === nothing

    # L4–L7. Under the saturated reading (see the box at the `:regionword` testset
    # above): the only one of the four that carries a transition is L5, and there the
    # fusion DOES fire — 2 terms via route `:fusion`. L4/L6/L7 have no transition at
    # all and decline.
    ldir = joinpath(@__DIR__, "data", "plan36")
    if isfile(joinpath(ldir, "L5.cd"))
        for (n, expected) in ((4, nothing), (5, 2), (6, nothing), (7, nothing))
            fd = load_circulardecorated(joinpath(ldir, "L$n.cd"))
            m = CircularMorphismGraph(fd.graph, 0, 5)
            rl = circular_rex_fusion_step(
                     CircularDecoratedMorphism(m, fill(one(SoergelPoly),
                                                       region_count(fd.graph))))
            if expected === nothing
                @test rl === nothing
            else
                @test rl !== nothing && rl.route === :fusion
                @test length(pairs_of(rl.combo)) == expected
                @test isempty(circular_degree_violations(circular_degree(fd.graph), rl.combo))
                for (d, _) in pairs_of(rl.combo)
                    @test euler(d.graph) == 2
                    @test isempty(check_wiring(d.graph))
                end
            end
        end
    end
end

@testset "circular_dot_fusion_step — the general D4 rule" begin
    @test CIRCULAR_DOT_FUSION_ENABLED[] == true   # no driver hook yet

    # THE POSITIVE FIXTURE: the window (`121` + 2-dot) — for the rex
    # fusion the NEGATIVE case (`:dot`), here the positive one: transition
    # word `121`, colour 2, fall `:dihedral`, and the cut curve is the tree
    # path to the dot's region extended by the dot's cap edge.
    NP, L, E = NodePort, Leaf, Edge
    g30 = CircularGraph(CircularWord([1, 2, 1, 2, 1, 2, 1]), CircularNode[circular_node([2])],
                   E[E(1, L(2), L(1)), E(2, L(3), L(7)), E(1, L(4), L(6)),
                     E(2, L(5), NP(1, 1))])
    fdm30 = CircularDecoratedMorphism(CircularMorphismGraph(g30, 1, 5),
                                      fill(one(SoergelPoly), region_count(g30)))
    T = circular_unreduced_transition(fdm30.m)
    @test (T.word, T.colour, T.kind, T.case) == ([1, 2, 1], 2, :dot, :dihedral)

    # `:edge` transitions are the rex fusion's case — the step declines there.
    @test circular_dot_fusion_step(rex_cycle_fdm([1, 2, 1], 2)) === nothing

    if isfile(circular_rex_cache_path([1, 2, 1], 2))
        r = circular_dot_fusion_step(fdm30)
        @test r !== nothing
        p = circular_path_edges(fdm30.m, T.region)
        @test r.curve == vcat(p.edges, T.edge) && r.dot == T.to
        @test length(pairs_of(r.combo)) == 4
        # the real safeguard of an INSERT operation ("check"):
        for (d, _) in pairs_of(r.combo)
            @test euler(d.graph) == 2
            @test isempty(check_wiring(d.graph))
            @test all(isone, d.region_labels) && isone(d.outer_label)
        end
        # version 1 is label-free: one non-trivial label and the step declines
        lab = fill(one(SoergelPoly), length(fdm30.region_labels))
        lab[1] = alpha(1)
        @test circular_dot_fusion_step(
                  CircularDecoratedMorphism(fdm30.m, lab)) === nothing
    end

    # `reduce = true` cleans up every term with reduce_circular_combo —
    # on the figure (`id_121321 ⊗ 3-dot`) the dots at nodes vanish
    # (one :merge/C5 step each), while the term count stays 5.
    if isfile(circular_rex_cache_path([1, 2, 1, 3, 2, 1], 3))
        fm50 = circular_morphism(tensor(light_leaf_up([1, 2, 1, 3, 2, 1]),
                                        dot_morphism(3)))
        fdm50 = CircularDecoratedMorphism(fm50, fill(one(SoergelPoly),
                                                     region_count(fm50.graph)))
        r50 = circular_dot_fusion_step(fdm50; reduce = true)
        @test r50 !== nothing && length(pairs_of(r50.combo)) == 5
        dots_on(g) = count(e -> e.a isa NodePort && e.b isa NodePort &&
                                ((arm_count(g.nodes[e.a.node]) == 1) ⊻
                                 (arm_count(g.nodes[e.b.node]) == 1)), g.edges)
        @test all(dots_on(d.graph) == 0 for (d, _) in pairs_of(r50.combo))
    end
end

