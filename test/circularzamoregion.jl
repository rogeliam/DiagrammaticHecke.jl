# test/circularzamoregion.jl — Zamo as a MULTI-REGION pattern
# (src/circular/rules/CircularZamoRegion.jl). Both Zamos are implemented in general.
#
# ⚠️ COST: the (2,9) SURGERY goes through `apply_zamo_rule_combo` and derives Z2 once
# from the two X/Y reductions (~40 s warm, as in test/zamotermrules.jl). The MATCHER
# itself is cheap — it does not touch `zamo_term_rules()`. Hence the matcher tests come
# first.

using Test
using DiagrammaticHecke

@testset "Zamo as a multi-region pattern (CircularZamoRegion.jl)" begin

    # ---- M1: the census that unblocked the point -------------------------------
    # There are six inner regions, not seven. All FOUR Zamo sides, not only the
    # two left ones.
    @testset "M1 — six inner regions, all of boundary-word length 3" begin
        for (nm, m) in (("flip(Zamo(8,1))", flip(Zamo(8, 1))),
                        ("Zamo(1,8)",       Zamo(1, 8)),
                        ("flip(Zamo(9,2))", flip(Zamo(9, 2))),
                        ("Zamo(2,9)",       Zamo(2, 9)))
            g = circular(m.graph)
            inner = [w for w in circular_boundary_words(g)
                     if !isempty(w.letters) && is_interior(w)]
            @test length(g.nodes) == 7
            @test length(inner) == 6                       # NOT 7
            @test all(length(w.letters) == 3 for w in inner)
            # the six regions determine the cluster completely
            corners = unique(vcat([[l.corner.vertex[2] for l in w.letters]
                                 for w in inner]...))
            @test sort(corners) == collect(1:7)
        end
    end

    # ---- M2: do the two Zamos have the SAME pattern? ---------------------------
    # The first measurement: NO. (1,8) is a 4-cycle with two pendants, (2,9) a
    # hexagon.
    @testset "M2 — (1,8) and (2,9) have DIFFERENT patterns" begin
        P18 = zamo_region_pattern((1, 8))
        P29 = zamo_region_pattern((2, 9))
        @test P18.len == P29.len == 3
        @test length(P18.colours) == length(P29.colours) == 6
        @test length(P18.nodes) == length(P29.nodes) == 7

        s18 = zamo_region_signature(flip(Zamo(8, 1)))
        s29 = zamo_region_signature(flip(Zamo(9, 2)))
        @test s18 != s29
        # (2,9) is a HEXAGON: six edges between inner regions, each
        # region has exactly two inner neighbours. (1,8) likewise has six
        # edges, but two regions with only ONE inner neighbour (the pendants).
        function inner_degrees(m)
            g = circular(m.graph)
            iset = Set(w.region for w in circular_boundary_words(g)
                       if !isempty(w.letters) && is_interior(w))
            adj, _ = circular_region_adjacency(g)
            return sort([count(x -> x[1] in iset, adj[R]) for R in iset])
        end
        @test inner_degrees(flip(Zamo(9, 2))) == [2, 2, 2, 2, 2, 2]     # hexagon
        @test inner_degrees(flip(Zamo(8, 1))) == [1, 1, 2, 2, 3, 3]     # cycle + pendants
        # the edge colours separate the two hexagons: LHS {2,3}, RHS {1,2}
        @test sort(unique(s29.edge_colours)) == [2, 3]
        @test sort(unique(zamo_region_signature(Zamo(2, 9)).edge_colours)) == [1, 2]
    end

    # ---- M3: the pattern finds itself in its own fixture, exactly once ---------
    @testset "M3 — self-match with the right anchor" begin
        for p in ((1, 8), (2, 9))
            L = flip(Zamo(p[2], p[1]))
            ms = circular_zamo_region_matches(L; pair = p)
            @test length(ms) == 1
            plain = find_zamo_matches(L.graph; pair = p)
            @test !isempty(plain)
            # the anchor candidates of the region version cover the plain anchors
            @test Set(unique(m.anchor for m in plain)) ⊆ Set(ms[1].anchors)
            @test length(ms[1].regions) == 6
            @test sort(collect(values(ms[1].nodemap))) == collect(1:7)
        end
    end

    # ---- M4: freedom from self-matches (the termination argument) -----------------
    # No pattern may sit in a right-hand side — otherwise the rule would hit its own
    # output. For (1,8) that is not trivial: LHS and RHS have the SAME shape, and only
    # the edge colour separates them.
    @testset "M4 — no pattern in a right-hand side" begin
        for p in ((1, 8), (2, 9)), R in (Zamo(1, 8), Zamo(2, 9))
            @test isempty(circular_zamo_region_matches(R; pair = p))
        end
    end

    # ---- M5: the oracle measurement at Zamo(1,1) -------------------------------
    # The whole 14-step cycle: 17 inner regions of length 3 (the bare
    # COUNT 6 is thus only a pre-filter), plain finds three (1,8) sites
    # and no (2,9) site.
    @testset "M5 — region matcher against find_zamo_matches on Zamo(1,1)" begin
        Z = Zamo(1, 1)
        inner3 = count(w -> !isempty(w.letters) && is_interior(w) &&
                            length(w.letters) == 3, circular_boundary_words(circular(Z.graph)))
        @test inner3 == 17                     # the 6 alone would be no guard

        for p in ((1, 8), (2, 9))
            ms = circular_zamo_region_matches(Z; pair = p)
            plain = find_zamo_matches(Z.graph; pair = p)
            @test length(ms) == length(plain)
            @test Set(sort(collect(values(m.nodemap))) for m in ms) ==
                  Set(sort(collect(values(m.phi))) for m in plain)
        end
        @test length(circular_zamo_region_matches(Z; pair = (1, 8))) == 3
        @test isempty(circular_zamo_region_matches(Z; pair = (2, 9)))
    end

    # ---- M6: the right-hand side against the plain oracle ----------------------
    # This is where the Z2 derivation is incurred (~40 s the first time).
    keys_of(c) = sort([string(circular_canonical_key(t.graph), "|", co)
                       for (t, co) in pairs_of(c)])

    @testset "M6 — right-hand side == apply_zamo_rule_combo" begin
        Z = Zamo(1, 1)
        for m in circular_zamo_region_matches(Z; pair = (1, 8))
            o = apply_zamo_rule_combo(Z.graph, m.anchors[1]; pair = (1, 8))
            @test o !== nothing
            @test length(pairs_of(o)) == 1                 # Z1: one term
        end
        r = circular_zamo_region_step(Z.graph; pair = (1, 8))
        @test r !== nothing
        @test keys_of(r) == keys_of(apply_zamo_rule_combo(
            Z.graph, circular_zamo_region_matches(Z; pair = (1, 8))[1].anchors[1];
            pair = (1, 8)))

        # (2,9) on its own left-hand side — three terms, +1/+1/−1
        L = flip(Zamo(9, 2))
        ms = circular_zamo_region_matches(L; pair = (2, 9))
        r2 = circular_zamo_region_step(L.graph; pair = (2, 9))
        o2 = apply_zamo_rule_combo(L.graph, ms[1].anchors[1]; pair = (2, 9))
        @test r2 !== nothing && o2 !== nothing
        @test length(pairs_of(r2)) == 3
        @test keys_of(r2) == keys_of(o2)

        # no site: nothing
        @test circular_zamo_region_step(Zamo(1, 2).graph; pair = (1, 8)) === nothing
    end

    # ---- M7: termination in the PIPELINE ---------------------------------------
    # M4 checks the fixtures; M7 checks what the driver actually computes: on the END
    # TERMS of `reduce_to_circular_leave` the region version must not match. ⚠ Without
    # the rotation guard `colour_ok`, P1 and the region version can oscillate
    # (9→8→9→9→…) into a StackOverflow — this affects 26 of the 196 pairs
    # `Zamo(i,j)`. `Zamo(4,3)` (13 end terms, likewise 0 self-matches) is too
    # expensive for the suite at ~154 s.
    @testset "M7 — no self-match on the end terms of Zamo(1,7)" begin
        fm = circular_morphism(Zamo(1, 7))
        fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
        c = reduce_to_circular_leave(fdm; maxrounds = 100, node_growth_limit = 15)
        term_list = first.(collect(pairs_of(c)))
        @test length(term_list) == 1                       # one end term, 9 nodes
        @test length(term_list[1].graph.nodes) == 9
        for dd in term_list, p in ((1, 8), (2, 9))
            @test circular_zamo_region_step(dd.graph; pair = p) === nothing
        end
    end

    # ---- the morphism version ---------------------------------------------------
    @testset "the MorphismGraph version forwards the cuts" begin
        Z = Zamo(1, 1)
        res = circular_zamo_region_step(Z; pair = (1, 8))
        @test res !== nothing
        c, c1, c2 = res
        @test (c1, c2) == (Z.cut1, Z.cut2)
        @test length(pairs_of(c)) == 1
    end

    # ==== Zamo in general ================

    # ---- Z1: Zamo triangle ≠ Zamo region ---------------------------------------
    # The fixtures stay at 6 (the definition cuts away nothing healthy); zamo_triangles
    # over-counts the Zamo(1,8) pipeline end term (8 vs zamo_regions' 6), Zamo(1,7)
    # (8 vs 5), and Zamo(1,10) (13 vs 10) — table in the block comment of
    # CircularZamoRegion.jl.
    @testset "Z1 — zamo_triangles / zamo_regions on the fixtures" begin
        for m in (flip(Zamo(8, 1)), Zamo(1, 8), flip(Zamo(9, 2)), Zamo(2, 9))
            g = circular(m.graph)
            @test length(zamo_triangles(g)) == 6
            @test length(zamo_regions(g)) == 6
            @test length(zamo_region_components(g)) == 1
        end
        # A Zamo triangle at a [1,1,3,3] node (P1-merged) does not count: the Zamo(1,7)
        # end term (test M7 below measures it too).
        @test length(zamo_regions(Zamo(1, 2))) == 0
    end

    # ---- Z2: the backward patterns ----------------------------------------------
    # The RHS pattern exactly once in Zamo(1,8) resp. Zamo(2,9), nowhere in the LHS.
    @testset "Z2 — backward patterns Z1⁻¹/Z2⁻¹" begin
        for p in ((1, 8), (2, 9))
            Pi = zamo_region_pattern(p; inverse = true)
            @test Pi.len == 3 && length(Pi.colours) == 6 && length(Pi.nodes) == 7
            R = Zamo(p[1], p[2])
            @test length(circular_zamo_region_matches(R; pair = p, inverse = true)) == 1
            @test isempty(circular_zamo_region_matches(flip(Zamo(p[2], p[1])); pair = p,
                                                  inverse = true))
        end
        @test isempty(circular_zamo_region_matches(Zamo(2, 9); pair = (1, 8), inverse = true))
        @test isempty(circular_zamo_region_matches(Zamo(1, 8); pair = (2, 9), inverse = true))
        # The terms of the backward direction, computed from the forward rule: Z1⁻¹ one
        # term (+1), Z2⁻¹ three terms with signs +1, −1, +1.
        _, t18 = zamo_inverse_terms((1, 8))
        @test [s for (_, s) in t18] == [1]
        _, t29 = zamo_inverse_terms((2, 9))
        @test length(t29) == 3 && t29[1][2] == 1 && sort([s for (_, s) in t29]) == [-1, 1, 1]
        # Z1⁻¹ on Zamo(1,8) gives exactly rev(Zamo(8,1)) — and forwards, the Z1 pattern
        # is found on it again (the backward direction is the inverse).
        inv = circular_zamo_region_step(circular(Zamo(1, 8).graph); pair = (1, 8), inverse = true)
        @test inv !== nothing && length(pairs_of(inv)) == 1
        back = first(pairs_of(inv))[1].graph
        @test length(circular_zamo_region_matches(back; pair = (1, 8))) == 1
        # morphism key (up to rotation of the leaves, as in ZamoTermRules.jl)
        L = flip(Zamo(8, 1))
        @test DiagrammaticHecke._zamo_mkeyf(back, Zamo(1, 8).cut1) ==
              DiagrammaticHecke._zamo_mkeyf(circular(L.graph), L.cut1)
    end

    # ---- Z3: the choice of direction in the pipeline ---------------------------
    # Zamo(1,7): 5 Zamo regions (< 7) ⇒ forwards only, no backward attempt;
    # Zamo(1,10): 10 (≥ 7), forwards blocked ⇒ Z1⁻¹ applies, then CIRCULAR_RULES — result:
    # 6 terms with ≤ 2 Zamo regions.
    function _end_terms(i, j)
        fm = circular_morphism(Zamo(i, j))
        fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
        return reduce_to_circular_leave(fdm; maxrounds = 100, node_growth_limit = 15)
    end
    @testset "Z3 — direction: Zamo(1,7) forwards only, Zamo(1,10) backwards" begin
        zamo_inverse_log_clear!()
        c7 = _end_terms(1, 7)
        @test isempty(zamo_inverse_log())
        t7 = first.(collect(pairs_of(c7)))
        @test length(t7) == 1 && length(zamo_regions(t7[1].graph)) == 5

        zamo_inverse_log_clear!()
        c10 = _end_terms(1, 10)
        log = zamo_inverse_log()
        @test length(log) == 1 && log[1].status === :accepted
        @test log[1].zamo_regions == 10 && log[1].fires == [:circular_rules]
        t10 = first.(collect(pairs_of(c10)))
        @test length(t10) == 6
        @test all(length(zamo_regions(dd.graph)) <= 2 for dd in t10)
        # What the step is for: the 10 Zamo regions are gone (the line above),
        # and the figure comes apart into few nodes. The node counts are pinned
        # as a multiset so that a change in the rule set is visible here rather
        # than hidden behind a loose bound.
        @test sort([length(dd.graph.nodes) for dd in t10]) == [4, 5, 5, 5, 8, 8]
        zamo_inverse_log_clear!()
    end

    # ---- Z4: the guard -------------------------------------------------------
    # `circular_zamo_direction_step` directly: on a Zamo LHS the forward direction applies
    # (`inverse == false`); on `Zamo(1,8)` alone (6 regions) nothing applies —
    # and no backward direction below 7. The acceptance criterion
    # `zamo_nonzamo_rule_fires` reports `:p1` on the BARE Zamo sides (the `[1,3,1,3]`
    # nodes of the fixtures carry parallel 1/3 edges, and P1 is in
    # the pipeline BEFORE Zamo — which is why `Zamo(1,8)` as input becomes a 9-node end
    # term); on an end term of the pipeline it reports `nothing`.
    @testset "Z4 — guard: no backward direction below 7 regions, acceptance criterion" begin
        L = flip(Zamo(8, 1))
        r = circular_zamo_direction_step(circular(L.graph), L.cut1, L.cut2; log = false)
        @test r !== nothing && r[2] == false
        R = Zamo(1, 8)
        @test circular_zamo_direction_step(circular(R.graph), R.cut1, R.cut2; log = false) === nothing
        @test zamo_nonzamo_rule_fires(circular_decorated(circular(R.graph)), R.cut1, R.cut2) === :p1
        c7 = _end_terms(1, 7)
        dd = first(pairs_of(c7))[1]
        fm = circular_morphism(Zamo(1, 7))
        @test zamo_nonzamo_rule_fires(dd, fm.cut1, fm.cut2) === nothing
    end

    # ---- Z5: the measurement switch ------------------------------------------
    # `CIRCULAR_ZAMO_ENABLED[] = false` turns BOTH Zamo versions off: the
    # The step sequence of `reduce_to_circular_leave(circular_morphism(zamo_lhs()))` then contains
    # no `:zamo_*` step, and the reduction still terminates. With the switch on (the
    # default) a Zamo step appears — unchanged.
    @testset "Z5 — measurement switch CIRCULAR_ZAMO_ENABLED" begin
        _phases(fm) = begin
            ph = Symbol[]
            fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly),
                                                region_count(fm.graph)))
            c = reduce_to_circular_leave(fdm; maxrounds = 100, node_growth_limit = 15,
                                    trace = (d, p, r, co) -> push!(ph, p))
            (c, ph)
        end
        @test CIRCULAR_ZAMO_ENABLED[]                       # default: ON
        # zamo_rhs() is the fixture whose reduction really does take a Zamo step
        # (zamo_lhs() runs without one: only :d4_round/:after_reduce_circular).
        fm = circular_morphism(zamo_rhs())
        c_on, ph_on = _phases(fm)
        @test any(startswith(String(p), "zamo") for p in ph_on)
        CIRCULAR_ZAMO_ENABLED[] = false
        try
            c_off, ph_off = _phases(fm)
            @test !any(startswith(String(p), "zamo") for p in ph_off)
            @test !isempty(collect(pairs_of(c_off)))   # terminates Zamo-free
        finally
            CIRCULAR_ZAMO_ENABLED[] = true
        end
    end
end

