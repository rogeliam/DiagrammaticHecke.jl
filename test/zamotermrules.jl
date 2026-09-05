# test/zamotermrules.jl — the Zamo rules as whole-term replacement + sited + the
# Pipeline step (src/rules/ZamoTermRules.jl).
#
# ⚠️ COST: the first `zamo_term_rules()` call derives Z2 from the two
# X/Y reductions (~40 s warm) — once per process, cached afterwards.

using Test
using DiagrammaticHecke
const CWZ = DiagrammaticHecke

@testset "Zamo whole-term rules (ZamoTermRules.jl)" begin

    @testset "Z1 alone is cheap" begin
        t0 = time()
        r1 = zamo_term_rules(z2 = false)
        @test length(r1) == 1
        @test r1[1].name === :zamo_1_8
        @test length(r1[1].rhs) == 1 && r1[1].rhs[1][2] == 1
    end

    rules = zamo_term_rules()          # the one-off Z2 derivation happens here

    @testset "the rule inventory and the termination control" begin
        @test [R.name for R in rules] == [:zamo_1_8, :zamo_2_9]
        @test [s for (_, s) in rules[2].rhs] == [1, 1, -1]
        @test all(length(g.nodes) == 7 for (g, _) in rules[2].rhs)
        # no right-hand side carries a left-hand key again (otherwise zamo_term_rules
        # itself would throw — this is only the visual check on top)
        for R in rules, (g, _) in R.rhs
            @test !any(S.lhs == CWZ._zamo_mkeyf(g, 0) for S in rules)
        end
    end

    @testset "zamo_pass: Z1/Z2 on their own left-hand side" begin
        for (nm, m, nterms) in ((:zamo_1_8, flip(Zamo(8, 1)), 1),
                                (:zamo_2_9, flip(Zamo(9, 2)), 3))
            cL = CircularComboR(circular_decorated(circular(m.graph)))
            cR, hit_names = zamo_pass(cL, m.cut1)
            @test hit_names == [nm]
            @test length(collect(pairs_of(cR))) == nterms
            _, again = zamo_pass(cR, m.cut1)
            @test isempty(again)                 # never hits its own output
        end
    end

    @testset "sited: apply_zamo_rule_combo" begin
        # Z2 on its own left-hand side == zamo_pass (keys + coefficients)
        mZ2 = flip(Zamo(9, 2))
        a = first(m.anchor for m in find_zamo_matches(mZ2.graph; pair = (2, 9)))
        combo = apply_zamo_rule_combo(mZ2.graph, a; pair = (2, 9))
        ref, _ = zamo_pass(CircularComboR(circular_decorated(circular(mZ2.graph))), mZ2.cut1)
        ks(c) = sort([(string(CWZ._zamo_mkeyf(t.graph, mZ2.cut1)), co)
                      for (t, co) in pairs_of(c)], by = first)
        @test ks(combo) == ks(ref)
        @test all(isempty(check_wiring(t.graph)) for (t, _) in pairs_of(combo))

        # Z1 EMBEDDED (Zamo(1,1), 14 nodes): bit-identical at every anchor to the
        # plain rule apply_zamo_rule, after the circular conversion.
        Z11 = Zamo(1, 1)
        a11 = unique(m.anchor for m in find_zamo_matches(Z11.graph; pair = (1, 8)))
        @test !isempty(a11)
        for a in a11
            c = apply_zamo_rule_combo(Z11.graph, a; pair = (1, 8))
            ts = collect(pairs_of(c))
            @test length(ts) == 1 && isone(ts[1][2])
            plain = apply_zamo_rule(Z11.graph, a; pair = (1, 8))
            @test circular_canonical_key(ts[1][1].graph) == circular_canonical_key(circular(plain))
        end

        # no false alarm: Z2 does not match in Zamo(1,1) (a different boundary word)
        @test isempty(find_zamo_matches(Z11.graph; pair = (2, 9)))
        @test apply_zamo_rule_combo(Z11.graph, 1; pair = (2, 9)) === nothing
    end

    @testset "pipeline: reduce_to_circular_leave fires the Zamo step" begin
        # rev(Zamo(8,1)) is a CIRCULAR_RULES fixed point; only the Zamo step replaces it —
        # and the result is exactly Zamo(1,8).
        for (m, pairname, nterms) in ((flip(Zamo(8, 1)), :zamo_1_8, 1),
                                      (flip(Zamo(9, 2)), :zamo_2_9, 3))
            fm = circular_morphism(m)
            fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
            c = reduce_to_circular_leave(fdm)
            @test length(collect(pairs_of(c))) == nterms
        end
        # The result is the circular-leaf NORMAL FORM of Zamo(1,8) — after the Zamo step the
        # pipeline reduces further, and in Zamo(1,8) P1 merges parallel 1/3 strands into
        # [1,1,3,3] nodes. So BOTH sides are sent through the pipeline and then compared;
        # on Zamo(1,8) itself the Zamo step does not fire (an RHS key).
        leave(m) = begin
            fm = circular_morphism(m)
            fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
            (reduce_to_circular_leave(fdm), fm.cut1)
        end
        cL, cutL = leave(flip(Zamo(8, 1)))
        cR, cutR = leave(Zamo(1, 8))
        tL, coL = only(pairs_of(cL)); tR, coR = only(pairs_of(cR))
        @test isone(coL) && isone(coR)
        # ⚠️ NOT bit-identical: both sides are 9-node normal forms with two P1-merged
        # [1,1,3,3] nodes and an identical node multiset, but P1's choice of pair depends
        # on the leaf numbering — the P1 ambiguity: two normal forms of the SAME
        # morphism. What is tested are the robust invariants. Arm lists are
        # determined only up to cyclic rotation — so normalise them.
        cyc(a) = minimum([circshift(a, r) for r in 0:(length(a) - 1)])
        nk(g) = sort([(nd.kind, cyc(CWZ.arms(nd))) for nd in g.nodes], by = string)
        @test nk(tL.graph) == nk(tR.graph)
        @test length(tL.graph.nodes) == 9
        # [1,1,3,3] = 4-armed, colours {1,1,3,3}, equal colours ADJACENT
        # (the crossing [1,3,1,3] has the same multiset, but alternating)
        ist1133(nd) = arm_count(nd) == 4 && sort(CWZ.arms(nd)) == [1, 1, 3, 3] &&
                      any(CWZ.arms(nd)[i] == CWZ.arms(nd)[mod1(i + 1, 4)] for i in 1:4)
        @test count(ist1133, tL.graph.nodes) == 2
        @test circular_degree(tL.graph) == circular_degree(tR.graph)
    end
end

