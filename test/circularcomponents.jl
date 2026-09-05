# circularcomponents.jl — freely floating connected components
# (src/circular/CircularComponents.jl).

@testset "circular_connected_components" begin
    # Node 3 is a single dot of colour 2 on leaf 1; leaf 2 hangs directly on node 3's
    # only edge (a dot has just 1 arm — the edge runs from leaf 1 to the dot, and the
    # rest of the strand needs NO node of its own, see the `leaf_colour`/
    # `join_graphs` convention: a dot caps off the strand). Nodes 1/2 are the isolated
    # barbell component, hanging on no leaf.
    n1 = circular_node([1]); n2 = circular_node([1]); n3 = circular_node([2])
    g = CircularGraph(CircularWord([2]), [n1, n2, n3],
        [Edge(2, Leaf(1), NodePort(3, 1)),
         Edge(1, NodePort(1, 1), NodePort(2, 1))])
    comps = circular_connected_components(g)
    @test Set(Set.(comps)) == Set([Set([1, 2]), Set([3])])
end

# ⚠ A polynomial that belongs to NO cell any more does not become a scalar factor but
# stays in the `outer_label`. `circular_extract_floating_components` therefore returns
# `factor == 1` and writes the value outside. The tests here therefore check the PRODUCT
# `factor * outer_label` resp. `coeff * outer_label` — it is the same value either way,
# and it is exactly that value which is the mathematical statement.

@testset "circular_extract_floating_components: isolated barbell next to a strand" begin
    # A strand of colour 2 (via a trivalent-free intermediate node, so as to have
    # genuine leaf edges) PLUS an isolated barbell component (two degree-1 dots of
    # colour 1) — `_fr_barbell`'s case, via the general pre-step rather than the
    # incremental rule.
    n1 = circular_node([1]); n2 = circular_node([1])
    g = CircularGraph(CircularWord([2, 2]), [n1, n2],
        [Edge(2, Leaf(1), Leaf(2)),
         Edge(1, NodePort(1, 1), NodePort(2, 1))])
    fdm = CircularDecoratedMorphism(CircularMorphismGraph(g, 0, 2),
                                fill(one(SoergelPoly), region_count(g)))
    fnew, factor = circular_extract_floating_components(fdm)

    @test factor * fnew.outer_label == alpha(1)   # the same α as _fr_barbell
    @test length(fnew.m.graph.nodes) == 0          # barbell node gone
    @test all(isone, fnew.region_labels)           # nothing wrongly labelled
    @test bottom(fnew.m) == [2, 2] && top(fnew.m) == [2, 2]

    # Comparison with the incremental route: reduce_circular on the same diagram
    # (without the floating-components pre-step) must give the same α via `_fr_barbell`.
    combo, _ = reduce_circular(circular_decorated(g))
    @test length(collect(pairs_of(combo))) == 1
    (dd, c) = first(pairs_of(combo))
    @test c * dd.outer_label *
          (isempty(dd.region_labels) ? one(SoergelPoly) :
           reduce(*, dd.region_labels; init = one(SoergelPoly))) == alpha(1)
end

@testset "no floating component: no-op" begin
    # a single dot of colour 1 on leaf 1 — it hangs on the boundary, hence
    # NOT floating freely.
    n = circular_node([1])
    g = CircularGraph(CircularWord([1]), [n], [Edge(1, Leaf(1), NodePort(1, 1))])
    fdm = CircularDecoratedMorphism(CircularMorphismGraph(g, 0, 1),
                                fill(one(SoergelPoly), region_count(g)))
    fnew, factor = circular_extract_floating_components(fdm)
    @test factor == one(SoergelPoly)
    @test fnew === fdm                             # returned unchanged
end

@testset "regression: pairing(dl_f, dl_e) on 121212 — the same sum, now clean" begin
    x = [1, 2, 1, 2, 1, 2]
    t = dl(x)
    e = [1, 1, 1, 1, 1, 1]
    f = [0, 0, 0, 0, 0, 0]
    dl_e = dl_entries(t)[findfirst(en -> en.e == e && en.f == e, dl_entries(t))]
    dl_f = dl_entries(t)[findfirst(en -> en.e == f && en.f == f, dl_entries(t))]

    h = compose(dl_f.morphism, dl_e.morphism)
    fm = circular_morphism(h)
    fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))

    # the new step applies ALREADY on the raw, unreduced diagram: the 6-armed braid with
    # 6 dots (floating free, no leaf, no other node) is isolated and evaluated in one
    # go.
    fnew, factor = circular_extract_floating_components(fdm)
    @test factor * fnew.outer_label == alpha(1)^2 * alpha(2) + alpha(1) * alpha(2)^2
    @test length(fnew.m.graph.nodes) == 7          # only the 6-armed braid half remains

    p = pairing(dl_f, dl_e)
    @test length(p.coeffs) == 1
    (_, coeff) = first(p.coeffs)
    @test coeff == alpha(1)^2 * alpha(2) + alpha(1) * alpha(2)^2

    for (_, fd) in p.leaves
        @test euler(fd.graph) == 2
        @test isempty(check_wiring(fd.graph))
    end
end

# ---- labels outside and on regions without a boundary gap ------
#
# On the CIRCULAR side there is NO label slot for the outer region —
# `region_count == cs.ncells` does not count the outer face (`cell_of_face[outer] == 0`).
# A label on a gapless region becomes a scalar — the same treatment the component's
# own value gets.
@testset "circular: a label on a gapless region becomes a scalar instead of a throw" begin
    NP = DiagrammaticHecke.NodePort
    b1 = circular_node([1]); b2 = circular_node([1])
    g = CircularGraph(CircularWord([2, 2]), [b1, b2],
        [Edge(2, Leaf(1), Leaf(2)),
         Edge(1, NP(1, 1), NP(2, 1))])
    @test region_count(g) == 3
    @test isempty(regions(g)[3].gaps)          # the phantom region around the barbell

    labels = fill(one(SoergelPoly), region_count(g))
    labels[3] = alpha(1)
    fdm = CircularDecoratedMorphism(CircularMorphismGraph(g, 0, 2), labels)

    fnew, factor = circular_extract_floating_components(fdm)
    # α₁ from the barbell itself TIMES α₁ from the label of the dying region
    @test factor * fnew.outer_label == alpha(1)^2
    @test length(fnew.m.graph.nodes) == 0
    @test all(isone, fnew.region_labels)

    # and the same value through the whole pipeline
    combo = reduce_to_circular_leave(CircularDecoratedMorphism(CircularMorphismGraph(g, 0, 2), labels))
    @test length(collect(pairs_of(combo))) == 1
    (dd2, coeff) = first(pairs_of(combo))
    @test coeff * dd2.outer_label == alpha(1)^2
end

@testset "circular_canonical_key separates label variants" begin
    gm = circular_morphism(light_leaf([1, 2, 1], [0, 1, 1])).graph
    lab0 = fill(one(SoergelPoly), region_count(gm))
    lab1 = copy(lab0); lab1[end] = alpha(1)
    d0, d1 = circular_decorated(gm, lab0), circular_decorated(gm, lab1)
    @test circular_canonical_key(d0) != circular_canonical_key(d1)
    @test d0 != d1
    # two terms differing ONLY by a scalar therefore do not collapse — the canonical
    # key distinguishes them by design.
end

