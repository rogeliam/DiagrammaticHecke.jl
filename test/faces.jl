# Focused checks for the cell recognition (Faces.jl): the cells attribute is
# filled by the plain 3-argument constructor, and the counts match the
# outside-in picture (dot merges neighbour gaps; node legs cut off sectors).

@testset "cells / faces (step 2)" begin
    # attribute is there and consistent with the API
    g = dot(1)
    @test g.cells.ncells == face_count(g) == 1
    @test g.cells.gap_cell == [1]

    # barbell: internal edge, empty boundary → one cell
    @test face_count(barbell_graph()) == 1

    # trivalent: three legs to the boundary → three cells, all gaps distinct
    t = trivalent(1)
    @test face_count(t) == 3
    @test length(unique(t.cells.gap_cell)) == 3
    @test all(!iszero, t.cells.sector_cell[1])

    # a single strand (chord) splits the disc into two cells
    chord = WordGraph(CircularWord([1, 1]), Node[], Edge[Edge(1, Leaf(1), Leaf(2))])
    @test face_count(chord) == 2

    # two parallel strands → three cells; the middle one spans gaps 2 and 4
    str = WordGraph(CircularWord([1, 1, 1, 1]), Node[],
                    Edge[Edge(1, Leaf(1), Leaf(2)), Edge(1, Leaf(3), Leaf(4))])
    @test face_count(str) == 3
    @test str.cells.gap_cell[2] == str.cells.gap_cell[4]
    @test length(unique(str.cells.gap_cell)) == 3

    # a dot MERGES its two neighbour gaps: two dots side by side → one cell
    dots = WordGraph(CircularWord([1, 1]),
                     Node[Node(:dot, [1], 0), Node(:dot, [1], 0)],
                     Edge[Edge(1, Leaf(1), NodePort(1, 1)),
                          Edge(1, Leaf(2), NodePort(2, 1))])
    @test face_count(dots) == 1
    @test dots.cells.gap_cell == [1, 1]

    # braid nodes: one cell per boundary gap (tree with 4 / 6 arms)
    @test face_count(braid(1, 3; m = 2)) == 4
    @test face_count(braid(1, 2; m = 3)) == 6

    # crossed needle: the lens between the two internal edges is a cell of its own
    # (it touches NO boundary gap; needle_a_graph's parallel wiring is twisted)
    nx = WordGraph(CircularWord([1, 1]),
                   Node[Node(:trivalent, [1], 0), Node(:trivalent, [1], 0)],
                   Edge[Edge(1, NodePort(1, 1), NodePort(2, 2)),
                        Edge(1, NodePort(1, 2), NodePort(2, 1)),
                        Edge(1, Leaf(1), NodePort(1, 3)),
                        Edge(1, Leaf(2), NodePort(2, 3))])
    @test face_count(nx) == 3
    @test count(c -> isempty(c.gaps), inner_faces(nx)) == 1     # the lens

    # a free circle adds one cell (its interior)
    loop = WordGraph(CircularWord([1]),
                     Node[Node(:dot, [1], 0)],
                     Edge[Edge(1, Leaf(1), NodePort(1, 1)), Edge(2, Circle(2), Circle(2))])
    @test face_count(loop) == 2

    # face_of_port: the two sides of the edge at a port.
    #
    # ⚡ There is NO absolute "left" and "right" — by "left" one only means the
    # marker; otherwise there is only further away or closer in.
    # `:left`/`:right` merely name the two sides of a dart; WHICH one is called which
    # depends on the orientation and is meaningless without a marker.
    #
    # What is checked is what lasts: the edge at leaf k separates the cells of gaps
    # k-1 and k, and the two sides are different.
    @test face_of_port(str, Leaf(1)) != face_of_port(str, Leaf(1); side = :right)
    @test Set([face_of_port(str, Leaf(1)), face_of_port(str, Leaf(1); side = :right)]) ==
          Set([gap_cell(str, 1), gap_cell(str, 2)])
    # And the substantive claim, WITH a marker: one of the two sides lies closer to
    # the marker than the other.
    let fm = circular_morphism(morphism_graph(str, 0, 2))
        dist = circular_region_distances(fm)
        @test dist[face_of_port(str, Leaf(1))] != dist[face_of_port(str, Leaf(1); side = :right)]
    end
    @test face_of_port(chord, Leaf(1)) != face_of_port(chord, Leaf(1); side = :right)
    @test_throws ArgumentError face_of_port(str, Leaf(1); side = :up)

    # Euler sanity: V − E + F = 2 (sphere; boundary arcs count as edges)
    for gg in (dot(1), trivalent(1), braid(1, 2; m = 3), str, nx, dots)
        V = length(gg.word) + length(gg.nodes)
        E = length(gg.edges) + length(gg.word)
        @test V - E + (face_count(gg) + 1) == 2
    end

    # MorphismGraph forwards to its graph (cuts do not change the cells)
    @test face_count(morphism_graph(braid(1, 3; m = 2), 0, 2)) == 4
end

# A second, parallel notion of a cell: the boundary is a WALL. `face_count` counts
# faces of the embedding (a face may touch the boundary several times, the parts being
# connected on the outside), `region_count` counts what one sees in the picture as
# separate areas. See the comment block in diagram/Faces.jl.
@testset "regions (solid boundary)" begin
    # Where every leaf is a wall, the two notions agree.
    for gg in (dot(1), trivalent(1), braid(1, 3; m = 2), braid(1, 2; m = 3))
        @test region_count(gg) == face_count(gg)
    end

    # Two parallel strands (chords 1—2 and 3—4). Gaps 2 and 4 BOTH lie in the band
    # between the chords and are connected there in the interior — the inner path runs
    # only over leaf–leaf chords, across no node. So ONE region, NOT 4 — that would be
    # the error of the pure boundary walk.
    str = WordGraph(CircularWord([1, 1, 1, 1]), Node[],
                    Edge[Edge(1, Leaf(1), Leaf(2)), Edge(1, Leaf(3), Leaf(4))])
    @test face_count(str) == 3
    @test region_count(str) == 3
    @test boundary_regions(str) == [[1], [2, 4], [3]]

    # A dot is NOT a wall: it merges its two neighbouring gaps.
    dots = WordGraph(CircularWord([1, 1]),
                     Node[Node(:dot, [1], 0), Node(:dot, [1], 0)],
                     Edge[Edge(1, Leaf(1), NodePort(1, 1)),
                          Edge(1, Leaf(2), NodePort(2, 1))])
    @test boundary_regions(dots) == [[1, 2]]
    @test region_count(dots) == 1

    # A leaf–leaf chord IS a wall (no dot at its end).
    chord = WordGraph(CircularWord([1, 1]), Node[], Edge[Edge(1, Leaf(1), Leaf(2))])
    @test boundary_regions(chord) == [[1], [2]]

    # Inner cells without a boundary gap still count (the lens in needle_a).
    nx = WordGraph(CircularWord([1, 1]),
                   Node[Node(:trivalent, [1], 0), Node(:trivalent, [1], 0)],
                   Edge[Edge(1, NodePort(1, 1), NodePort(2, 2)),
                        Edge(1, NodePort(1, 2), NodePort(2, 1)),
                        Edge(1, Leaf(1), NodePort(1, 3)),
                        Edge(1, Leaf(2), NodePort(2, 3))])
    @test region_count(nx) == 3          # 2 boundary regions + 1 lens
    @test length(boundary_regions(nx)) == 2

    # Regions REFINE faces, never the other way round.
    for gg in (dot(1), trivalent(1), braid(1, 2; m = 3), str, nx, dots, chord)
        @test region_count(gg) >= face_count(gg)
    end

    # ⚡ The NODE RULE, pinned to the three cases with which it was fixed. The
    # boundary word is `1 1 2 1 1 2` throughout, so 6 gaps; only the wiring inside
    # differs.
    #
    # These three are why the pure boundary walk does not suffice: the second and the
    # third case BOTH have a face with the gaps {3,6}, and in both the face cycle leads
    # through the interior from one to the other. What separates them is only the
    # question of WHAT the path runs through.
    dl121(e) = circular_double_leaf([1, 2, 1], e, [1, 2, 1], e)

    # (a) six dots: no leaf is a wall ⇒ everything is ONE area.
    @test region_count(dl121([0, 0, 0])) == 1
    @test boundary_regions(dl121([0, 0, 0])) == [[1, 2, 3, 4, 5, 6]]

    # (b) two leaf–leaf chords (Leaf1—Leaf3, Leaf6—Leaf4): face 2 has the gaps {3,6}
    # and NO sectors — it is the band between the chords, and gaps 3 and 6 are its two
    # ends. The path between them runs across no node ⇒ they merge ⇒ 3 regions (not
    # 4).
    b = dl121([1, 0, 1])
    @test face_count(b.graph) == 3
    @test region_count(b) == 3
    @test boundary_regions(b) == [[1, 2], [3, 6], [4, 5]]

    # (c) the same gaps {3,6} in face 2, but the path runs through the trivalents
    # n2/n4 (degree 3) ⇒ these separate ⇒ 4 regions. The counterpart to (b):
    # combinatorially the same cycle type, a different result.
    #
    # The boundary separates in the face walk itself (orientation block in
    # diagram/Faces.jl), which is why `face_count == 4` stands here and NOT 2 — the 2
    # would arise because two areas separate in the picture were connected around
    # the outside of the disc.
    c = dl121([1, 0, 0])
    @test face_count(c.graph) == 4
    @test region_count(c) == 4
    @test boundary_regions(c) == [[1, 2], [3], [4, 5], [6]]

    # `regions(g)`: numbered through, with the face id that polynomial labels and BFS
    # distances still hang on. Several regions may share the same `cell` — precisely
    # the case where ONE number may stand for several visible areas.
    #
    # The ORDER of the regions is meaningless — it follows the cell numbering. So the
    # SET of gap groups is compared; `boundary_regions` still gives the version sorted
    # by smallest gap.
    rs = regions(c)
    @test length(rs) == region_count(c) == 4
    @test [R.id for R in rs] == 1:4
    @test Set(Set.(R.gaps for R in rs)) == Set(Set.([[1, 2], [3], [4, 5], [6]]))
    @test boundary_regions(c) == [[1, 2], [3], [4, 5], [6]]
    # A region is exactly one cell, so `R.cell == R.id`. ([1, 2, 1, 2] would mean two
    # regions sharing a face each — exactly the defect this guards against.)
    @test [R.cell for R in rs] == [R.id for R in rs]
    @test all(R -> 1 <= R.cell <= face_count(c.graph), rs)

    # Inner regions without a boundary gap are ordinary cells (the lens in needle_a):
    # two boundary sectors and the lens between them. The order does not matter.
    @test Set(Set.(R.gaps for R in regions(nx))) == Set(Set.([[1], [2], Int[]]))
    @test count(R -> isempty(R.gaps), regions(nx)) == 1

    # The triggering case: the double leaf on 12321. Without the right orientation, 5
    # faces would stand against 9 regions — that would be the defect, not the
    # mathematics. With the right orientation the two coincide.
    fm = circular_double_leaf([1, 2, 3, 2, 1], [1, 0, 1, 0, 0],
                         [1, 2, 3, 2, 1], [1, 0, 1, 0, 0])
    @test face_count(fm.graph) == 9
    @test region_count(fm) == 9
    @test boundary_regions(fm) == [[1, 2], [3, 4], [5], [6, 7], [8, 9], [10]]
    @test length(regions(fm)) == 9

    # Regions REFINE faces: every region lies in exactly ONE face. (The test for the
    # node rule — it may only merge what lies in the same face anyway, never across
    # face boundaries.)
    for gg in (b, c, fm)
        gc = gg.graph.cells.gap_cell
        for reg in boundary_regions(gg)
            @test length(unique(gc[k] for k in reg)) == 1
        end
    end

    # The thin and the circular picture see the same topology.
    m = light_leaf([3, 2, 1], [1, 0, 1])
    @test region_count(m) == region_count(circular_morphism(m))

    # Cuts do not change the regions.
    @test region_count(morphism_graph(braid(1, 3; m = 2), 0, 2)) ==
          region_count(braid(1, 3; m = 2))
end

@testset "cell display: numbers + marker distances" begin
    # cell_number_svg puts every cell id into the SVG
    svg = cell_number_svg(trivalent(1))
    @test all(occursin(">$i</text>", svg) for i in 1:3)

    # a MorphismGraph keeps the black marker AND gets numbered cells
    m = morphism_graph(braid(1, 3; m = 2), 0, 2)
    svgm = cell_number_svg(m)
    @test occursin("fill=\"#000\"", svgm)
    @test all(occursin(">$i</text>", svgm) for i in 1:4)

    # distances on two parallel strands: marker at gap 1 → its cell has distance 0
    str = WordGraph(CircularWord([1, 1, 1, 1]), Node[],
                    Edge[Edge(1, Leaf(1), Leaf(2)), Edge(1, Leaf(3), Leaf(4))])
    ms = morphism_graph(str, 1, 3)          # bottom = 11, top = 11
    dist = cell_distances(ms)
    @test dist[gap_cell(str, 1)] == 0       # the cell touching the marker
    @test dist[gap_cell(str, 2)] == 1       # middle cell is adjacent
    @test dist[gap_cell(str, 3)] == 2       # far cell one more step
    @test maximum(dist) == 2

    # Zamo(1,10): exactly one cell touches the marker; distances are BFS levels
    z = Zamo(1, 10)
    dz = cell_distances(z)
    @test count(==(0), dz) == 1
    @test all(d >= 0 for d in dz)
    @test maximum(dz) == 6
    svgd = cell_distance_svg(z)
    @test occursin("fill=\"#000\"", svgd) && occursin(">0</text>", svgd)
end

# ⚡ Guards against the orientation error.
#
# The face walk presupposes a PLANAR embedding; the criterion for one is
# `V − E + F = 2`. For everything out of `light_leaf`/`double_leaf`, the walk must
# stay on a genuinely planar embedding (χ = 2) — not just the hand-wired generators
# checked above, which are a special case and do not exercise the general
# construction.
@testset "the embedding is planar: V − E + F = 2" begin
    χ(g, deg) = begin
        V = length(g.word) + length(g.nodes)
        E = count(e -> !(e.a isa DiagrammaticHecke.Circle || e.b isa DiagrammaticHecke.Circle),
                  g.edges) + length(g.word)
        t = DiagrammaticHecke._trace_generic(g.word, g.nodes, g.edges, deg)
        V - E + length(t.cycles)
    end

    # the hand-wired generators and fixtures
    for gg in (dot(1), dot(2), trivalent(1), trivalent(3), braid(1, 2), braid(2, 3),
               braid(1, 3; m = 2), join_trivalents(1), DiagrammaticHecke.unit_graph(),
               DiagrammaticHecke.needle_a_graph(), DiagrammaticHecke.needle_b_graph(),
               DiagrammaticHecke.barbell_graph(), DiagrammaticHecke.r3_braid_dot(),
               DiagrammaticHecke.r4_braid_dot())
        @test χ(gg, DiagrammaticHecke._node_degree) == 2
        @test χ(circular(gg), DiagrammaticHecke._circular_degree) == 2
    end

    # and the computational world: light leaves and double leaves, plain and circular
    for x in ([1, 2, 1], [1, 3, 1, 3]), bits in 0:(2^length(x) - 1)
        e = [(bits >> (i - 1)) & 1 for i in 1:length(x)]
        g = light_leaf(x, e).graph
        @test χ(g, DiagrammaticHecke._node_degree) == 2
        @test χ(circular(g), DiagrammaticHecke._circular_degree) == 2
    end

    # needle_a: left, right and the lens in the middle = 3 faces, NOT 2 — a non-planar
    # embedding would give the wrong answer there. Likewise barbell: the disc minus an
    # edge is ONE area.
    @test face_count(DiagrammaticHecke.needle_a_graph()) == 3
    @test face_count(DiagrammaticHecke.barbell_graph()) == 1
end

# GLOBAL crossing-freeness of the Zamo diagrams: none of them uses a non-planar
# braid.
#
# Over ALL 196 pairs `Zamo(i,j)`, i,j ∈ 1:14, `euler == 2` throughout and
# `check_wiring` is EMPTY throughout. This holds at the level of the EMBEDDING:
# `check_wiring` walks the boundary in leaf order and checks whether a node's slots
# appear in cyclic order (`:leaf_ring`, plus `:node_ring`/`:reversed_node` for
# node↔node); the complete global test is `euler == 2`. What this does NOT say:
# whether the Tutte renderer also DRAWS the planar embedding without crossings —
# that is a layout question, not a wiring one.
@testset "Zamo(i,j) are globally embedded without crossings" begin
    unsauber = Tuple{Int,Int,Int,Int}[]
    for i in 1:14, j in 1:14
        g = Zamo(i, j).graph
        e, vs = euler(g), check_wiring(g)
        (e == 2 && isempty(vs)) || push!(unsauber, (i, j, e, length(vs)))
    end
    @test unsauber == Tuple{Int,Int,Int,Int}[]

    # the four cases checked individually, with their node counts
    for (i, j, nv) in ((1, 1, 14), (1, 8, 7), (1, 10, 9), (1, 14, 13))
        g = Zamo(i, j).graph
        @test length(g.nodes) == nv
        @test is_planar_embedding(g)
        @test isempty(check_wiring(g))
        @test is_planar_braid(g)     # the LOCAL slot test, for comparison
    end
end

# `euler`/`is_planar_embedding` are PUBLIC (diagram/WiringCheck.jl). The same formula
# also appears as `_circular_is_planar_embedding`
# (circular/rules/CircularInverseRules.jl), an alias — which this testset checks too.
@testset "euler / is_planar_embedding are public" begin
    g = DiagrammaticHecke.unit_graph()
    @test euler(g) == 2
    @test is_planar_embedding(g)
    @test euler(circular(g)) == 2
    # the private name points at the same answer
    @test DiagrammaticHecke._circular_is_planar_embedding(circular(g)) == is_planar_embedding(circular(g))
    # The same number as the formula V = |nodes|+|word|, E = |edges|+|word|,
    # F = face_count+1 — as long as there IS a boundary (free circles cancel in
    # both versions).
    alt(h) = (length(h.nodes) + length(h.word)) - (length(h.edges) + length(h.word)) +
             (face_count(h) + 1)
    for h in (g, DiagrammaticHecke.needle_a_graph())
        @test euler(h) == alt(h)
    end
    # ⚠️ WITHOUT A BOUNDARY that formula is WRONG: with no boundary there is
    # no separate outer face, the tracer already counts the one cycle as a cell, and
    # `face_count + 1` counts it a second time. The barbell thus comes to 3 instead of
    # 2 by that formula, wrongly suggesting a non-planar diagram. `euler` agrees with
    # the χ of the testset above, which takes `length(t.cycles)` directly.
    bar = DiagrammaticHecke.barbell_graph()
    @test euler(bar) == 2 && alt(bar) == 3
    # MorphismGraph forwards (the cuts change nothing about the cells)
    m = light_leaf([1, 2, 1], [1, 0, 1])
    @test euler(m) == euler(m.graph) == 2
end

# check_wiring computes the HOUSE CONVENTION directly (boundary ccw, arms cw) instead
# of reading it off V−E+F indirectly. Four kinds: loose leaf, leaf↔node, node↔node,
# global. The derivation is in the comment block above the function in
# diagram/WiringCheck.jl.
@testset "check_wiring — the convention computed directly" begin
    N = DiagrammaticHecke.Node; E = DiagrammaticHecke.Edge
    LP = DiagrammaticHecke.Leaf; NP = DiagrammaticHecke.NodePort

    # clean: the hand-wired generators, plain AND circular
    for gg in (dot(1), trivalent(1), braid(1, 2), braid(1, 3; m = 2),
               join_trivalents(1), DiagrammaticHecke.unit_graph(),
               DiagrammaticHecke.needle_a_graph(), DiagrammaticHecke.barbell_graph())
        @test isempty(check_wiring(gg))
        @test isempty(check_wiring(circular(gg)))
    end

    # LOOSE LEAF: a boundary word of length 2 but only one leaf wired. `is_wired` does
    # not see that (it checks node slots only), `check_wiring` does.
    los = WordGraph(CircularWord([1, 1]), [N(:dot, [1], 0)],
                    E[E(1, LP(1), NP(1, 1))])
    ds = check_wiring(los)
    @test any(v -> v.kind === :leaf_degree, ds)
    @test occursin("leaf 2", first(v.detail for v in ds if v.kind === :leaf_degree))

    # LEAF↔NODE: the R8 fixture `gA`. The braid node hangs on the
    # leaves 5, 1, 4 — in slot order therefore FORWARDS along the leaf ring, whereas the
    # convention requires backwards. Here V−E+F = −2.
    gA = WordGraph(CircularWord([3, 1, 1, 3, 1]),
        [N(:trivalent, [1], 0), N(:braid, [1, 3], 2)],
        E[E(1, LP(2), NP(1, 1)), E(1, LP(3), NP(1, 2)),
          E(1, LP(5), NP(2, 1)), E(3, LP(1), NP(2, 2)),
          E(1, NP(1, 3), NP(2, 3)), E(3, LP(4), NP(2, 4))])
    vs = check_wiring(gA)
    @test euler(gA) == -2
    @test any(v -> v.kind === :leaf_ring && v.node == 2, vs)
    @test any(v -> v.kind === :euler, vs)

    # NODE↔NODE: the same claim as a MEASUREMENT — reverse a node's arm sequence and if
    # V−E+F rises, it hung mirrored. `braid(1,3;m=2)` is clean; reversing one node there
    # by hand makes check_wiring report exactly that node.
    ok = join_trivalents(1)
    @test isempty(check_wiring(ok))
    d = DiagrammaticHecke._node_degree(ok.nodes[1])
    f(p) = p isa NP && p.node == 1 ? NP(1, d + 1 - p.slot) : p
    kaputt = WordGraph(ok.word, ok.nodes,
                       E[E(e.colour, f(e.a), f(e.b)) for e in ok.edges])
    ws = check_wiring(kaputt)
    @test any(v -> v.kind === :reversed_node && v.node == 1, ws)

    # the table outputs give something readable
    @test occursin("no violations", check_wiring_table(ok))
    @test occursin("⚠", check_wiring_table(gA))
end

