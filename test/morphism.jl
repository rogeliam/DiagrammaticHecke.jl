# morphism layer: braid moves, generators, light/double leaves, Zamolodchikov,
# and the diagram operations degree/hflip/shift/flip/associativity.

@testset "braid moves & reex graph (A₃)" begin
    # same element / canonical short-lex (the A₃ normal form of algebra/Coxeter.jl)
    @test same_element([1,2,3,2,1], [1,3,2,3,1])
    @test canonical_word([1,3,2,3,1]) == [1,2,3,2,1]
    @test !same_element([1,2], [2,1])
    @test is_reduced([1,2,1,3,2,1])
    @test !is_reduced([1,1])

    # braid_to_word: 12321 -> 13231 via 232->323 at (2,4)
    mv = braid_to_word([1,2,3,2,1], [1,3,2,3,1])
    @test mv == [(2,4)]
    @test replay_braid_moves([1,2,3,2,1], mv) == [1,3,2,3,1]
    @test_throws ErrorException braid_to_word([1,2,3,2,1], [1,2,3])   # diff element

    # braid_to_end_with: reorder to end in a chosen right descent
    w, m = braid_to_end_with([1,2,3,2,1], 3)
    @test w[end] == 3
    @test replay_braid_moves([1,2,3,2,1], m) == w
    @test same_element(w, [1,2,3,2,1]) && is_reduced(w)
    @test braid_to_end_with([1,2,3,2,1], 1) == ([1,2,3,2,1], Tuple{Int,Int}[])  # already ends in 1
    @test_throws ErrorException braid_to_end_with([1,2,3,2,1], 2)   # 2 not a right descent

    # apply_braid_move directly
    @test apply_braid_move([1,2,1], 1, 3) == [2,1,2]
    @test apply_braid_move([1,3], 1, 2) == [3,1]

    # braid_move_sites: 121321 has the 121->212 block (1,3) and the 13 commutation (3,4)
    @test braid_move_sites([1,2,1,3,2,1]) == [(1,3), (3,4)]
    @test braid_move_sites([1,2,1]) == [(1,3)]

    # reex graph: 121->{121,212}; w₀=121321 has 16 reduced words
    @test Set(reduced_words([1,2,1])) == Set([[1,2,1],[2,1,2]])
    @test Set(reduced_words([1,3]))  == Set([[1,3],[3,1]])
    words, edges = reex_graph([1,2,1,3,2,1])
    @test length(words) == 16
    @test all(same_element(w, [1,2,1,3,2,1]) && is_reduced(w) for w in words)
    # every edge is a genuine single braid move between listed words
    @test all(apply_braid_move(words[a], i, j) == words[b] for (a,b,(i,j)) in edges)
end

@testset "braid move as a morphism (join a vertex onto a map)" begin
    f = light_leaf_up([1,2,3,2,1])            # a map with top = 12321
    # one braid move = one new braid node; bottom unchanged, top = v'
    g = braid_top(f, 2, 4)                     # 232 -> 323
    @test bottom(g) == [1,2,3,2,1]
    @test top(g) == [1,3,2,3,1]
    @test length(g.graph.nodes) == 1
    @test g.graph.nodes[1].kind === :braid
    @test is_planar_braid(g.graph)

    # braid_move_morphism standalone: bottom = v, top = v'
    bm = braid_move_morphism([1,2,3,2,1], 2, 4)
    @test bottom(bm) == [1,2,3,2,1] && top(bm) == [1,3,2,3,1]

    # a whole sequence stacks one node per move
    moves = braid_to_word([1,2,1,3,2,1], [3,2,1,3,2,3])
    h = braid_top_to(light_leaf_up([1,2,1,3,2,1]), [3,2,1,3,2,3])
    @test top(h) == [3,2,1,3,2,3]
    @test length(h.graph.nodes) == length(moves)
    @test all(nd.kind === :braid for nd in h.graph.nodes)
end

@testset "tensor of morphisms + identity_morphism" begin
    t1 = tensor(identity_strand(1), identity_strand(2))
    @test bottom(t1) == [1, 2] && top(t1) == [1, 2]

    t2 = tensor(identity_morphism(), identity_strand(1))
    @test bottom(t2) == bottom(identity_strand(1)) && top(t2) == top(identity_strand(1))
    t2b = tensor(identity_strand(1), identity_morphism())
    @test bottom(t2b) == bottom(identity_strand(1)) && top(t2b) == top(identity_strand(1))

    t3 = tensor(dot_morphism(1), identity_strand(2))
    @test bottom(t3) == [1, 2] && top(t3) == [2]

    # chaining ≥ 2 tensors must NOT reverse the piece order (a naive "reverse both
    # tops" rule only happens to work when tops have length ≤ 1)
    cur = identity_morphism()
    for i in [1, 2, 1]
        cur = tensor(cur, identity_strand(i))
    end
    @test bottom(cur) == [1, 2, 1] && top(cur) == [1, 2, 1]
end

@testset "merge/cap generators" begin
    @test bottom(merge_morphism(1)) == [1, 1] && top(merge_morphism(1)) == [1]
    @test bottom(cap_morphism(1)) == [1, 1] && top(cap_morphism(1)) == Int[]

    c = compose(tensor(identity_strand(1), identity_strand(1)), merge_morphism(1))
    @test bottom(c) == [1, 1] && top(c) == [1]

    # split_morphism: dual of merge_morphism, bottom = [i], top = [i,i]
    @test bottom(split_morphism(1)) == [1] && top(split_morphism(1)) == [1, 1]
    @test bottom(split_morphism(2)) == [2] && top(split_morphism(2)) == [2, 2]

    sc = compose(split_morphism(1), merge_morphism(1))
    @test bottom(sc) == [1] && top(sc) == [1]

    cs = compose(merge_morphism(1), split_morphism(1))
    @test bottom(cs) == [1, 1] && top(cs) == [1, 1]
end

@testset "subexpressions / decorations / defect / expressed_word" begin
    @test length(subexpressions(3)) == 8
    @test Set(subexpressions(2)) == Set([[0,0],[1,0],[0,1],[1,1]])

    # hand-checked (D0/D1 convention): word=121
    @test decorations([1,2,1], [1,0,0]) == [:U1, :U0, :D0]
    @test expressed_word([1,2,1], [1,0,0]) == [1]
    @test defect([1,2,1], [1,0,0]) == 0

    @test decorations([1,2,1], [1,0,1]) == [:U1, :U0, :D1]
    @test expressed_word([1,2,1], [1,0,1]) == Int[]
    @test defect([1,2,1], [1,0,1]) == 1
end

@testset "light_leaf (full U0/U1/D0/D1 induction)" begin
    # trivial 1-letter cases coincide with the plain generators
    l1 = light_leaf([1], [1])
    @test bottom(l1) == bottom(identity_strand(1)) && top(l1) == top(identity_strand(1))
    l0 = light_leaf([1], [0])
    @test bottom(l0) == bottom(dot_morphism(1)) && top(l0) == top(dot_morphism(1))

    # word=121, e=100 -> U1,U0,D0, top=[1] (no braiding needed: already ends in 1)
    l = light_leaf([1,2,1], [1,0,0])
    @test bottom(l) == [1,2,1] && top(l) == [1]

    # word=121, e=101 -> U1,U0,D1, top=ε (needs a real cap, hits the join_graphs
    # chain-fuse case)
    l2 = light_leaf([1,2,1], [1,0,1])
    @test bottom(l2) == [1,2,1] && top(l2) == Int[]

    # a D0 case that needs NO braiding: word=2131, e=1100 — at k=4 (letter 1) the
    # current top is [2,1] and already ends in 1, so the merge acts directly.
    l3 = light_leaf([2,1,3,1], [1,1,0,0])
    @test bottom(l3) == [2,1,3,1]
    @test top(l3) == expressed_word([2,1,3,1], [1,1,0,0])   # == [2,1]

    # a D1 case that genuinely BRAIDS first: word=1212, e=1111 — at k=4 (letter 2)
    # the current top is [1,2,1] (short-lex, ends in 1), so braid_to_end_with must
    # reorder it to [2,1,2] (ends in 2) before the cap absorbs the last two strands.
    l4 = light_leaf([1,2,1,2], [1,1,1,1])
    @test bottom(l4) == [1,2,1,2]
    @test top(l4) == expressed_word([1,2,1,2], [1,1,1,1])   # == [2,1]
    @test any(nd.kind === :braid for nd in l4.graph.nodes)
    @test is_planar_braid(l4.graph)

    # every subexpression of a length-4 word gives a consistent light leaf
    # (the @assert inside light_leaf is the real check; this just exercises all of them)
    for e in subexpressions(4)
        m = light_leaf([2,1,3,1], e)
        @test bottom(m) == [2,1,3,1]
        @test top(m) == expressed_word([2,1,3,1], e)
    end
end

@testset "double_leaf / light_leaves / double_leaves" begin
    # a double leaf of a leaf with itself is a genuine endomorphism x -> x
    dl = double_leaf([1,2,1], [1,1,1], [1,2,1], [1,1,1])
    @test bottom(dl) == [1,2,1] && top(dl) == [1,2,1]

    # domain mismatch (different expressed elements) errors
    @test_throws ErrorException double_leaf([1,2,1], [1,1,1], [1,2,1], [1,0,0])

    # light_leaves: all e with expressed_word == z, degree filter counts correctly
    lls = light_leaves([1,2,1], [1])
    @test Set(e for (e,_) in lls) == Set([[1,0,0], [0,0,1]])
    @test length(light_leaves([1,2,1], [1]; degree = 0)) == 1
    @test length(light_leaves([1,2,1], [1]; degree = 2)) == 1

    # double_leaves x->x: total count and a fixed-degree slice, cross-checked by
    # hand against the 8 subexpressions of 121 (see the decorations test above)
    dls = double_leaves([1,2,1], [1,2,1])
    @test length(dls) == 12
    @test all(d.degree == defect([1,2,1], d.e) + defect([1,2,1], d.f) for d in dls)
    @test length(double_leaves([1,2,1], [1,2,1]; degree = 0)) == 2
end

@testset "flip on a degenerate (empty top/bottom) MorphismGraph" begin
    # MorphismGraph's cut1==cut2 shortcut can only ever mean "bottom = whole boundary,
    # top = ε" — flipping such a morphism must produce bottom = ε, top = whole
    # boundary (not a silent no-op). Needed by every double_leaf through the identity
    # element z = ε.
    fd = flip(dot_morphism(1))
    @test bottom(fd) == Int[] && top(fd) == [1]
    @test bottom(flip(fd)) == bottom(dot_morphism(1)) && top(flip(fd)) == top(dot_morphism(1))

    l = light_leaf([1,2,1], [0,0,0])         # top = ε
    fl = flip(l)
    @test bottom(fl) == Int[] && top(fl) == [1,2,1]
end

@testset "Zamolodchikov cycle (A₃)" begin
    ws = zamo_words()
    @test length(ws) == 14                        # 14 distinct words, cyclic
    @test ws[1] == [1,2,1,3,2,1]                  # start = canonical short-lex w₀
    @test ws[8] == [3,2,1,3,2,3]                  # opposite corner of the hexagon
    @test all(same_element(w, [1,2,1,3,2,1]) && is_reduced(w) for w in ws)  # all w₀
    # every consecutive pair is a single braid move, INCLUDING the closing 14 → 1
    for k in 1:14
        a, b = ws[k], ws[mod1(k+1, 14)]
        @test any(apply_braid_move(a, i, j) == b for (i, j) in braid_move_sites(a))
    end
    # Zamo(i,j): forward walk, bottom = ws[i], top = ws[j]; i>j wraps; Zamo(i,i)=cycle
    f = Zamo(1, 8)
    @test bottom(f) == [1,2,1,3,2,1] && top(f) == [3,2,1,3,2,3]   # LHS half, upward
    @test length(f.graph.nodes) == 7
    @test length(Zamo(1, 1).graph.nodes) == 14    # Zamo(i,i) = whole cycle, NOT id
    @test bottom(Zamo(1,1)) == top(Zamo(1,1)) == [1,2,1,3,2,1]
    @test length(Zamo(1, 2).graph.nodes) == 1     # one step
    @test length(Zamo(10, 6).graph.nodes) == 10   # i>j wraps: mod(6-10,14)=10
    @test bottom(Zamo(10,6)) == ws[10] && top(Zamo(10,6)) == ws[6]
    # LHS = RHS relation: both 121321 → 321323, endpoints agree (necessary check)
    @test bottom(zamo_lhs()) == bottom(zamo_rhs()) == [1,2,1,3,2,1]
    @test top(zamo_lhs()) == top(zamo_rhs()) == [3,2,1,3,2,3]
    @test zamo_relation_endpoints_agree()
    @test length(zamo_lhs().graph.nodes) == 7 && length(zamo_rhs().graph.nodes) == 7
end

@testset "degree on diagrams (Soergel grading, dots − trivalents)" begin
    # generators: a plain dot / trivalent
    @test dot_count(dot(1)) == 1 && trivalent_count(dot(1)) == 0
    @test degree(dot(1)) == 1
    @test dot_count(trivalent(1)) == 0 && trivalent_count(trivalent(1)) == 1
    @test degree(trivalent(1)) == -1
    # :braid nodes contribute 0
    bm = braid_move_morphism([1,2,3,2,1], 2, 4)
    @test dot_count(bm.graph) == 0 && trivalent_count(bm.graph) == 0
    @test degree(bm) == 0

    # degree(a::SoergelPoly, m) = 2*degree(a) + degree(m)
    p = alpha(1) * alpha(2)          # degree 4
    @test degree(p, dot(1)) == 2 * 4 + 1

    # the sharp check: diagrammatic degree(d.morphism) == combinatorial d.degree
    # for every double leaf of every A3 word pair of length <= 3
    words = [w for w in a3_words() if length(w) <= 3]
    checked = 0
    for x in words, y in words
        for d in double_leaves(x, y)
            checked += 1
            @test degree(d.morphism) == d.degree
        end
    end
    @test checked > 0
end

@testset "hflip (horizontal mirror, left↔right)" begin
    words = [w for w in a3_words() if length(w) <= 3]

    # canonical-key + cuts + boundaries equality helper for MorphismGraph
    morph_eq(a::MorphismGraph, b::MorphismGraph) =
        canonical_key(a.graph) == canonical_key(b.graph) &&
        a.cut1 == b.cut1 && a.cut2 == b.cut2 &&
        bottom(a) == bottom(b) && top(a) == top(b)

    checked = 0
    commute_checked = 0
    for x in words, y in words
        for d in double_leaves(x, y)
            m = d.morphism
            hm = hflip(m)
            checked += 1
            # the contract: bottom stays bottom, top stays top,
            # only each word's letter order reverses.
            @test bottom(hm) == reverse(bottom(m))
            @test top(hm) == reverse(top(m))
            @test morph_eq(hflip(hm), m)                    # involutive
            @test is_planar_braid(hm.graph)                 # catches a wrong slot offset
            @test degree(hm) == degree(m)                   # dots/trivalents preserved
            # flip(hflip(m)) == hflip(flip(m)): the two mirrors commute. Checked
            # WITHOUT excluding the cut1==cut2 degenerate boundary (bottom=ε or
            # top=ε, n>1) — `flip` must handle this case correctly (see its
            # docstring).
            commute_checked += 1
            @test morph_eq(flip(hm), hflip(flip(m)))
        end
    end
    @test checked > 0
    @test commute_checked > 0

    # the degenerate flip(dot_morphism(1)) case
    fd = flip(dot_morphism(1))
    hfd = hflip(fd)
    @test bottom(hfd) == reverse(bottom(fd)) && top(hfd) == reverse(top(fd))
    @test morph_eq(hflip(hfd), fd)
    @test is_planar_braid(hfd.graph)
    @test degree(hfd) == degree(fd)
    # commutation must hold here too (this IS a cut1==cut2 degenerate case)
    @test morph_eq(flip(hfd), hflip(flip(fd)))
end

@testset "shift_right / shift_left (bottom↔top letter reinterpretation)" begin
    words = [w for w in a3_words() if length(w) <= 3]

    # exact equality helper: shift_right/shift_left must NEVER touch the graph
    # (only cut2 moves), so require the underlying WordGraph to be IDENTICAL
    # (not just canonical_key-equal), plus matching cuts/boundaries.
    same_morph(a::MorphismGraph, b::MorphismGraph) =
        a.graph == b.graph && a.cut1 == b.cut1 && a.cut2 == b.cut2

    checked = 0
    for x in words, y in words
        for d in double_leaves(x, y)
            m = d.morphism
            blen = length(bottom(m)); tlen = length(top(m))
            blen == 0 && tlen == 0 && continue    # nothing to shift either way
            checked += 1

            # graph untouched by shift in either direction
            if blen > 0
                r = shift_right(m)
                @test canonical_key(r.graph) == canonical_key(m.graph)
                @test bottom(r) == bottom(m)[1:end-1]
                @test top(r) == vcat(top(m), bottom(m)[end])
                @test is_planar_braid(r.graph)
                @test degree(r) == degree(m)
                @test same_morph(shift_left(r), m)     # inverse
            end
            if tlen > 0
                l = shift_left(m)
                @test canonical_key(l.graph) == canonical_key(m.graph)
                @test top(l) == top(m)[1:end-1]
                @test bottom(l) == vcat(bottom(m), top(m)[end])
                @test is_planar_braid(l.graph)
                @test degree(l) == degree(m)
                @test same_morph(shift_right(l), m)    # inverse
            end
            # k-fold == repeated single steps
            if blen >= 2
                @test same_morph(shift_right(m, 2), shift_right(shift_right(m)))
            end
            if tlen >= 2
                @test same_morph(shift_left(m, 2), shift_left(shift_left(m)))
            end
            # shifting more than available errors (conservative choice, no
            # cyclic wrap-around)
            @test_throws ArgumentError shift_right(m, blen + 1)
            @test_throws ArgumentError shift_left(m, tlen + 1)
        end
    end
    @test checked > 0

    # degenerate single cases named in the task
    d = dot_morphism(1)                      # bottom=[1], top=ε
    @test_throws ArgumentError shift_left(d)  # top empty, nothing to shift
    rd = shift_right(d)                       # bottom empties -> _EMPTY_BOTTOM_CUT sentinel
    @test bottom(rd) == Int[] && top(rd) == [1]
    @test rd.cut1 == DiagrammaticHecke._EMPTY_BOTTOM_CUT && rd.cut2 == DiagrammaticHecke._EMPTY_BOTTOM_CUT
    @test same_morph(shift_left(rd), d)        # inverse round-trips out of the sentinel

    fd = flip(dot_morphism(1))                # bottom=ε, top=[1] (sentinel input)
    @test_throws ArgumentError shift_right(fd) # bottom empty, nothing to shift
    lfd = shift_left(fd)
    @test bottom(lfd) == [1] && top(lfd) == Int[]
    @test same_morph(shift_right(lfd), fd)

    im = identity_morphism()                  # ε → ε
    @test_throws ArgumentError shift_right(im)
    @test_throws ArgumentError shift_left(im)

    cm = cap_morphism(1)                      # bottom=[1,1], top=ε
    @test_throws ArgumentError shift_left(cm)
    rc = shift_right(cm)
    @test bottom(rc) == [1] && top(rc) == [1]
    @test same_morph(shift_left(rc), cm)

    # interaction with flip/hflip (a property of the implementation, not enforced
    # by the code): flip(shift_right(m))
    # and shift_left(flip(m)) agree on bottom/top WORDS (not on the exact graph —
    # flip remaps the graph, shift_right never does), because flip's word-order
    # preservation composes with "last of bottom -> last of top" the same way on
    # both sides. hflip does NOT commute with shift_right/shift_left in general:
    # hflip reverses letter order, so "last letter" under hflip corresponds to the
    # ORIGINAL FIRST letter, not the last — shift_right/shift_left always act on
    # the last letter, so there is no simple commutation with hflip except for
    # palindromic words. Both are checked below on double_leaf([1,2,1],[1,3,2]).
    m = double_leaf([1,2,1],[0,0,0],[1,3,2],[0,0,0])
    a = flip(shift_right(m)); b = shift_left(flip(m))
    @test bottom(a) == bottom(b) && top(a) == top(b)     # flip commutes (words agree)
    c = hflip(shift_right(m)); e = shift_right(hflip(m))
    @test !(bottom(c) == bottom(e) && top(c) == top(e))  # hflip does NOT commute
end

@testset "associativity (A1)" begin
    # A1: three
    # constructions of the SAME morphism 11 -> 11 (two trivalents, no dots) must
    # agree. NOT a reduction rule (both sides already irreducible) — this test
    # only documents the identity, nothing is wired into GRAPH_RULES.

    # A: the light-leaf route
    A = double_leaf([1,1],[1,0],[1,1],[1,0])
    @test bottom(A) == [1,1] && top(A) == [1,1]
    @test dot_count(A.graph) == 0 && trivalent_count(A.graph) == 2

    # B: "first-with-first, last-with-last, then connect" directly from generators
    m = merge_morphism(1)
    B = compose(m, flip(m))
    @test bottom(B) == [1,1] && top(B) == [1,1]
    @test dot_count(B.graph) == 0 && trivalent_count(B.graph) == 2

    # C: same graph as B, both cuts rotated by +1 (the "shift bottom/top by one"
    # reading; NOT shift_right/shift_left, which move only cut2 and change the
    # bottom/top word lengths)
    n = length(B.graph.word)
    C = MorphismGraph(B.graph, mod(B.cut1 + 1, n), mod(B.cut2 + 1, n))
    @test bottom(C) == [1,1] && top(C) == [1,1]
    @test dot_count(C.graph) == 0 && trivalent_count(C.graph) == 2

    # canonical_key equality (documented comparator, diagram/Combo.jl)
    @test canonical_key(A.graph) == canonical_key(B.graph)
    @test canonical_key(B.graph) == canonical_key(C.graph)
    @test canonical_key(A.graph) == canonical_key(C.graph)

    # stronger manual isomorphism evidence: A and B even share IDENTICAL raw
    # nodes/edges (no relabelling needed here, so the known canonical_key gap —
    # array-index node fingerprints — does not come into play).
    @test A.graph.nodes == B.graph.nodes
    @test A.graph.edges == B.graph.edges
end
