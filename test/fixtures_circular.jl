# test/fixtures_circular.jl — shared PLANAR fixtures for circular.jl/circularregionrules.jl (the
# C7/C8 and C6 self-loop fixtures defined locally in those files are non-planar as EMBEDDINGS).
#
# House convention: the boundary runs COUNTER-clockwise (leaves in word order), the arms
# at a node run clockwise (slot order). A node at several consecutive leaves ⇒ the leaves
# run BACKWARDS along the slot ring; two nodes across two edges (a bigon) ⇒ the edges
# attach counter-oriented. The test: `euler(g) == 2` and
# `isempty(check_wiring(g))`.
#
# Every fixture below was verified by brute force over all colour-faithful slot
# assignments (the same node kinds/colours/topology as the non-planar fixture —
# only the slot assignment differs).

# C7 fixture (test/circular.jl "C7 _fr_braid_back (R5)"): two back-to-back m=3 braids,
# Word [1,2,1,1,2,1]. Internal edges a4–b4/a5–b3/a6–b2, b-legs to leaves 4–6.
# Identical to the local `_planar_braid_back_fixture` from
# test/circularregionrules.jl (used there as the C7 oracle), here as a plain `WordGraph` (not
# `circular(...)`) so both files can share it.
function _planar_braid_back_wordgraph()
    NP = DiagrammaticHecke.NodePort
    N = DiagrammaticHecke.Node
    nodes = [N(:braid, [1, 2], 3), N(:braid, [1, 2], 3)]
    edges = Edge[
        Edge(1, Leaf(1), NP(1, 3)), Edge(2, Leaf(2), NP(1, 2)), Edge(1, Leaf(3), NP(1, 1)),
        Edge(2, NP(1, 4), NP(2, 4)), Edge(1, NP(1, 5), NP(2, 3)), Edge(2, NP(1, 6), NP(2, 2)),
        Edge(1, Leaf(4), NP(2, 1)), Edge(2, Leaf(5), NP(2, 6)), Edge(1, Leaf(6), NP(2, 5)),
    ]
    return WordGraph(CircularWord([1, 2, 1, 1, 2, 1]), nodes, edges)
end

# C8-Fixtur (test/circular.jl "C8 _fr_braid_relation (R9)"): braid{2,1}–trivalent{1}–
# braid{1,2}, word [2,1,2,1,2,1,1] (the same node/edge topology as the fixture in
# test/circularweight.jl, ONLY the slot assignment differs). ⚠️ Leaf 7 is
# wired as colour 2, not colour 1 like the last letter of the word — that is NOT an
# error but the same peculiarity as in that fixture (the colour balance at the third
# node [3 slots per colour] admits no other choice for this topology: 2 leaves + bigon +
# trivalent edge would otherwise need 4 colour-1 slots where only 3 exist).
function _planar_braid_relation_wordgraph()
    NP = DiagrammaticHecke.NodePort
    N = DiagrammaticHecke.Node
    nodes = [N(:braid, [2, 1], 3), N(:trivalent, [1], 0), N(:braid, [1, 2], 3)]
    edges = Edge[
        Edge(2, Leaf(1), NP(1, 1)), Edge(1, Leaf(2), NP(1, 6)), Edge(2, Leaf(3), NP(1, 5)),
        Edge(1, Leaf(4), NP(2, 1)), Edge(2, Leaf(5), NP(3, 2)), Edge(1, Leaf(6), NP(3, 1)),
        Edge(2, Leaf(7), NP(3, 6)), Edge(1, NP(1, 4), NP(2, 2)), Edge(1, NP(2, 3), NP(3, 3)),
        Edge(1, NP(1, 2), NP(3, 5)), Edge(2, NP(1, 3), NP(3, 4)),
    ]
    return WordGraph(CircularWord([2, 1, 2, 1, 2, 1, 1]), nodes, edges)
end

# The C6 self-loop fixture (test/circular.jl, "slot counting at a braid SELF-LOOP"): braid
# [2,3,2,3,2,3] with dot(2) at slot 1 and self-loop(3) slot2<->slot6, leaves at slots
# 3/4/5 — wired BACKWARDS (leaf1→slot5, leaf2→slot4, leaf3→slot3) instead of forwards as
# in that fixture. Identical to the `gsl` fixture in test/circularregionrules.jl.
function _planar_braid_selfloop_fixture()
    NP = DiagrammaticHecke.NodePort
    return CircularGraph(CircularWord([2, 3, 2]),
        [circular_node([2, 3, 2, 3, 2, 3]), circular_node([2])],
        Edge[Edge(2, NP(1, 1), NP(2, 1)), Edge(3, NP(1, 2), NP(1, 6)),
             Edge(2, Leaf(1), NP(1, 5)), Edge(3, Leaf(2), NP(1, 4)),
             Edge(2, Leaf(3), NP(1, 3))])
end

