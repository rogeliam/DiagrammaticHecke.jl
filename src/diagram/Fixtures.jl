# diagram/Fixtures.jl — small named example graphs (the exact structures used
# throughout tests and notebooks): diagram data, not drawing code.

"unit: a trivalent whose third edge ends in a dot; the other two go to boundary 11."
function unit_graph()
    WordGraph(CircularWord([1, 1]),
        [Node(:trivalent, [1], 0), Node(:dot, [1], 0)],
        [Edge(1, Leaf(1), NodePort(1, 1)),
         Edge(1, Leaf(2), NodePort(1, 2)),
         Edge(1, NodePort(1, 3), NodePort(2, 1))])
end

"needle A: two trivalents joined by TWO edges; the remaining leg of each → boundary 11. = 0."
function needle_a_graph()
    WordGraph(CircularWord([1, 1]),
        [Node(:trivalent, [1], 0), Node(:trivalent, [1], 0)],
        # The two edges of the bigon attach in OPPOSITE slot order — that is
        # geometry, valid in either orientation. Wired parallel instead,
        # face_count came out 2 instead of 3.
        [Edge(1, NodePort(1, 1), NodePort(2, 2)),      # two internal joins
         Edge(1, NodePort(1, 2), NodePort(2, 1)),
         Edge(1, Leaf(1), NodePort(1, 3)),             # one leg each to the boundary
         Edge(1, Leaf(2), NodePort(2, 3))])
end

"needle B: one trivalent, two of its edges short-circuited (self-loop); third → boundary 1. = 0."
function needle_b_graph()
    WordGraph(CircularWord([1]),
        [Node(:trivalent, [1], 0)],
        [Edge(1, NodePort(1, 1), NodePort(1, 2)),      # self-loop on slots 1,2
         Edge(1, Leaf(1), NodePort(1, 3))])
end

"R3 braid(1,3) with a dot on a 1-slot → boundary 133 (or on a 3-slot → 113)."
function r3_braid_dot(dotcolour::Int = 1)
    # braid slots 1,2,3,4 carry colours 1,3,1,3. cap one slot of `dotcolour`.
    slotcol = [1, 3, 1, 3]
    cap = findfirst(==(dotcolour), slotcol)
    nodes = [Node(:braid, [1, 3], 2), Node(:dot, [dotcolour], 0)]
    edges = Edge[Edge(dotcolour, NodePort(1, cap), NodePort(2, 1))]
    # leaf ring runs ccw against the cw arm order: walk the slots backwards.
    lk = 0; word = Int[]
    for s in reverse(1:4)
        s == cap && continue
        lk += 1
        push!(word, slotcol[s])
        push!(edges, Edge(slotcol[s], Leaf(lk), NodePort(1, s)))
    end
    WordGraph(CircularWord(word), nodes, edges)
end

"R4 braid(1,2;m=3) with a dot on a 1-slot → boundary 21212 (or 2-slot → 12121)."
function r4_braid_dot(dotcolour::Int = 1)
    slotcol = [1, 2, 1, 2, 1, 2]
    cap = findfirst(==(dotcolour), slotcol)
    nodes = [Node(:braid, [1, 2], 3), Node(:dot, [dotcolour], 0)]
    edges = Edge[Edge(dotcolour, NodePort(1, cap), NodePort(2, 1))]
    # leaf ring runs ccw against the cw arm order: walk the slots backwards.
    lk = 0; word = Int[]
    for s in reverse(1:6)
        s == cap && continue
        lk += 1
        push!(word, slotcol[s])
        push!(edges, Edge(slotcol[s], Leaf(lk), NodePort(1, s)))
    end
    WordGraph(CircularWord(word), nodes, edges)
end

"barbell: two dots joined by one edge, no boundary. = α_i · (empty diagram)."
barbell_graph() = WordGraph(EMPTY,
    [Node(:dot, [1], 0), Node(:dot, [1], 0)],
    [Edge(1, NodePort(1, 1), NodePort(2, 1))])
