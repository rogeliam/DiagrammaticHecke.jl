# circular/CircularGraphDisplay.jl — the wiring check (is_wired), console display
# (show, show_circulargraph) and the named 4n-cluster fixture
# (general_braid_cluster, general_braid_cluster_swapped). Included after
# CircularGraphKey.jl.
#
# ---- wiring check ---------------------------------------------------------------

"""
    is_wired(g::CircularGraph) -> Bool

Every slot `1..deg` of every node is occupied exactly once, and the edge
colour matches the arm colour at both ends. Analogue of `is_planar_braid`, but
STRONGER: with an open arm count, slot completeness plus colour fidelity is
the entire planarity contract a `CircularNode` can even enter into.
"""
function is_wired(g::CircularGraph)
    used = [falses(arm_count(nd)) for nd in g.nodes]
    for e in g.edges
        for p in (e.a, e.b)
            if p isa NodePort
                nd = g.nodes[p.node]
                (1 <= p.slot <= arm_count(nd)) || return false
                used[p.node][p.slot] && return false        # occupied twice
                used[p.node][p.slot] = true
                arm_colour(nd, p.slot) == e.colour || return false
            elseif p isa Circle
                p.colour == e.colour || return false
            end
            # A leaf carries no colour of its own in the port — just like
            # `is_planar_braid`, we don't check leaves against the word.
        end
    end
    return all(all(u) for u in used)
end

# ---- display ----------------------------------------------------------------

Base.show(io::IO, nd::CircularNode) =
    print(io, "CircularNode(:", nd.kind, ", [", join(nd.arms, ","), "])")

Base.show(io::IO, g::CircularGraph) =
    print(io, "CircularGraph(", isempty(letters(g.word)) ? "ε" : join(letters(g.word)),
          ", ", length(g.nodes), " nodes, ", length(g.edges), " edges)")

"""
    show_circulargraph(g::CircularGraph; color = true) -> String

Console dump analogous to `show_wordgraph` (diagram/Graph.jl). Reuses
`_port_str` and `_wg_paint` from there — ports are identical between the two
worlds.
"""
function show_circulargraph(g::CircularGraph; color::Bool = true)
    io = IOBuffer()
    lets = letters(g.word)
    println(io, "CircularGraph  (boundary length ", length(lets), ")")
    print(io, "  boundary: ")
    if isempty(lets)
        print(io, "ε")
    else
        for (k, c) in enumerate(lets)
            print(io, _wg_paint(c, string(c); on = color))
            k < length(lets) && print(io, " ")
        end
    end
    println(io)
    println(io, "  nodes:")
    isempty(g.nodes) && println(io, "    (none)")
    for (idx, nd) in enumerate(g.nodes)
        println(io, "    n", idx, " = ", nd.kind, "(", join(nd.arms, ","), ")",
                "   [deg ", arm_count(nd), "]")
    end
    println(io, "  edges:")
    isempty(g.edges) && println(io, "    (none)")
    for e in g.edges
        line = "    " * _port_str(e.a) * " ── " * _port_str(e.b)
        println(io, _wg_paint(e.colour, line; on = color), "   (colour ",
                _wg_paint(e.colour, string(e.colour); on = color), ")")
    end
    return String(take!(io))
end

# ---- the `4n` CLUSTER as a named fixture -------------------------------------
#
# This figure lives ONCE here; notebooks and tests fetch it from here.

"""
    general_braid_cluster(n::Int; s = 1, t = 2) -> CircularGraph

The **`4n` cluster**: `n` `st`-braids in a circle, neighbouring braids each
connected via a `t`-trivalent, the inner `s`-arms meeting at the centre.
Boundary word `(s t s t)^n`.

Nodes: `1..n` the braids, `n+1..2n` the trivalents, `2n+1` the centre (an
`n`-armed `s`-node; for `n = 3` itself a plain trivalent, for `n = 2` it is
absent — the two inner `s`-arms are connected directly).

KEPT AS AN EXAMPLE ONLY: this term is always simplified by the general
parallel rule, so it needs no dedicated rule of its own. The figure stays
because it is a good demonstration that such a diagram DOES simplify.
`euler(g) == 2` and `check_wiring(g)` empty for every `n ≥ 2` — anyone
building their own copy should check exactly that.
"""
function general_braid_cluster(n::Int; s::Int = 1, t::Int = 2)
    n >= 2 || throw(ArgumentError("general_braid_cluster: n >= 2, got $n"))
    word  = CircularWord(repeat([s, t, s, t], n))
    nodes = CircularNode[]
    for _ in 1:n; push!(nodes, circular_node([s, t, s, t, s, t])); end
    for _ in 1:n; push!(nodes, circular_node([t, t, t])); end
    n >= 3 && push!(nodes, circular_node(fill(s, n)))
    Z = 2n + 1
    edges = Edge[]
    for i in 1:n
        a = 4 * (i - 1) + 1
        push!(edges, Edge(s, Leaf(a + 2), NodePort(i, 1)))
        push!(edges, Edge(t, Leaf(a + 1), NodePort(i, 2)))
        push!(edges, Edge(s, Leaf(a),     NodePort(i, 3)))
        ip = mod1(i - 1, n)
        push!(edges, Edge(t, NodePort(i, 4), NodePort(n + ip, 2)))
        push!(edges, Edge(t, NodePort(i, 6), NodePort(n + i, 1)))
        push!(edges, Edge(t, Leaf(a + 3), NodePort(n + i, 3)))
    end
    if n == 2
        push!(edges, Edge(s, NodePort(1, 5), NodePort(2, 5)))
    else
        for i in 1:n
            push!(edges, Edge(s, NodePort(i, 5), NodePort(Z, n + 1 - i)))
        end
    end
    return CircularGraph(word, nodes, edges)
end

"""
    general_braid_cluster_swapped(n::Int; s = 1, t = 2, shift = 1) -> CircularGraph

The colour-swapped sibling: centre coloured `t`, the connecting trivalents
coloured `s`, same boundary word `(st)^{2n}`.
"""
function general_braid_cluster_swapped(n::Int; s::Int = 1, t::Int = 2, shift::Int = 1)
    g = general_braid_cluster(n; s = t, t = s)
    N = 4n
    rot(p) = p isa Leaf ? Leaf(mod1(p.k + shift, N)) : p
    edges = Edge[Edge(e.colour, rot(e.a), rot(e.b)) for e in g.edges]
    word  = CircularWord([isodd(k) ? s : t for k in 1:N])
    return CircularGraph(word, g.nodes, edges)
end
