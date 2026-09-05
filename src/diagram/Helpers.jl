# diagram/Helpers.jl — incidence / rewiring helpers on a WordGraph.
# Shared core: used by the circular rules (via ZamoRules) and by morphism/.

# Do two ports refer to the same endpoint?
_same_port(a::Leaf, b::Leaf) = a.k == b.k
_same_port(a::NodePort, b::NodePort) = a.node == b.node && a.slot == b.slot
_same_port(::Port, ::Port) = false

# The edges of `g` incident to node index `ni` (as (edge_index, the-port-on-ni,
# the-other-port)).
function _edges_at_node(g::WordGraph, ni::Int)
    hits = Tuple{Int,Port,Port}[]
    for (ei, e) in enumerate(g.edges)
        if e.a isa NodePort && e.a.node == ni
            push!(hits, (ei, e.a, e.b))
        elseif e.b isa NodePort && e.b.node == ni
            push!(hits, (ei, e.b, e.a))
        end
    end
    return hits
end

# Rebuild `g` with the nodes at indices `dead` removed. Every surviving NodePort is
# renumbered to the new node indexing; edges given in `keep_edges` (already using
# OLD node indices) are renumbered too. Leaves are untouched. Returns a WordGraph.
function _delete_nodes(g::WordGraph, dead::Set{Int}, keep_edges::Vector{Edge})
    # old index -> new index for surviving nodes
    remap = Dict{Int,Int}()
    newnodes = Node[]
    for (ni, nd) in enumerate(g.nodes)
        ni in dead && continue
        push!(newnodes, nd)
        remap[ni] = length(newnodes)
    end
    fixport(p::Leaf) = p
    fixport(p::Circle) = p
    fixport(p::NodePort) = NodePort(remap[p.node], p.slot)
    newedges = [Edge(e.colour, fixport(e.a), fixport(e.b)) for e in keep_edges]
    return WordGraph(g.word, newnodes, newedges)
end

