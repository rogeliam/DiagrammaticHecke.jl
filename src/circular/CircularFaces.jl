# circular/CircularFaces.jl — cells for CircularGraph.
#
# NO tracer of its own: the face walk from diagram/Faces.jl is generic over the
# node type (`_trace_generic` / `_trace_cells_generic` take a degree function),
# here it is just called with `_circular_degree`. This way the two diagram worlds
# cannot drift apart — the face walk is the most subtle code in the package,
# and a second copy would silently drift.
#
# The only precondition for this: the slots `1..deg` of a `CircularNode` ARE its
# clockwise cyclic order (as with the plain `Node`). That holds by construction,
# see circular/CircularGraph.jl.

"""
    CircularGraph(word, nodes, edges)

The usual constructor: additionally computes the cells (`cells` field), via
the same generic tracer as `WordGraph`.
"""
function CircularGraph(word::CircularWord, nodes::Vector{CircularNode}, edges::Vector{Edge})
    _sc_audit(word, nodes, edges)          # only in audit mode, see `slot_colour_audit`
    return CircularGraph(word, nodes, edges,
                    _trace_cells_generic(word, nodes, edges, _circular_degree))
end

# The cell API — the same bodies as for WordGraph, just via `g.cells`.
face_count(g::CircularGraph) = g.cells.ncells
inner_faces(g::CircularGraph) = _inner_faces_from(g.cells)
face_of_port(g::CircularGraph, p::Port; side::Symbol = :left) =
    _face_of_port(g.cells, p; side = side)
gap_cell(g::CircularGraph, k::Int) = g.cells.gap_cell[k]
sector_cell(g::CircularGraph, node::Int, slot::Int) = g.cells.sector_cell[node][slot]
show_faces(g::CircularGraph) = _show_faces(g.cells)

# Regions under a solid boundary (a second, parallel cell notion — see the
# comment block in diagram/Faces.jl). In the circular picture a dot is a 1-armed
# node, `_circular_degree` gives that directly.
region_count(g::CircularGraph) =
    _region_count(g.word, g.nodes, g.edges, _circular_degree, g.cells)
boundary_regions(g::CircularGraph) =
    _boundary_regions(g.word, g.nodes, g.edges, _circular_degree)
regions(g::CircularGraph) = _regions(g.word, g.nodes, g.edges, _circular_degree, g.cells)

# Euler test and convention test — the same bodies as for WordGraph, just with
# `_circular_degree`. Defined in diagram/Faces.jl; only the type method sits here,
# so `euler`/`check_wiring` compute the same thing in both diagram worlds (no
# second code path, as with the tracer already).
euler(g::CircularGraph) = _euler_generic(g.word, g.nodes, g.edges, _circular_degree)
check_wiring(g::CircularGraph) = _check_wiring_generic(g.word, g.nodes, g.edges, _circular_degree)

# Face walks as vertex sequences (label anchors for the renderer) — the same
# `_face_walks_generic` core as `_face_walks` (diagram/Faces.jl), only the degree
# function differs.
_circular_face_walks(g::CircularGraph) =
    _face_walks_generic(g.word, g.nodes, g.edges, _circular_degree)
