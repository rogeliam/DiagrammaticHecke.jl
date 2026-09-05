# CORE — the shared diagram/word layer under the circular pipeline.
# Faces.jl — cell detection.
#
# Idea (outside in, as the peel view suggests): between any two neighbouring letters
# of the boundary word there is a GAP, behind which a region of its own may extend.
# A dot merges its two neighbouring gaps; if several leaves run into the same
# braid/trivalent node, the sector between neighbouring legs becomes its own cell.
# Technically this is a face tracing over directed darts (a combinatorial map)
# starting at the boundary arcs — robust also for non-peelable graphs (e.g. two
# parallel strands without nodes, lens cells as in needle_a).
#
# ROTATION SYSTEM (derived from the existing conventions, no new data field needed):
#   * nodes: slots 1..deg are the cyclic order CLOCKWISE (trivalent 1,2,3;
#     braid 1..2m; dot trivial).
#   * boundary: leaves form the boundary polygon in word order; at leaf k the
#     clockwise order of outgoing darts is
#       [boundary arc k→k-1, graph edge, boundary arc k→k+1].
#     The strand sits BETWEEN the arcs because it points into the disc.
#
# ⚡ ORIENTATION — the house convention:
#
#   The BOUNDARY runs counter-clockwise (leaves in word order),
#   the ARMS at a node run clockwise (slot order).
#
# The advantage: at EVERY connection — node↔node as well as node↔boundary — the same
# rule holds, namely that one runs cw at one end and ccw at the other. With ccw
# throughout this would split into two rules (boundary nodes co-oriented, inner nodes
# counter-oriented). The choice is mathematically free (cw and ccw are mirror
# images); the uniformity is not.
#
# Consequences to check a wiring against:
#   * a node at several LEAVES: its slots run BACKWARDS along the leaf ring.
#   * two nodes joined by TWO edges (bigon): the two edges attach in OPPOSITE
#     directions.
#
# If boundary and nodes are NOT in this relation, the rotation system is inconsistent
# and the face walk runs on a surface of positive genus. Test: test/faces.jl, testset
# "the embedding is planar: V − E + F = 2".
#
# The boundary being a wall by itself is also why a separate, gap-based REGION notion
# is unnecessary: region == cell (see `region_count` below).
#
# FACE WALK: with σ = clockwise rotation, φ(d) = σ_target(d)(reverse(d)); the cycles
# of φ are the faces (face to the left of the dart). The outer face is the cycle
# through the forward boundary arcs (k→k+1) and is dropped; with an empty boundary
# (n = 0) every cycle is an inner cell. Free circles (Circle ports) have no vertices:
# each splits the region containing it and counts as one extra cell (its interior).
#
# The result is stored in the `cells` field of the WordGraph (see diagram/Graph.jl);
# the 3-argument constructor below is the convenient public entry point (it
# computes `cells` automatically).

# ---- the tracing core --------------------------------------------------------

struct _Dart
    src::Tuple{Symbol,Int}      # (:leaf,k) or (:node,i)
    dst::Tuple{Symbol,Int}
    rev::Int                    # id of the opposite dart
end

_vertex(p::Leaf) = (:leaf, p.k)
_vertex(p::NodePort) = (:node, p.node)
_vertex(p::Circle) = error("Circle ports have no vertex (free circles)")

"""
    WordGraph(word, nodes, edges)

The usual constructor: additionally computes the cells (`cells` field).
"""
function WordGraph(word::CircularWord, nodes::Vector{Node}, edges::Vector{Edge})
    _sc_audit(word, nodes, edges)          # audit mode only, see `slot_colour_audit`
    return WordGraph(word, nodes, edges, _trace_cells(word, nodes, edges))
end

# The tracer is GENERIC over the node type: the ONLY thing it needs from a node is
# its DEGREE. `deg` maps a node to its slot count (`_node_degree` for `Node`,
# `_circular_degree` for `CircularNode` — see circular/CircularFaces.jl). It is shared
# between the two diagram worlds (WordGraph and CircularGraph): the face walk is the
# subtlest code in the package (dart pairing, σ order at a leaf, outer face, circle
# counting), and two separate copies would silently drift apart.
function _trace_generic(word::CircularWord, nodes::AbstractVector,
                        edges::Vector{Edge}, deg)
    darts = _Dart[]
    slot_dart = Dict{Tuple{Int,Int}, Int}()   # (node, slot) → dart id
    leaf_dart = Dict{Int, Int}()              # leaf k → dart id of its edge
    port_dart = Dict{Port, Int}()
    ncircles = 0

    for e in edges
        if e.a isa Circle || e.b isa Circle
            ncircles += 1
            continue
        end
        v1, v2 = _vertex(e.a), _vertex(e.b)
        i, j = length(darts) + 1, length(darts) + 2
        push!(darts, _Dart(v1, v2, j))
        push!(darts, _Dart(v2, v1, i))
        port_dart[e.a] = i; port_dart[e.b] = j
        e.a isa NodePort && (slot_dart[(e.a.node, e.a.slot)] = i)
        e.b isa NodePort && (slot_dart[(e.b.node, e.b.slot)] = j)
        e.a isa Leaf && (leaf_dart[e.a.k] = i)
        e.b isa Leaf && (leaf_dart[e.b.k] = j)
    end

    n = length(word)
    # boundary arc: arc[k] = dart k → k+1 (cyclic); the reverse arc points outwards
    arc = Vector{Int}(undef, n)
    for k in 1:n
        i, j = length(darts) + 1, length(darts) + 2
        push!(darts, _Dart((:leaf, k), (:leaf, mod1(k + 1, n)), j))
        push!(darts, _Dart((:leaf, mod1(k + 1, n)), (:leaf, k), i))
        arc[k] = i
    end

    # rotation system σ
    σ = Dict{Tuple{Symbol,Int}, Vector{Int}}()
    for k in 1:n
        back = darts[arc[mod1(k - 1, n)]].rev   # dart k → k-1
        fwd  = arc[k]                            # dart k → k+1
        # The strand lies BETWEEN the boundary arcs, i.e. inside the disc, so the
        # order at a leaf is clockwise — the same sense as the slots at a node. See
        # the orientation block above.
        σ[(:leaf, k)] = haskey(leaf_dart, k) ? Int[back, leaf_dart[k], fwd] :
                                               Int[back, fwd]
    end
    for (ni, nd) in enumerate(nodes)
        lst = Int[]
        for s in 1:deg(nd)
            haskey(slot_dart, (ni, s)) && push!(lst, slot_dart[(ni, s)])
        end
        σ[(:node, ni)] = lst
    end

    # face walk: cycles of φ
    pos = Dict{Int, Int}()
    for lst in values(σ), (i, id) in enumerate(lst)
        pos[id] = i
    end
    nxt = Vector{Int}(undef, length(darts))
    for (id, d) in enumerate(darts)
        lst = σ[d.dst]
        nxt[id] = lst[mod1(pos[d.rev] + 1, length(lst))]
    end
    seen = falses(length(darts))
    cycles = Vector{Int}[]
    outer_face = 0
    # Outer face = the cycle of the BACKWARD boundary arcs. With the strand between
    # the arcs a forward arc turns inwards; it is the backward arcs that close up
    # around the disc (φ(arc[k].rev) = arc[k-1].rev).
    outer_dart = n >= 1 ? darts[arc[1]].rev : 0
    for id in 1:length(darts)
        seen[id] && continue
        cyc = Int[]
        j = id
        while !seen[j]
            seen[j] = true
            push!(cyc, j)
            j = nxt[j]
        end
        push!(cycles, cyc)
        if outer_dart != 0 && any(==(outer_dart), cyc)
            outer_face = length(cycles)
        end
    end

    # cell ids: cycles in order of their smallest dart id (canonical), circles last
    cell_of_face = zeros(Int, length(cycles))
    cid = 0
    for (fi, _) in enumerate(cycles)
        fi == outer_face && continue
        cell_of_face[fi] = (cid += 1)
    end
    face_of_dart = zeros(Int, length(darts))
    for (fi, cyc) in enumerate(cycles), d in cyc
        face_of_dart[d] = fi
    end

    return (darts = darts, slot_dart = slot_dart, port_dart = port_dart, arc = arc,
            cycles = cycles, outer_face = outer_face, cell_of_face = cell_of_face,
            face_of_dart = face_of_dart, ncircles = ncircles, ncells = cid + ncircles)
end

function _trace_cells_generic(word::CircularWord, nodes::AbstractVector,
                              edges::Vector{Edge}, deg)
    t = _trace_generic(word, nodes, edges, deg)
    darts, arc = t.darts, t.arc
    cell_of_face, face_of_dart = t.cell_of_face, t.face_of_dart
    n = length(word)

    # Gap k lies in the face of the INWARD-pointing boundary arc; with the cw leaf
    # order that is `arc[k]` itself.
    gap_cell = [cell_of_face[face_of_dart[arc[k]]] for k in 1:n]

    sector_cell = Vector{Vector{Int}}(undef, length(nodes))
    for (ni, nd) in enumerate(nodes)
        d = deg(nd)
        sc = zeros(Int, d)
        wired = [(s, t.slot_dart[(ni, s)]) for s in 1:d if haskey(t.slot_dart, (ni, s))]
        for (i, (s, _)) in enumerate(wired)
            d_next = wired[mod1(i + 1, length(wired))][2]
            # sector after slot s (cw) = face left of the dart at the next slot
            sc[s] = cell_of_face[face_of_dart[d_next]]
        end
        sector_cell[ni] = sc
    end

    port_cell = Dict{Port, NTuple{2, Int}}()
    for (p, id) in t.port_dart
        l = cell_of_face[face_of_dart[id]]
        r = cell_of_face[face_of_dart[darts[id].rev]]
        port_cell[p] = (l, r)
    end

    return Cells(t.ncells, gap_cell, sector_cell, port_cell)
end

# WordGraph-specific wrappers — all existing call sites stay unchanged.
_trace(word::CircularWord, nodes::Vector{Node}, edges::Vector{Edge}) =
    _trace_generic(word, nodes, edges, _node_degree)
_trace_cells(word::CircularWord, nodes::Vector{Node}, edges::Vector{Edge}) =
    _trace_cells_generic(word, nodes, edges, _node_degree)

# face walks as vertex sequences (for label anchors in render/Cells.jl):
# (walks, cell_of_face, ncells, ncircles); walk = src vertices of the cycle's darts.
function _face_walks_generic(word::CircularWord, nodes::AbstractVector,
                            edges::Vector{Edge}, deg)
    t = _trace_generic(word, nodes, edges, deg)
    walks = [[t.darts[d].src for d in cyc] for cyc in t.cycles]
    return walks, t.cell_of_face, t.ncells, t.ncircles
end

_face_walks(g::WordGraph) = _face_walks_generic(g.word, g.nodes, g.edges, _node_degree)

# ---- the API -----------------------------------------------------------------

"""
    Cell

An inner cell: `id`, the boundary gaps it touches (`gaps`) and the node sectors
(`sectors` = pairs `(node, slot)`; the sector lies clockwise after the slot). Cells
with neither gaps nor sectors are the interiors of free circles.
"""
struct Cell
    id::Int
    gaps::Vector{Int}
    sectors::Vector{Tuple{Int, Int}}
end

"""
    face_count(g::WordGraph) -> Int

Number of inner cells of the diagram (a stored attribute, O(1)).
"""
face_count(g::WordGraph) = g.cells.ncells

"""
    inner_faces(g::WordGraph) -> Vector{Cell}

The inner cells of the diagram, canonically numbered (smallest dart id of the face
cycle, free circles last).
"""
inner_faces(g::WordGraph) = _inner_faces_from(g.cells)

# The three bodies below read ONLY `Cells` — factored out so that CircularGraph
# (circular/CircularFaces.jl) shares exactly this logic.
function _inner_faces_from(cs::Cells)
    out = [Cell(i, Int[], Tuple{Int,Int}[]) for i in 1:cs.ncells]
    for (k, c) in enumerate(cs.gap_cell)
        c >= 1 && push!(out[c].gaps, k)
    end
    for (ni, scs) in enumerate(cs.sector_cell), (s, c) in enumerate(scs)
        c >= 1 && push!(out[c].sectors, (ni, s))
    end
    return out
end

"""
    face_of_port(g::WordGraph, p::Port; side = :left) -> Union{Int, Nothing}

Cell id on the side `side` (`:left` / `:right`) of the edge at port `p`; `nothing`
if the outer face lies there.
"""
face_of_port(g::WordGraph, p::Port; side::Symbol = :left) =
    _face_of_port(g.cells, p; side = side)

function _face_of_port(cs::Cells, p::Port; side::Symbol = :left)
    haskey(cs.port_cell, p) ||
        throw(ArgumentError("port $p is not attached to any edge"))
    l, r = cs.port_cell[p]
    side === :left  && return l == 0 ? nothing : l
    side === :right && return r == 0 ? nothing : r
    throw(ArgumentError("side must be :left or :right, not :$side"))
end

"""
    gap_cell(g::WordGraph, k::Int) -> Int

Cell id behind the gap between leaf `k` and `k+1` (boundary arc k→k+1).
"""
gap_cell(g::WordGraph, k::Int) = g.cells.gap_cell[k]

"""
    sector_cell(g::WordGraph, node::Int, slot::Int) -> Int

Cell id of the sector clockwise after `slot` at node `node`.
"""
sector_cell(g::WordGraph, node::Int, slot::Int) = g.cells.sector_cell[node][slot]

# A MorphismGraph is just a WordGraph with two cuts, and cells do not depend on the
# cuts — so simply forward.
face_count(m::MorphismGraph) = face_count(m.graph)
inner_faces(m::MorphismGraph) = inner_faces(m.graph)
face_of_port(m::MorphismGraph, p::Port; side::Symbol = :left) =
    face_of_port(m.graph, p; side = side)

"""
    show_faces(g::WordGraph) -> String

Compact console dump of the inner cells (gaps + sectors).
"""
show_faces(g::WordGraph) = _show_faces(g.cells)

function _show_faces(cs::Cells)
    io = IOBuffer()
    println(io, "cells: ", cs.ncells)
    for c in _inner_faces_from(cs)
        gaps = isempty(c.gaps) ? "—" : join(c.gaps, ",")
        secs = isempty(c.sectors) ? "—" :
               join(("n$n.$s" for (n, s) in c.sectors), ",")
        println(io, "  [", c.id, "] gaps {", gaps, "}  sectors {", secs, "}")
    end
    return String(take!(io))
end

# ---- gap groups along a SOLID boundary --------------------------------------
#
# Since the boundary separates in the face walk itself, region == cell and
# `region_count == face_count`. The grouping helpers below compute the same
# boundary regions independently, via a union-find over boundary gaps.
#
# THE RULE, in two parts:
#
#   (a) AT THE BOUNDARY: neighbouring gaps merge if the leaf between them is not a
#       wall. The edge at leaf k is a wall, except when
#         * it ends in a 1-armed node (dot) — the strand stops there;
#         * the leaf has no edge at all — then nothing separates there.
#       A leaf–leaf chord IS a wall (both ends are leaves, no dot).
#
#   (b) THROUGH THE INTERIOR: two gaps of the same face can also be connected
#       without being neighbours on the boundary. What matters is WHAT the inner
#       path runs through: through a NODE of degree >= 2 it separates the two gaps;
#       over leaf–leaf chords only they belong together. A chord merely cuts off an
#       ear, whereas a real node with its >= 2 arms genuinely divides the disc.
#
# Part (b) is what a pure boundary walk misses. The case that showed it:
# `circular_double_leaf([1,2,1], [1,0,1], [1,2,1], [1,0,1])`, two chords `Leaf1—Leaf3`
# and `Leaf6—Leaf4` — face 2 has gaps {3,6} and NO sectors, i.e. it is the band
# between the chords, so gaps 3 and 6 are its two ends and one region inside; the
# boundary walk counted 4, correct is 3. Compare `[1,0,0]`: there face 2 has the SAME
# gaps {3,6}, but the path runs through the trivalents n2/n4 ⇒ separated ⇒ 4.
# Combinatorially both faces are the same cycle type; only the node rule tells them
# apart.

"true if the edge at leaf `k` does NOT interrupt the boundary there"
function _leaf_is_transparent(nodes::AbstractVector, edges::Vector{Edge}, k::Int, deg)
    for e in edges
        for (p, q) in ((e.a, e.b), (e.b, e.a))
            if p isa Leaf && p.k == k
                # dot (one arm) => transparent; any larger node => wall.
                q isa NodePort && return deg(nodes[q.node]) == 1
                return false            # leaf–leaf chord: wall
            end
        end
    end
    return true                         # no edge at all: nothing separates here
end

"""
    _rim_gap_groups(word, nodes, edges, deg) -> Vector{Int}

Group id per boundary gap: `out[k] == out[l]` means gaps `k` and `l` belong to the
same region along a solid boundary. Union-find over the `n` gaps, merged (a) across
every transparent leaf and (b) across every inner path meeting no node of degree
>= 2 — see the comment block above.
"""
function _rim_gap_groups(word::CircularWord, nodes::AbstractVector,
                         edges::Vector{Edge}, deg)
    n = length(word)
    n == 0 && return Int[]
    parent = collect(1:n)
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    merge_gaps(a, b) = ((a, b) = (find(a), find(b)); a != b && (parent[a] = b))

    # (a) at the boundary: leaf k sits between gap k-1 (arc k-1→k) and gap k.
    for k in 1:n
        _leaf_is_transparent(nodes, edges, k, deg) && merge_gaps(mod1(k - 1, n), k)
    end

    # (b) through the interior: along each inner face cycle. The INWARD-pointing
    # boundary arcs (`arc[k].rev`) mark the gaps; two consecutive ones in the cycle
    # merge if the path between them visits no node of degree >= 2.
    t = _trace_generic(word, nodes, edges, deg)
    gap_of_dart = Dict(t.darts[t.arc[k]].rev => k for k in 1:n)
    for (fi, cyc) in enumerate(t.cycles)
        fi == t.outer_face && continue
        L = length(cyc)
        at = [i for i in 1:L if haskey(gap_of_dart, cyc[i])]
        length(at) < 2 && continue
        for a in 1:length(at)
            i, j = at[a], at[mod1(a + 1, length(at))]
            blocked = false
            p = mod1(i + 1, L)
            while p != j
                v = t.darts[cyc[p]].src
                if v[1] === :node && deg(nodes[v[2]]) >= 2
                    blocked = true
                    break
                end
                p = mod1(p + 1, L)
            end
            blocked || merge_gaps(gap_of_dart[cyc[i]], gap_of_dart[cyc[j]])
        end
    end

    return [find(k) for k in 1:n]
end

"""
    region_count(g) -> Int

Number of regions — **identical to `face_count(g)`**, because a region is exactly
one cell (see `Region`). The name is kept so the call sites (label vectors,
distances) keep working.
"""
function _region_count(word::CircularWord, nodes::AbstractVector,
                       edges::Vector{Edge}, deg, cs::Cells)
    return cs.ncells
end

region_count(g::WordGraph) =
    _region_count(g.word, g.nodes, g.edges, _node_degree, g.cells)
region_count(m::MorphismGraph) = region_count(m.graph)

"""
    boundary_regions(g) -> Vector{Vector{Int}}

The boundary regions as lists of gap indices, in order of their smallest gap.
Together with the gapless inner cells they make up `region_count(g)`.
"""
function _boundary_regions(word::CircularWord, nodes::AbstractVector,
                           edges::Vector{Edge}, deg)
    n = length(word)
    n == 0 && return Vector{Int}[]
    # The boundary separates in the face walk itself (cw leaf order, see the
    # orientation block above). Two gaps belong together exactly when they lie in
    # the SAME CELL, read directly from `cs.gap_cell`.
    cs = _trace_cells_generic(word, nodes, edges, deg)
    out = Dict{Int, Vector{Int}}()
    for k in 1:n
        push!(get!(out, cs.gap_cell[k], Int[]), k)
    end
    return sort(collect(values(out)), by = first)
end

boundary_regions(g::WordGraph) =
    _boundary_regions(g.word, g.nodes, g.edges, _node_degree)
boundary_regions(m::MorphismGraph) = boundary_regions(m.graph)

"""
    Region

**A region is EXACTLY ONE CELL**: only faces are used — every
dart and every leaf-to-leaf connection is a possible separation. `id == cell`,
and `gaps` are the boundary gaps this cell touches.

The type is kept so the call sites (`R.gaps`, `R.cell`, `regions(g)`) keep working —
a thin view on `Cell`, not a refinement.

The cell assignment comes out of the tracer without heuristics (`face_of_dart` ->
`cell_of_face`); the boundary arcs are genuine darts there.
"""
struct Region
    id::Int
    gaps::Vector{Int}
    cell::Int
end

"""
    regions(g) -> Vector{Region}

The regions, numbered like the cells: **one region per cell**, `id == cell`, so
`length(regions(g)) == region_count(g) == face_count(g)`.
"""
function _regions(word::CircularWord, nodes::AbstractVector,
                  edges::Vector{Edge}, deg, cs::Cells)
    return [Region(c.id, c.gaps, c.id) for c in _inner_faces_from(cs)]
end

regions(g::WordGraph) = _regions(g.word, g.nodes, g.edges, _node_degree, g.cells)
regions(m::MorphismGraph) = regions(m.graph)
