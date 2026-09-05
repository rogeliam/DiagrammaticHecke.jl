# circular/CircularBoundaryWord.jl — the BOUNDARY WORD of a region.
#
# WHY. The C rules match on node shape (`kind`/`arms`) — exactly what merging
# makes variable: the same local situation looks different depending on how
# far the surrounding nodes have merged, and every rule needs case splits per
# variant (duplicated predicates, unreachable guards). A region's BOUNDARY WORD — the cyclic
# sequence of its boundary edges with the corners between them — is instead
# INVARIANT under merging: a bigon stays a bigon no matter how far its
# corners have merged into circular nodes. Rules that match on boundary words
# (circular/rules/CircularRegionRules.jl) inherit that invariance for free.
#
# DATA SOURCE is the verified face tracer (diagram/Faces.jl, §ROTATIONSSYSTEM:
# boundary ccw, arms cw; region == cell): `_circular_dart_region`
# (circular/CircularRegion.jl) gives face cycles built from darts; here they are only
# translated into edge/corner view. NO second tracer.
#
# READING a cycle: letter j = the j-th dart of the face walk (the edge, seen
# from the region), corner j = the vertex BETWEEN dart j and dart j+1. A
# dot-capped arm inside a region is NOT a break: the face walk runs into the
# excursion (v -> dot -> v), so the capping edge appears TWICE in the boundary
# word (once per direction), with a corner at the dot in between. That is
# exactly how a region rule recognises "a dot hangs in this region" (pitchfork,
# CircularRegionRules.jl).
#
# FREE CIRCLES (circle ports) are invisible to the tracer (Faces.jl): their
# cells have NO face cycle and get an empty boundary word here.

"""
    CircularCorner

A CORNER in the boundary word of a region: the vertex between two consecutive
boundary edges. `vertex` is `(:node, v)` or `(:leaf, k)` (tracer convention,
diagram/Faces.jl). At a node, `slot_in`/`slot_out` are the two slots the face
walk enters/leaves through (`0` at a leaf). For a fully wired graph
(`is_wired`), `slot_out` is the cyclically next slot after `slot_in`,
clockwise.
"""
struct CircularCorner
    vertex::Tuple{Symbol,Int}
    slot_in::Int
    slot_out::Int
end

"""
    CircularBoundaryLetter

A LETTER of the boundary word: an edge on the region's boundary, plus the
corner AFTER it (see the file header for the reading convention). `edge` is
the index into `g.edges`; `edge == 0` (and `colour == 0`) means a BOUNDARY ARC
of the diagram — the region touches the diagram boundary there, between two
leaves.
"""
struct CircularBoundaryLetter
    edge::Int
    colour::Int
    corner::CircularCorner
end

"""
    CircularBoundaryWord

The boundary word of region `region`: the cyclic letter sequence in face-walk
order (region left of the dart). Cyclic — rotations of the same word describe
the same region; the matcher (CircularRegionRules.jl) tries all rotations, NO
reflection (same convention as `circular_canonical_key`: equality holds only up to
rotation, `flip` is a separate operation).
"""
struct CircularBoundaryWord
    region::Int
    letters::Vector{CircularBoundaryLetter}
end

Base.length(w::CircularBoundaryWord) = length(w.letters)

"""
    circular_boundary_words(g::CircularGraph) -> Vector{CircularBoundaryWord}

The boundary words of ALL regions of `g`, indexed like `regions(g)`
(`words[R].region == R`). One tracer run (`_circular_dart_region`), then per inner
face cycle the translation dart -> (edge, corner). Cells of free circles have
no cycle and get an empty boundary word (see file header).
"""
function circular_boundary_words(g::CircularGraph)
    dart_region, t = _circular_dart_region(g)

    # Edge per dart (0 = boundary arc). Circle edges have no darts.
    edge_of_dart = zeros(Int, length(t.darts))
    for (ei, e) in enumerate(g.edges)
        (e.a isa Circle || e.b isa Circle) && continue
        d = t.port_dart[e.a]
        edge_of_dart[d] = ei
        edge_of_dart[t.darts[d].rev] = ei
    end
    # Slot per node-dart (a dart leaves (:node,ni) through this slot).
    slot_of_dart = Dict{Int,Int}()
    for ((_, s), d) in t.slot_dart
        slot_of_dart[d] = s
    end

    nregs = region_count(g)
    words = [CircularBoundaryWord(R, CircularBoundaryLetter[]) for R in 1:nregs]
    for (fi, cyc) in enumerate(t.cycles)
        fi == t.outer_face && continue
        R = t.cell_of_face[fi]
        letters = CircularBoundaryLetter[]
        for (j, d) in enumerate(cyc)
            dn = cyc[mod1(j + 1, length(cyc))]
            ei = edge_of_dart[d]
            col = ei == 0 ? 0 : g.edges[ei].colour
            vtx = t.darts[d].dst
            slot_in  = get(slot_of_dart, t.darts[d].rev, 0)
            slot_out = get(slot_of_dart, dn, 0)
            push!(letters, CircularBoundaryLetter(ei, col, CircularCorner(vtx, slot_in, slot_out)))
        end
        words[R] = CircularBoundaryWord(R, letters)
    end
    return words
end

"""
    circular_region_boundary_word(g::CircularGraph, R::Int) -> CircularBoundaryWord

The boundary word of ONE region — a convenience wrapper around
[`circular_boundary_words`](@ref) (still runs the full tracer; the graphs are
small).
"""
function circular_region_boundary_word(g::CircularGraph, R::Int)
    ws = circular_boundary_words(g)
    (1 <= R <= length(ws)) || error("circular_region_boundary_word: region $R does not exist")
    return ws[R]
end

"""
    is_interior(w::CircularBoundaryWord) -> Bool

An INTERIOR region: no letter is a boundary arc, no corner is a leaf. The
local region rules (CircularRegionRules.jl) only match interior regions.
"""
is_interior(w::CircularBoundaryWord) =
    all(l -> l.edge != 0 && l.corner.vertex[1] === :node, w.letters)

# Console dump in the style of `show_circulargraph`.
function Base.show(io::IO, w::CircularBoundaryWord)
    parts = String[]
    for l in w.letters
        e = l.edge == 0 ? "∂" : "e$(l.edge)(c$(l.colour))"
        v = l.corner.vertex[1] === :node ? "n$(l.corner.vertex[2])" : "leaf$(l.corner.vertex[2])"
        push!(parts, "$e→$v")
    end
    print(io, "CircularBoundaryWord(R", w.region, ": ", join(parts, " "), ")")
end
