# render/CircularRegionTree.jl — the DUAL GRAPH overlay: the region-word BFS tree.
#
# An optional dual graph: it goes from the 0-region to all others along a tree and
# writes the next reflection on each edge.
#
# WHAT IS DRAWN. One dot per REACHABLE region, sitting exactly on the anchor the
# region NUMBER uses (`_circular_region_anchors`, render/CircularCells.jl) — so
# the overlay and `cell_numbers`/`cell_labels = :words` are guaranteed to agree
# on where a region is. One arrow per tree edge (parent region → child region),
# taken from the `parent_edge` vector of `circular_region_words`
# (circular/rules/CircularRegionWord.jl); the root (the MARKED region) has
# `parent_edge == 0` and gets a ring instead of an incoming arrow.
#
# THE EDGE LABEL is the SIMPLE reflection `s` — the colour of the crossed edge
# (`g.edges[parent_edge[R]].colour`). That is
# exactly the letter the region word appends at this step, so reading the arrow
# labels along the tree path from the mark to `R` spells out `words[R]`. The
# conjugated wall reflection `w·s·w⁻¹` is NOT drawn: it would need a normal
# form to display and says nothing the path does not already say.
#
# WHY IT IS AN OVERLAY, not its own picture: the whole point is to see WHICH
# regions of the drawn diagram the tree runs through. It is printed LAST in
# `circular_tutte_svg`, so it sits on top of the strands.
#
# THE MARKING IS REQUIRED. Region words exist only relative to the black mark,
# so `region_parent` is filled by `_circular_tutte_morphism` — a bare
# `CircularGraph` draws no tree (`region_parent === nothing`), the same rule as
# for `:distances`.

"The colour palette of the diagram, reused for the arrow labels (colour = letter)."
_crt_colour(c::Int) = c == 1 ? "#c0392b" : c == 2 ? "#2471a3" : c == 3 ? "#1e8449" : "#555"

"""
    _circular_region_tree_svg(g, parent, words; size, r, mark_between, …) -> String

The SVG fragment of the dual-graph overlay (file head): a dot per reachable
region at its `_circular_region_anchors` position, an arrow per tree edge
`parent[R] != 0` from the region on the OTHER side of that edge to `R`, labelled
with the edge's colour.

`parent`, `words` and `letters` come from
[`circular_region_tree_parents`](@ref): `parent[R]` is the REGION the tree step
into `R` came from and `letters[R]` the letters it appended. `parent[R] == 0`
holds for the root AND for unreachable regions — `words[R] === nothing` is what
tells them apart, so an unreachable region gets no dot at all.

A step may cross a NODE rather than an edge, and then it appends several
letters at once and there is no edge to point at. That is why the arrow runs
region to region and carries the step word as its label: at a `222222` node
every sector is one step from the marked region, so the overlay draws a STAR
out of the root, not a chain. The layout arguments MUST be the same
ones `circular_tutte_svg` draws with — otherwise the dots land in a different
embedding than the diagram (the trap documented in render/CircularCells.jl).
"""
function _circular_region_tree_svg(g::CircularGraph, parent::Vector{Int},
                                   words::Union{Nothing,Vector} = nothing;
                                   letters::Union{Nothing,Vector} = nothing,
                                   size::Int, r::Real,
                                   mark_between::MarkSpec = nothing,
                                   mark_between2::MarkSpec = nothing,
                                   bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                                   top_leaves::Union{Nothing,Vector{Int}} = nothing)
    anchors = _circular_region_anchors(g; mark_between = mark_between,
                                       mark_between2 = mark_between2,
                                       bottom_leaves = bottom_leaves,
                                       top_leaves = top_leaves)
    isempty(anchors) && return ""
    cx = cy = size / 2
    # The dots sit on the SAME anchors as the region labels, so they would land
    # exactly on the text. The whole overlay is therefore shifted down by a fixed
    # amount — geometry unchanged (every endpoint moves alike), text clear.
    DY = 13.0
    XY(p) = (cx + r * p[1], cy + r * p[2] + DY)
    io = IOBuffer()

    # Arrows first, dots on top of them.
    for R in eachindex(parent)
        P = parent[R]
        P == 0 && continue
        (haskey(anchors, R) && haskey(anchors, P)) || continue
        lab = letters === nothing || R > length(letters) || letters[R] === nothing ?
              Int[] : letters[R]
        isempty(lab) && continue
        (x1, y1) = XY(anchors[P]); (x2, y2) = XY(anchors[R])
        # The arrow takes the colour of the LAST letter of the step; a
        # node crossing appends several and the label spells them all.
        ccol = lab[end]
        col = _crt_colour(ccol)
        # Stop the arrow short of the target dot so the head stays visible.
        dx = x2 - x1; dy = y2 - y1; L = sqrt(dx^2 + dy^2)
        L < 1e-6 && continue
        ex = x2 - 9 * dx / L; ey = y2 - 9 * dy / L
        print(io, """<line x1="$x1" y1="$y1" x2="$ex" y2="$ey" stroke="$col" """ *
                  """stroke-width="2" stroke-dasharray="5,3" """ *
                  """marker-end="url(#crtarrow$(ccol))"/>""")
        # The label: the simple reflection, at the midpoint, on a white plate.
        mx = (x1 + x2) / 2; my = (y1 + y2) / 2
        print(io, """<text x="$mx" y="$my" text-anchor="middle" dy="-0.4em" """ *
                  """font-size="12" font-weight="bold" fill="$col" stroke="white" """ *
                  """stroke-width="3" paint-order="stroke">$(join(lab))</text>""")
    end
    reachable(R) = words === nothing || (R <= length(words) && words[R] !== nothing)
    for (R, p) in anchors
        reachable(R) || continue
        (x, y) = XY(p)
        root = R <= length(parent) && parent[R] == 0
        # The root — the MARKED region, distance 0 — gets a double ring.
        print(io, """<circle cx="$x" cy="$y" r="$(root ? 6 : 4)" fill="#111" """ *
                  """stroke="white" stroke-width="1.5"/>""")
        root && print(io, """<circle cx="$x" cy="$y" r="9.5" fill="none" """ *
                          """stroke="#111" stroke-width="1.5"/>""")
    end

    # Arrow heads, one marker per colour (SVG needs them defined in the doc).
    defs = IOBuffer()
    print(defs, "<defs>")
    for c in 1:3
        col = _crt_colour(c)
        print(defs, """<marker id="crtarrow$c" viewBox="0 0 10 10" refX="9" refY="5" """ *
                    """markerWidth="6" markerHeight="6" orient="auto-start-reverse">""" *
                    """<path d="M 0 0 L 10 5 L 0 10 z" fill="$col"/></marker>""")
    end
    print(defs, "</defs>")
    return String(take!(defs)) * String(take!(io))
end
