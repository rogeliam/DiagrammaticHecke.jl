# render/CellDistance.jl — cell numbers, cell adjacency/BFS distance from the black
# mark, and their SVG renderers, plus the outer region's label. Included after
# Cells.jl and CellLabels.jl.

# ---- 1. cell numbers ----------------------------------------------------------

"""
    cell_number_svg(g; size = 600, r = 0.34·size) -> String

The Tutte SVG of `g` (for a MorphismGraph including the black mark), with the
REGION number in every visibly separate area.
"""
function cell_number_svg(g::WordGraph; size::Int = 600, r::Real = 0.34 * size)
    base = tutte_svg(g; size = size, r = r)
    labels = _cell_labels_svg(g, collect(1:region_count(g)); size = size, r = r)
    return replace(base, "</svg>" => labels * "</svg>")
end

function cell_number_svg(m::MorphismGraph; size::Int = 600, r::Real = 0.34 * size)
    base = tutte_svg(m.graph; size = size, r = r, mark_between = _left_mark(m))
    labels = _cell_labels_svg(m.graph, collect(1:region_count(m));
                              size = size, r = r, mark_between = _left_mark(m))
    return replace(base, "</svg>" => labels * "</svg>")
end

"""
    display_cell_number(g)

Draws the diagram (Tutte) with the number of every interior cell.
"""
display_cell_number(g) = (display("image/svg+xml", cell_number_svg(g)); nothing)

# ---- 2. distance from the black mark -----------------------------------------

"""
    cell_adjacency(g::WordGraph) -> Vector{Set{Int}}

Adjacency of the interior cells: two cells are adjacent when an edge separates
them (left/right in `port_cell`).
"""
function cell_adjacency(g::WordGraph)
    adj = [Set{Int}() for _ in 1:face_count(g)]
    for (_, (l, r)) in g.cells.port_cell
        (l >= 1 && r >= 1 && l != r) || continue
        push!(adj[l], r)
        push!(adj[r], l)
    end
    return adj
end

"""
    cell_distances(m::MorphismGraph) -> Vector{Int}

BFS distance of every interior cell from the black mark (the gap between bottom
and top start): touching cell = 0, its neighbours = 1, and so on. Unreachable cells
(free circles) get -1. Without a mark: all -1.
"""
function cell_distances(m::MorphismGraph)
    g = m.graph
    n = length(g.word)
    (n == 0 || _left_mark(m) === nothing) && return fill(-1, face_count(g))
    gap = mod1(m.cut1, n)                       # the mark sits on the cut1 arc
    start = gap_cell(g, gap)
    adj = cell_adjacency(g)
    dist = fill(-1, face_count(g))
    dist[start] = 0
    queue = [start]
    while !isempty(queue)
        c = popfirst!(queue)
        for nb in adj[c]
            dist[nb] == -1 || continue
            dist[nb] = dist[c] + 1
            push!(queue, nb)
        end
    end
    return dist
end

"""
    cell_distance_svg(m::MorphismGraph; size, r) -> String

Like `cell_number_svg`, but with the distance from the mark instead of the
number.

The distances still hang on the FACE notion (`cell_distances` runs over
`port_cell`) and are looked up per region via `R.cell`. Two regions of the same
face therefore show the same number in two places.
"""
function cell_distance_svg(m::MorphismGraph; size::Int = 600, r::Real = 0.34 * size)
    base = tutte_svg(m.graph; size = size, r = r, mark_between = _left_mark(m))
    dz = cell_distances(m)
    labels = [(d = dz[R.cell]; d < 0 ? "·" : string(d)) for R in regions(m)]
    labels = _cell_labels_svg(m.graph, labels;
                              size = size, r = r, mark_between = _left_mark(m))
    return replace(base, "</svg>" => labels * "</svg>")
end

"""
    display_cell_distance(m::MorphismGraph)

Draws the morphism (Tutte + black mark) with the BFS distance of every cell from
the mark.
"""
display_cell_distance(m::MorphismGraph) =
    (display("image/svg+xml", cell_distance_svg(m)); nothing)

# ---- the OUTER region's label -----------------------------------------------
# render/CircularCells.jl is the only caller.

# Direction the outer label sits in (unit circle = boundary circle; magnitude
# > 1 means "outside"). Top left, where diagrams reach least often.
const _OUTER_LABEL_DIR = (-0.82, -0.82)

"The SVG text snippet for the outer region's label (empty if there's nothing to show)."
# `f` is UNTYPED: `DecoratedDiagram{P}` is generic in the coefficient
# ring, and the body only needs `isone` and `show`.
function _outer_label_svg(f; size::Int, r::Real, debug::Bool)
    txt = isone(f) ? "" : string(f)
    # No "∞" is printed unless the label is non-trivial, so the outer region's marker
    # stays out of debug pictures with a trivial label. `debug` stays as a parameter
    # (callers pass it through).
    isempty(txt) && return ""
    cx = cy = size / 2
    X = cx + r * _OUTER_LABEL_DIR[1]
    Y = cy + r * _OUTER_LABEL_DIR[2]
    return """<text x="$X" y="$Y" text-anchor="middle" dy="0.35em" """ *
           """font-size="17" font-style="italic" fill="#555" stroke="white" """ *
           """stroke-width="4" paint-order="stroke">$(_polynomial_label_svg(txt))</text>"""
end
