# render/Decorated.jl — cell labels (polynomials) on the Tutte SVG.
#
# `decorated_svg(d)` draws the Tutte SVG of the graph and writes the polynomial
# label of every inner cell into it (anchors as in render/Cells.jl). A label of
# `1` is dropped so trivial labels don't clutter the picture.
#
# THE OUTER REGION has no cell anchor — its label sits OUTSIDE the boundary
# circle, top left. In debug mode (`debug_labels()`) every region normally
# carries a number; since the outer region has none, it is called "∞" (inner
# cells count from 1).

# (`_outer_label_svg` and `_OUTER_LABEL_DIR` live in render/Cells.jl; the other
#  caller is render/CircularCells.jl.)

"""
    decorated_svg(d::DecoratedDiagram; size = 600, r = 0.34·size,
                  mark_between = nothing) -> String

The Tutte SVG of `d.graph`, with the polynomial label in every inner cell
(canonical order from `inner_faces`) and the label of the OUTER region outside
the boundary circle (top left, with "∞" as the region name in debug mode).
With `mark_between = (a, b)`, the arc between leaves `a` and `b` is rotated to
the left and marked with a black boundary marker (as in `tutte_svg`).
"""
function decorated_svg(d::DecoratedDiagram; size::Int = 600, r::Real = 0.34 * size,
                       mark_between::MarkSpec = nothing,
                       slot_numbers::Bool = true,
                       distances::Union{Nothing,Vector{Int}} = nothing)
    g = d.graph
    # `slot_numbers` is only passed through, so slides can turn off the debug
    # arm numbers as everywhere else.
    base = tutte_svg(g; size = size, r = r, mark_between = mark_between,
                     slot_numbers = slot_numbers)
    # Debug mode shows `(number, polynomial)` — gated on the FLAG, not on
    # whether `distances` happens to be passed.
    dbg = distances !== nothing || debug_labels()
    txts = _region_label_texts(regions(g), d.face_labels, distances; debug = dbg)
    labels = _cell_labels_svg(g, txts; size = size, r = r, mark_between = mark_between)
    labels *= _outer_label_svg(d.outer_label; size = size, r = r, debug = dbg)
    return replace(base, "</svg>" => labels * "</svg>")
end

# Image display as for `WordGraph`/`MorphismGraph` (render/tutte/Display.jl): SVG
# for browser-Jupyter, PNG for VS Code (which ignores `image/svg+xml`). Without
# these methods `DecoratedDiagram` has NO `show`/`showable` for image MIME types
# and `display(d)` falls back to the text form — cell labels (alpha_i) would then
# be invisible in VS Code.
Base.show(io::IO, ::MIME"image/svg+xml", d::DecoratedDiagram) = print(io, decorated_svg(d))
Base.show(io::IO, ::MIME"image/png", d::DecoratedDiagram) = write(io, _svg_to_png(decorated_svg(d)))
Base.showable(::MIME"image/svg+xml", ::DecoratedDiagram) = true
Base.showable(::MIME"image/png", ::DecoratedDiagram) = true

"""
    display_decorated(d::DecoratedDiagram; mark_between = nothing)

Draws the decorated diagram (Tutte) with the polynomial of every inner cell.
With `mark_between = (a, b)`, also draws the left boundary marker.
"""
display_decorated(d::DecoratedDiagram; mark_between::MarkSpec = nothing) =
    (display("image/svg+xml", decorated_svg(d; mark_between = mark_between)); nothing)

"""
    display_combo(c::DiagramComboR; mark_between = nothing)

Draws every term of a `DiagramComboR`: coefficient line + SVG. If
`mark_between` is given, it is used for every term's marker.
"""
function display_combo(c::DiagramComboR; mark_between::MarkSpec = nothing)
    isempty(c) && (println("0 (empty sum)"); return nothing)
    for (d, coeff) in pairs_of(c)
        println(coeff, " · ", d)
        display_decorated(d; mark_between = mark_between)
    end
    return nothing
end

"""
    display(c::DiagramComboR)

A linear combination draws like its terms: coefficient line + image per term.

`display_combo` is the explicit form and the only one that accepts
`mark_between`; `display(c)` is the usual way.

WHY a `Base.display` method and not a `Base.show` MIME method: a combo is
SEVERAL images plus text lines, so it cannot be delivered as ONE
`image/svg+xml` document. `Base.display` is the level at which several
individual `display` calls are allowed.
"""
Base.display(c::DiagramComboR) = display_combo(c)
