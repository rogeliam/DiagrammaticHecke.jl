# render/CircularCells.jl — cell labels in the circular-Tutte SVG.
#
# Cell anchors and labels for `CircularGraph`/`CircularDecorated`, via
# `_circular_face_walks` (CircularFaces.jl), `circular_tutte_positions`
# (CircularTutte.jl) and the anchor geometry of render/Cells.jl.

# Anchor = centroid of the cell polygon.
function _circular_cell_anchors(g::CircularGraph; mark_between::MarkSpec = nothing,
                           mark_between2::MarkSpec = nothing,
                           bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                           top_leaves::Union{Nothing,Vector{Int}} = nothing)
    walks, cell_of_face, _, _ = _circular_face_walks(g)
    leafxy, nodexy = circular_tutte_positions(g; mark_between = mark_between,
                                         mark_between2 = mark_between2,
                                         bottom_leaves = bottom_leaves, top_leaves = top_leaves)
    return _cell_anchors_from(walks, cell_of_face, inner_faces(g),
                              leafxy, nodexy, length(g.word))
end

"""
    _circular_region_anchors(g::CircularGraph; …) -> Dict{Int,Tuple}

Anchors per REGION for a `CircularGraph` (the renderer labels regions, not faces —
see [`_region_anchors`](@ref)). Passes
`bottom_leaves`/`top_leaves` through, like `_circular_cell_anchors`, and likewise
`mark_between2`: the grey marker enters via `_add_nudge!` in
`circular_tutte_positions`, so the anchors must come from the SAME positions as the
image, or the labels sit shifted.
"""
function _circular_region_anchors(g::CircularGraph; mark_between::MarkSpec = nothing,
                             mark_between2::MarkSpec = nothing,
                             bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                             top_leaves::Union{Nothing,Vector{Int}} = nothing)
    walks, cell_of_face, _, _ = _circular_face_walks(g)
    leafxy, nodexy = circular_tutte_positions(g; mark_between = mark_between,
                                         mark_between2 = mark_between2,
                                         bottom_leaves = bottom_leaves,
                                         top_leaves = top_leaves)
    return _region_anchors(regions(g), walks, cell_of_face, leafxy, nodexy,
                           length(g.word))
end

function _circular_cell_labels_svg(g::CircularGraph, labels; size::Int, r::Real,
                              mark_between::MarkSpec = nothing,
                              mark_between2::MarkSpec = nothing,
                              bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                              top_leaves::Union{Nothing,Vector{Int}} = nothing)
    anchors = _circular_region_anchors(g; mark_between = mark_between,
                                  mark_between2 = mark_between2,
                                  bottom_leaves = bottom_leaves, top_leaves = top_leaves)
    cx = cy = size / 2
    io = IOBuffer()
    for c in eachindex(labels)
        haskey(anchors, c) || continue
        txt = string(labels[c])
        isempty(txt) && continue
        (x, y) = anchors[c]
        X = cx + r * x; Y = cy + r * y
        svg_txt = _polynomial_label_svg(txt)
        # Font size 13, not 17 (make the font a
        # bit smaller, otherwise it can't be read in the images) — with region
        # WORDS in the cell (`cell_labels = :words`) size 17 runs over the
        # cell edges and becomes unreadable at the usual on-screen width.
        print(io, """<text x="$X" y="$Y" text-anchor="middle" dy="0.35em" """ *
                  """font-size="13" font-style="italic" fill="#555" stroke="white" """ *
                  """stroke-width="3" paint-order="stroke">$svg_txt</text>""")
    end
    return String(take!(io))
end

"""
    circular_decorated_svg(fd::CircularDecorated; size = 600, r = 0.34·size,
                      mark_between = nothing, region_based = false) -> String

The circular-Tutte SVG of `fd.graph`, with the polynomial label in each **region**
(order of `regions`) and the label of the OUTER region outside the boundary
circle (top left, shown as "∞" in debug mode). Analogue of `decorated_svg`
(render/Decorated.jl), but via `circular_tutte_svg`/`_circular_cell_anchors`.

`distances` can add distances. With `region_based=true` (the current
convention), `distances` refers to regions (length `region_count`); otherwise
to faces (length `face_count`) — this affects ONLY the distances, labels are
always region-indexed.

`mark_between2` is the grey boundary marker. It MUST be passed through for the
image to match `circular_tutte_svg(m)`: it does enter `circular_tutte_positions`
(`_add_nudge!(mark_between2)`, render/CircularTutte.jl) and shifts two leaves by
0.12 on the unit disc.
"""
function circular_decorated_svg(fd::CircularDecorated; size::Int = 600, r::Real = 0.34 * size,
                          mark_between::MarkSpec = nothing,
                          mark_between2::MarkSpec = nothing,
                          bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                          top_leaves::Union{Nothing,Vector{Int}} = nothing,
                          distances::Union{Nothing,Vector{Int}} = nothing,
                          # LEAF NUMBERS ALWAYS ON — same house convention as
                          # render/CircularTutte.jl; ARM NUMBERS hang off the debug
                          # switch there too.
                          leaf_numbers::Bool = true,
                          arm_numbers::Bool = _circular_debug_on(),
                          region_based::Bool = true)
    g = fd.graph
    # Image and labelling MUST use the same orientation — see the comment at
    # the `cell_numbers` spot in CircularTutte.jl.
    base = circular_tutte_svg(g; size = size, r = r, mark_between = mark_between,
                         mark_between2 = mark_between2,
                         bottom_leaves = bottom_leaves, top_leaves = top_leaves,
                         leaf_numbers = leaf_numbers, arm_numbers = arm_numbers)
    # Debug mode shows `(number, polynomial, distance)` — gated by the FLAG,
    # not by whether `distances` happens to be given.
    dbg = debug_circular() || debug_labels()
    txts = _region_label_texts(regions(g), fd.region_labels, distances;
                               debug = dbg,
                               region_based = region_based)
    labels = _circular_cell_labels_svg(g, txts;
                                  size = size, r = r, mark_between = mark_between,
                                  mark_between2 = mark_between2,
                                  bottom_leaves = bottom_leaves, top_leaves = top_leaves)
    # THE OUTER REGION has no region anchor; its label sits OUTSIDE the
    # boundary circle (top left). Literally the same building block as the
    # plain path — `_outer_label_svg` (render/Decorated.jl), shown as "∞" in
    # debug mode.
    labels *= _outer_label_svg(fd.outer_label; size = size, r = r, debug = dbg)
    return replace(base, "</svg>" => labels * "</svg>")
end

# The boundary marking always belongs in the image. `display(g::CircularGraph)`
# sets the CANONICAL mark `_circular_default_mark` (CircularTutte.jl, black between leaf
# n and 1); the decorated path here does the same, so `display(fd)` never
# draws without a mark, which would leave the ring orientation undetermined
# and the leaf ring possibly reversed. Both paths behave identically.
Base.show(io::IO, ::MIME"image/svg+xml", d::CircularDecorated) =
    print(io, circular_decorated_svg(d; mark_between = _circular_default_mark(d.graph),
                                leaf_numbers = true))
Base.show(io::IO, ::MIME"image/png", d::CircularDecorated) =
    write(io, _svg_to_png(circular_decorated_svg(d; mark_between = _circular_default_mark(d.graph),
                                            leaf_numbers = true)))
Base.showable(::MIME"image/svg+xml", ::CircularDecorated) = true
Base.showable(::MIME"image/png", ::CircularDecorated) = true

"""
    display_circular_decorated(fd::CircularDecorated; mark_between = nothing,
                          bottom_leaves = nothing, top_leaves = nothing,
                          distances = nothing, region_based = false)

Draws the decorated CircularGraph (Tutte) with the polynomial of each inner cell.

**Prefer `display(fd)`.** `CircularDecorated` has the `Base.show` MIME methods
(directly above), so a call without keywords is literally `display(fd)`; this
function is a thin pass-through to `circular_decorated_svg(::CircularDecorated)`
that also takes the keywords. `bottom_leaves`/`top_leaves` decide the
MIRRORING of the leaf ring and must be passed along whenever `fd` originates
from a morphism.

WATCH OUT for a common confusion: `display(fm)` on a `CircularMorphismGraph` draws
NO polynomial labels — only this decorated path does. The cell NUMBERS in
debug mode (`debug_circular()`) are a different thing from the polynomial LABELS.
"""
function display_circular_decorated(fd::CircularDecorated; mark_between::MarkSpec = nothing,
                               bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                               top_leaves::Union{Nothing,Vector{Int}} = nothing,
                               distances::Union{Nothing,Vector{Int}} = nothing,
                               region_based::Bool = false)
    if mark_between === nothing && bottom_leaves === nothing &&
       top_leaves === nothing && distances === nothing
        display(fd)          # uses the registered Base.show MIME method
    else
        display("image/svg+xml",
                circular_decorated_svg(fd; mark_between = mark_between,
                                  bottom_leaves = bottom_leaves, top_leaves = top_leaves,
                                  distances = distances, region_based = region_based))
    end
    return nothing
end

"""
    show_circular_combo(c::CircularComboR; mark_between = nothing,
                   bottom_leaves = nothing, top_leaves = nothing, cuts = nothing)

Draws each term of a `CircularComboR`: coefficient line + SVG (a sum should never be
left as a bare count). Analogue of `display_combo`
(render/Decorated.jl).

`bottom_leaves`/`top_leaves` fix the ORIENTATION of the leaf ring (they
determine the leaf-angle assignment, `_halfcircle_angles`); without them the
decorated path draws with a different orientation than `circular_tutte_svg(m)`.

`cuts = (cut1, cut2)` are the cuts of the surrounding morphism. A
`CircularDecorated` has none — they are needed to compute each term's region
distances, which in turn pick the one region a face's polynomial is shown in.
"""
function show_circular_combo(c::CircularComboR; mark_between::MarkSpec = nothing,
                        bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                        top_leaves::Union{Nothing,Vector{Int}} = nothing,
                        cuts::Union{Nothing,Tuple{Int,Int}} = nothing)
    isempty(c) && (println("0 (empty sum)"); return nothing)
    for (d, coeff) in pairs_of(c)
        println(coeff, " · ", d)
        if !(d isa CircularDecorated)
            display(d)
        elseif cuts !== nothing
            # THE CORRECT WAY: draw the term as a CircularDecoratedMorphism.
            # `circular_decorated_svg(::CircularDecoratedMorphism)`
            # (circular/CircularDecoratedRules.jl) gets marker, orientation AND distances
            # from the morphism itself.
            #
            # The region structure hangs off the TERM, not the original
            # morphism: D4 splits and merges regions. The CircularMorphismGraph is
            # therefore rebuilt from `cuts` for each term.
            display(CircularDecoratedMorphism(CircularMorphismGraph(d.graph, cuts[1], cuts[2]),
                                         d.region_labels, d.outer_label))
        else
            # Without cuts a `CircularDecorated` has NO orientation (no bottom/top,
            # no marking) — nodes can then appear rotated by up to 180°. Whoever
            # has the orientation should pass it (2-arg form below); without it,
            # a WARNING is raised here instead of a silently wrong picture.
            @warn("show_circular_combo without `cuts`: the term is drawn without " *
                  "orientation (nodes can appear rotated by up to 180°). " *
                  "Use `show_circular_combo(c, m)` with the surrounding CircularMorphismGraph.",
                  maxlog = 1)
            display_circular_decorated(d; mark_between = mark_between,
                                  bottom_leaves = bottom_leaves,
                                  top_leaves = top_leaves)
        end
    end
    return nothing
end

"""
    show_circular_combo(c::CircularComboR, m::CircularMorphismGraph)

Convenience overload: takes everything the decorated path needs directly from
`m` instead of requiring it at the call site — the usual situation in
examples, where `c` comes from `reduce_to_circular_leave(CircularDecoratedMorphism(m,
...))` and should carry the same marking and orientation as `m` (the reduction
moves neither marking nor cuts):

- `mark_between = _circular_left_mark(m)` — the black boundary marker;
- `bottom_leaves`/`top_leaves` — the ORIENTATION, as `_circular_tutte_morphism`
  gives it to `circular_tutte_svg`;
- `cuts = (m.cut1, m.cut2)` — from these, each term builds its own
  `CircularDecoratedMorphism`, which determines marker, orientation AND region
  distances itself.

The grey marker (`_circular_right_mark`) is passed through this path as well; see
`test/circular.jl` for the regression test.
"""
show_circular_combo(c::CircularComboR, m::CircularMorphismGraph) =
    show_circular_combo(c; mark_between = _circular_left_mark(m),
                   bottom_leaves = _bottom_leaves(m), top_leaves = _top_leaves(m),
                   cuts = (m.cut1, m.cut2))

"""
    display(c::CircularComboR)

Like `display(::DiagramComboR)` (render/Decorated.jl): coefficient line +
image per term. `show_circular_combo` remains as a name (it takes `mark_between`),
`display(c)` is the usual way.
"""
Base.display(c::CircularComboR) = show_circular_combo(c)
