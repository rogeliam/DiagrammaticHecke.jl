# render/CircularTutteSVG.jl — the SVG emitter for the circular-Tutte drawing.
#
# Included directly after CircularTutte.jl, so the layout
# (`circular_tutte_positions`, `_circular_leg_angles`, `_circular_node_radius`, the
# debug switch) is already defined.
#
# What lives here: `circular_tutte_svg`, `Base.show`/`Base.showable` for both MIME
# types, `display_circular_tutte` and `show_circular_step`.

# ---- 7. SVG --------------------------------------------------------------------

"""
    circular_tutte_svg(g::CircularGraph; size = 600, r = 0.34·size, iters = 400, repel = 0.006,
                 mark_between = nothing, mark_between2 = nothing,
                 bottom_leaves = nothing, top_leaves = nothing,
                 slot_numbers = true, edge_numbers = slot_numbers,
                 cell_numbers = false) -> String

Draws `g` as a planar Tutte embedding: boundary leaves as labelled discs on the
dashed outer circle, internal circular nodes at their barycentres, legs as cubic
BÉZIER curves from node boundary to node boundary (control points along the
slot angles from `_circular_leg_angles`, length `k = 0.4·|B-A|` of the secant). Free
circles and both boundary markers work as in `tutte_svg`
(render/tutte/Tutte.jl). A completely empty graph draws only the dashed circle
(early return, as in render/tutte/Tutte.jl's `_empty_diagram_svg`).

`arm_numbers` (default `false`) writes the arm number (slot 1..k) on each leg.

`slot_numbers` (default `true`) is the umbrella switch that `node_numbers` and
`edge_numbers` hang off. `node_numbers` (default = `slot_numbers`) writes the
NODE number centred INSIDE the node disc (white on the node colour; the disc
is enlarged if needed) — the same number as in `show_wiring`. For the full
WIRING (arm → target), use `wiring_table(g)` (render/tutte/Wiring.jl) as a
table next to the image rather than labels on the diagram itself.

`edge_numbers` (default = `slot_numbers && !arm_numbers`) writes the EDGE
number (index in `g.edges`) once per edge, directly on the node boundary in
slot direction at the edge's first node end — plain black, no badge. This is
the number rules and measurements use to name an edge
(`circular_2parallel_apply(g, ei, ej)`, `check_wiring`, `circular_region_adjacency`).

`cell_numbers` (default `false`) writes the REGION number in each region;
together with `slot_numbers` this is the debug overlay `display_circular_tutte`
draws under `debug_circular(true)`.

`leaf_numbers` (default `false`) writes the LEAF number `k` (the index in
`g.word`, as it also appears in `Leaf(k)` in an edge dump) just OUTSIDE the
dashed boundary circle next to each boundary leaf — lets you see which
`Leaf(k)` a given node arm connects to without first reading the edge list.

`cell_labels` (default `:numbers`) selects WHAT is shown in a region:
`:numbers` (the region number), `:distances` (the BFS distance to the black
marking), `:words` (the REGION WORD, `ε` at the mark), `:words_num` (`R|word`),
or `:none`. Only has an effect when `cell_numbers = true`.

`region_tree` (default `false`) overlays the DUAL GRAPH: the BFS tree of
`circular_region_words` from the marked region to every other one, one arrow per
tree edge, labelled with the simple reflection `s` crossed there — so the arrows
from the mark to a region spell out that region's region word. Like
`:distances`/`:words` it needs the marking, i.e. draw the MORPHISM, not
`m.graph`. See `render/CircularRegionTree.jl`.

**`:distances` requires the MARKING and thus a `CircularMorphismGraph`** —
`circular_region_distances` (circular/CircularRegion.jl) measures from the gap between the
bottom and top start. A bare `CircularGraph` has none, so `:distances` gives `·`
(unreachable) everywhere there. To see distances, draw the morphism, not
`m.graph`.
`edge_repel` (default `true`) makes edges AVOID foreign nodes instead
of running straight over them (render/tutte/Avoid.jl). At legs of a `:braid`
node this happens only within tight angle limits, so the slot order is
preserved. `edge_repel = false` draws every edge straight.

`dot_gap` (default `12` px = one dot diameter) is the clearance between the
edge of a dot and the edge of its neighbour (node or leaf disc).
"""
function circular_tutte_svg(g::CircularGraph; size::Int = 600, r::Real = 0.34 * size,
                       iters::Int = 400, repel::Real = 0.02,
                       clip::Real = 0.92,
                       mark_between::MarkSpec = nothing,
                       mark_between2::MarkSpec = nothing,
                       bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                       top_leaves::Union{Nothing,Vector{Int}} = nothing,
                       slot_numbers::Bool = true,
                       node_numbers::Bool = slot_numbers,
                       arm_numbers::Bool = false,
                       # EDGE NUMBERS YIELD TO ARM NUMBERS: both sit on the node
                       # boundary; when arm numbers are on (debug mode), edge
                       # numbers are dropped to avoid clutter. An explicitly
                       # passed `edge_numbers = true` still applies.
                       edge_numbers::Bool = slot_numbers && !arm_numbers,
                       cell_numbers::Bool = false,
                       cell_labels::Symbol = :numbers,
                       leaf_numbers::Bool = false,
                       distances::Union{Nothing,Vector{Int}} = nothing,
                       region_words::Union{Nothing,Vector} = nothing,
                       region_tree::Bool = false,
                       region_parent::Union{Nothing,Vector{Int}} = nothing,
                       region_letters::Union{Nothing,Vector} = nothing,
                       region_labels = nothing,
                       edge_repel::Bool = true, dot_gap::Real = 2 * _DOT_R)
    # A free monochrome circle (`Edge(c, Circle(c), Circle(c))`, no leaves/nodes
    # otherwise) must NOT take the empty-diagram fast path: that would skip
    # drawing the free-circle ring itself (only the dashed boundary would show).
    # Analogue of the same fix in render/tutte/Tutte.jl.
    isempty(g.word) && isempty(g.nodes) && isempty(g.edges) &&
        mark_between === nothing && mark_between2 === nothing &&
        return _empty_diagram_svg(size)

    leafxy, nodexy = circular_tutte_positions(g; iters = iters, repel = repel,
                                         clip = clip,
                                        mark_between = mark_between, mark_between2 = mark_between2,
                                        bottom_leaves = bottom_leaves, top_leaves = top_leaves,
                                        dot_gap = dot_gap)
    n = length(g.word)
    baseθ(k) = pi/2 + 2pi * (k - 1) / max(n, 1)
    rot, reflect = _ring_transform(baseθ, n, mark_between, bottom_leaves, top_leaves)
    cx = cy = size / 2

    hc_sv = _halfcircle_angles(n, bottom_leaves, top_leaves)
    # same case distinction as in `circular_tutte_positions`
    structured_sv = !(mark_between !== nothing &&
                      !(mark_between isa Int && mark_between == 0) &&
                      (bottom_leaves === nothing || top_leaves === nothing ||
                       isempty(bottom_leaves) || isempty(top_leaves)))
    function _marker_pt(mb::MarkSpec)
        # same as the plain renderer: under half-circle assignment, derive from
        # the flank angles
        if structured_sv && !(mb isa Symbol)
            θh = _halfcircle_marker_angle(hc_sv, n, mb)
            θh !== nothing && return (cx + r * cos(θh), cy + r * sin(θh))
        end
        θraw = _mark_angle(baseθ, n, mb)
        θraw === nothing && return nothing
        θ = mb isa Symbol ? θraw : reflect * (θraw + rot)
        return (cx + r * cos(θ), cy + r * sin(θ))
    end
    X(p) = cx + r * p[1];  Y(p) = cy + r * p[2]
    pos(port::Leaf) = (X(leafxy[port.k]), Y(leafxy[port.k]))
    pos(port::NodePort) = (X(nodexy[port.node]), Y(nodexy[port.node]))

    θangles = _circular_leg_angles(g, leafxy, nodexy)      # in unit-disc coordinates
    # Leg direction of a slot in SVG pixel coordinates (the scaling is
    # isotropic: X/Y both scale with the same `r`, so the angle is unchanged
    # by the unit-disc → pixel transition).
    legdir(v, s) = (cos(θangles[v][s]), sin(θangles[v][s]))

    io = IOBuffer()
    print(io, """<svg width="$size" height="$size" viewBox="0 0 $size $size" xmlns="http://www.w3.org/2000/svg" font-family="DejaVu Sans Mono, Consolas, Liberation Mono, Courier New, monospace">""")
    print(io, """<circle cx="$cx" cy="$cy" r="$r" fill="none" stroke="#e2e7ee" stroke-dasharray="3" stroke-width="1"/>""")

    # ---- Legs: one cubic Bézier per edge with at least one node end. No
    # bundling/fan needed (see module docstring): different slots get different
    # departure angles automatically via the monotonicity repair.
    drawable = [e for e in g.edges if !(e.a isa Circle || e.b isa Circle)]
    for e in drawable
        (x1, y1) = pos(e.a); (x2, y2) = pos(e.b)
        col = get(_RGB, e.colour, "#666")
        # Start point/direction: if `a` is a node slot, start exactly at the node
        # boundary along its slot angle; a leaf end has no "boundary" of its own —
        # there the anchor stays the leaf position itself (as in the plain renderer).
        if e.a isa NodePort
            nd = g.nodes[e.a.node]
            rad = _circular_node_radius(nd)
            (dx, dy) = legdir(e.a.node, e.a.slot)
            (x1, y1) = (x1 + rad * dx, y1 + rad * dy)
        else
            dx = (x2 - x1); dy = (y2 - y1)
            L = hypot(dx, dy); (dx, dy) = L > 1e-9 ? (dx / L, dy / L) : (0.0, 0.0)
        end
        if e.b isa NodePort
            nd = g.nodes[e.b.node]
            rad = _circular_node_radius(nd)
            (dx2, dy2) = legdir(e.b.node, e.b.slot)
            (x2, y2) = (x2 + rad * dx2, y2 + rad * dy2)
        else
            dx2 = (x1 - x2); dy2 = (y1 - y2)
            L2 = hypot(dx2, dy2); (dx2, dy2) = L2 > 1e-9 ? (dx2 / L2, dy2 / L2) : (0.0, 0.0)
        end
        secL = hypot(x2 - x1, y2 - y1)
        k = 0.4 * secL
        # If this edge connects to a leaf OR a dot-node (arm_count == 1) at one
        # end, draw it straight (no cubic curvature) so dots/leaves appear to
        # go straight out of the node. A dot is NOT a Leaf but a CircularGraph node
        # with arm_count == 1 — without this special case it would fall into the
        # cubic branch below and get visibly bulged by k = 0.4*secL despite the
        # short dot distance.
        a_is_stub = e.a isa Leaf || (e.a isa NodePort && arm_count(g.nodes[e.a.node]) == 1)
        b_is_stub = e.b isa Leaf || (e.b isa NodePort && arm_count(g.nodes[e.b.node]) == 1)
        # EDGE REPULSION (render/tutte/Avoid.jl):
        # foreign nodes (dots included) are obstacles, the curve avoids them.
        obstacles = _circular_obstacles(g, e, nodexy, X, Y, edge_repel)
        if a_is_stub || b_is_stub
            (qx, qy, moved) = _avoid_quad(x1, y1, (x1+x2)/2, (y1+y2)/2, x2, y2, obstacles)
            if moved
                print(io, """<path d="M $x1 $y1 Q $qx $qy $x2 $y2" stroke="$col" stroke-width="3.4" fill="none"/>""")
            else
                print(io, """<path d=\"M $x1 $y1 L $x2 $y2\" stroke=\"$col\" stroke-width=\"3.4\" fill=\"none\"/>""")
            end
        else
            c1x = x1 + k * dx; c1y = y1 + k * dy
            c2x = x2 + k * dx2; c2y = y2 + k * dy2
            (c1x, c1y, c2x, c2y, _) = _avoid_cubic(x1, y1, c1x, c1y, c2x, c2y, x2, y2, obstacles)
            print(io, """<path d="M $x1 $y1 C $c1x $c1y $c2x $c2y $x2 $y2" stroke="$col" stroke-width="3.4" fill="none"/>""")
        end
    end

    # free monochrome circles (Edge(c, Circle(c), Circle(c))) — as in render/tutte/Tutte.jl.
    freecircles = [e.colour for e in g.edges if e.a isa Circle && e.b isa Circle]
    for (j, c) in enumerate(freecircles)
        col = get(_RGB, c, "#666")
        fx = cx + (j - (length(freecircles) + 1) / 2) * 34.0
        fy = cy
        print(io, """<circle cx="$fx" cy="$fy" r="14" fill="none" stroke="$col" stroke-width="3.4"/>""")
    end

    # internal circular nodes: degree 1 ⇒ SOLID dot (as in render/tutte/Tutte.jl),
    # otherwise a filled disc with degree-dependent radius.
    for (v, nd) in enumerate(g.nodes)
        (x, y) = (X(nodexy[v]), Y(nodexy[v]))
        rad = _circular_node_radius(nd)
        col = _circular_node_fill(nd)
        if arm_count(nd) == 1
            # degree 1 = dot generator: SOLID (see render/tutte/Tutte.jl) —
            # hollow reads as a hole when enlarged.
            print(io, """<circle cx="$x" cy="$y" r="$rad" fill="$col" stroke="$col" stroke-width="3" class="dot"/>""")
        else
            print(io, """<circle cx="$x" cy="$y" r="$rad" fill="$col"/>""")
        end
    end

    # Draw a small loop for nodes with a self-connection (R10 case): when two
    # slots of a node are directly connected to each other (an edge between two
    # NodePorts of the same node), draw a loop at the node boundary. Not only
    # trivalent nodes (arm_count == 3) can have this — e.g. `merge_at_edge` can
    # produce a 4-armed node with two consecutive connected slots.
    # `_circular_self_loop_slots` returns `nothing` when no self-connection exists,
    # so the check below applies to any degree.
    for (v, nd) in enumerate(g.nodes)
        # compute this node's centre and radius fresh each iteration (not
        # reused from the previous loop iteration)
        (x, y) = (X(nodexy[v]), Y(nodexy[v]))
        rad = _circular_node_radius(nd)
        if arm_count(nd) >= 2
            s_pair = _circular_self_loop_slots(g, v)
            if s_pair !== nothing
                s, t = s_pair
                a1 = θangles[v][s]
                a2 = θangles[v][t]
                # Start/end point at the node boundary — same construction as the
                # normal legs above, just with both ends at the same node `v`. A
                # cubic Bézier rather than a free-standing circle, so the loop
                # actually starts AT the node and curves back to it (a free
                # `<circle>` would float behind the node boundary, ignoring
                # a1/a2 entirely).
                (dx1, dy1) = (cos(a1), sin(a1))
                (dx2, dy2) = (cos(a2), sin(a2))
                (lx1, ly1) = (x + rad * dx1, y + rad * dy1)
                (lx2, ly2) = (x + rad * dx2, y + rad * dy2)
                # Control points along the respective slot direction, pointing
                # outward; the arm length is large enough for a clearly visible
                # loop (~2.5-3·rad, here 2.8·rad).
                bulge = 2.8 * rad
                c1x = lx1 + bulge * dx1; c1y = ly1 + bulge * dy1
                c2x = lx2 + bulge * dx2; c2y = ly2 + bulge * dy2
                # use the self-connecting edge's own colour (not a hardcoded
                # grey) so the loop reads as part of the same coloured strand
                # as the node's other legs.
                loop_colour = only(e.colour for e in g.edges
                                    if e.a isa NodePort && e.b isa NodePort &&
                                       e.a.node == v && e.b.node == v &&
                                       ((e.a.slot == s && e.b.slot == t) ||
                                        (e.a.slot == t && e.b.slot == s)))
                loop_col = get(_RGB, loop_colour, "#666")
                print(io, """<path d="M $lx1 $ly1 C $c1x $c1y $c2x $c2y $lx2 $ly2" stroke="$loop_col" stroke-width="3.4" fill="none"/>""")
            end
        end
    end

    # Slot numbers (debug aid): per node, write the running arm number
    # 1..arm_count just outside the node boundary on the respective leg. This
    # lets you read off from the DRAWING which slot is wired where — colour
    # sequences alone don't reveal that (two nodes can have the same arm
    # sequence but different wiring). Toggled via `arm_numbers`; the TARGET of
    # each leg is not shown here, only in `wiring_table`
    # (render/tutte/Wiring.jl).
    # ⚠ Default is OFF (`arm_numbers`, not `slot_numbers`): showing both the arm
    # number and the edge number on every leg reads as double-labelling. The
    # edge number (below) is what rules and measurements use to name a strand.
    # Switch `arm_numbers = true` for the slot numbering, or read `wiring_table(g)`.
    if arm_numbers
        for (v, nd) in enumerate(g.nodes)
            arm_count(nd) == 1 && continue          # dots have only one arm
            (x, y) = (X(nodexy[v]), Y(nodexy[v]))
            rad = _circular_node_radius(nd)
            for s in 1:arm_count(nd)
                (dx, dy) = (cos(θangles[v][s]), sin(θangles[v][s]))
                # Offset: at `rad + 24` the arm numbers sat so close to the edge
                # numbers (which sit on the node boundary at `rad`) that the two
                # digit rows ran into each other and looked like one number.
                # `rad + 42` visibly separates the two rows.
                lx = x + (rad + 42) * dx
                ly = y + (rad + 42) * dy
                print(io, """<circle cx="$lx" cy="$ly" r="6" fill="#ffffff" fill-opacity="0.85" stroke="none"/>""")
                print(io, """<text x="$lx" y="$(ly + 3)" text-anchor="middle" font-size="9" fill="#8a8f9a">$s</text>""")
            end
        end
    end

    # EDGE NUMBERS ON THE NODE. The index of the edge in `g.edges` is the number
    # used by every rule and measurement to name an edge
    # (`circular_2parallel_apply(g, ei, ej)`, `check_wiring`, the adjacency lists).
    #
    # Drawn once PER EDGE, not per port (so it does not double up with the arm
    # number as two numbers per leg), directly ON the node boundary in slot
    # direction (radius `rad`), plain black, no badge — the node itself carries
    # the colour.
    if edge_numbers
        for (i, e) in enumerate(g.edges)
            port = e.a isa NodePort ? e.a : (e.b isa NodePort ? e.b : nothing)
            port === nothing && continue     # leaf-leaf edge: no node, no place to draw
            v = port.node
            (x, y) = (X(nodexy[v]), Y(nodexy[v]))
            rad = _circular_node_radius(g.nodes[v])
            (dx, dy) = (cos(θangles[v][port.slot]), sin(θangles[v][port.slot]))
            lx = x + (rad + 11) * dx
            ly = y + (rad + 11) * dy
            print(io, """<text x="$lx" y="$(ly + 3.5)" text-anchor="middle" font-size="10" """ *
                      """font-weight="bold" fill="#000000">$i</text>""")
        end
    end

    # NODE NUMBERS (debug aid). The index `v` from `g.nodes` is drawn centred in
    # the node — the same number `show_wiring` uses as "node v" in the wiring
    # table and that every rule uses to name its positions (`(tv, bi)`). Tied to
    # `slot_numbers` by default (`node_numbers = slot_numbers`), the same debug
    # switch as the arm numbers; individually toggleable.
    if node_numbers
        for (v, nd) in enumerate(g.nodes)
            (x, y) = (X(nodexy[v]), Y(nodexy[v]))
            rad = _circular_node_radius(nd)
            # The number sits INSIDE the node. To fit a one- or two-digit
            # number, the disc is enlarged to at least `9 + 2·digits` px — the
            # node colour stays, the digit is white on top.
            digits = length(string(v))
            r2 = max(rad, 8.0 + 2.5 * digits)
            r2 > rad && print(io, """<circle cx="$x" cy="$y" r="$r2" fill="$(_circular_node_fill(nd))" stroke="none"/>""")
            fs = digits >= 2 ? 11 : 12
            print(io, """<text x="$x" y="$(y + 0.36 * fs)" text-anchor="middle" """ *
                      """font-size="$fs" font-weight="bold" fill="#ffffff">$v</text>""")
        end
    end

    # Boundary leaves (labelled discs, same style as render/tutte/Tutte.jl) —
    # including the red "✗" marking for a LOOSE leaf (`_circular_leaf_colour` returns
    # 0 there, and colour 0 doesn't exist). `check_wiring` reports the same
    # condition as `:leaf_degree`.
    for k in 1:n
        c = _circular_leaf_colour(g, k)
        (x, y) = (X(leafxy[k]), Y(leafxy[k]))
        if c == 0
            print(io, """<circle cx="$x" cy="$y" r="10" fill="#ffecec" stroke="#c62828" stroke-width="3" stroke-dasharray="3 2"/>""")
            print(io, """<text x="$x" y="$y" text-anchor="middle" dy="0.35em" font-size="12" fill="#c62828">✗</text>""")
        else
            col = get(_RGB, c, "#666")
            print(io, """<circle cx="$x" cy="$y" r="10" fill="white" stroke="$col" stroke-width="2.8"/>""")
            print(io, """<text x="$x" y="$y" text-anchor="middle" dy="0.35em" font-size="12" fill="$col">$c</text>""")
        end
    end

    # LEAF NUMBERS (debug aid): the index k from the word, just OUTSIDE the
    # boundary circle, radially outward from the centre (analogous to the arm
    # numbers above, but at the leaf instead of the node), outside the leaf
    # disc, which already shows the colour, not the index.
    if leaf_numbers
        for k in 1:n
            (x, y) = (X(leafxy[k]), Y(leafxy[k]))
            (dx, dy) = (x - cx, y - cy)
            L = hypot(dx, dy); (dx, dy) = L > 1e-9 ? (dx / L, dy / L) : (0.0, -1.0)
            lx = x + 15 * dx; ly = y + 15 * dy
            print(io, """<text x="$lx" y="$(ly + 3.5)" text-anchor="middle" font-size="10" fill="#334">$k</text>""")
        end
    end

    (pt = _marker_pt(mark_between)) !== nothing &&
        print(io, """<circle cx="$(pt[1])" cy="$(pt[2])" r="6.5" fill="#000"/>""")
    (pt2 = _marker_pt(mark_between2)) !== nothing &&
        print(io, """<circle cx="$(pt2[1])" cy="$(pt2[2])" r="6.5" fill="#999"/>""")

    # Cell numbers (debug aid): each REGION gets its number. Regions (solid
    # boundary), not faces — a face touching the boundary at several separate
    # places is labelled at each.
    #
    # `bottom_leaves`/`top_leaves` MUST be passed through: they determine the
    # leaf-angle assignment (`_halfcircle_angles`). Without them the anchors are
    # computed in a DIFFERENT layout than the one actually drawn, and the
    # numbers land in the wrong cell.
    if cell_numbers && cell_labels !== :none
        cell_labels in (:numbers, :distances, :labels, :labels_dist, :words, :words_num) ||
            throw(ArgumentError("cell_labels must be :numbers, :distances, :labels, " *
                                ":labels_dist, :words, :words_num or :none, " *
                                "not :$cell_labels"))
        # `:distances` without given values means: there is no marking (bare
        # CircularGraph) — then `·` is shown everywhere, as for an unreachable region
        # in the decorated path.
        txts = if cell_labels === :numbers
            collect(1:region_count(g))
        elseif cell_labels === :words || cell_labels === :words_num
            # The REGION WORD in the cell. Like
            # `:distances` this needs the MARKING, so `region_words` is filled by
            # `_circular_tutte_morphism`; a bare CircularGraph shows `·`.
            # `ε` for the marked region, `·` for unreachable ones. `:words_num`
            # prefixes the region number (`R|word`), for the cases where the
            # picture is also referred to by region number.
            [begin
                 w = region_words === nothing ? nothing : region_words[k]
                 wtxt = w === nothing ? "·" : (isempty(w) ? "ε" : join(w))
                 cell_labels === :words ? wtxt : string(k, "|", wtxt)
             end for k in 1:region_count(g)]
        elseif cell_labels === :labels
            # `(region number, label)` — the label only if it's NOT 1. Without
            # given labels it falls back to the bare number.
            region_labels === nothing ? collect(1:region_count(g)) :
                [isone(region_labels[k]) ? string(k) :
                 string(k, ",", region_labels[k]) for k in 1:region_count(g)]
        elseif cell_labels === :labels_dist
            # `number|dDISTANCE` and, if != 1, `,label`. The distance is the
            # number `find_circular_d4_match`/`find_circular_d4_node_match` decide on;
            # without a marking (bare CircularGraph) it shows `·`.
            [begin
                 dtxt = distances === nothing ? "·" :
                        (distances[k] < 0 ? "·" : string(distances[k]))
                 base = string(k, "|d", dtxt)
                 region_labels === nothing || isone(region_labels[k]) ?
                     base : string(base, ",", region_labels[k])
             end for k in 1:region_count(g)]
        elseif distances === nothing
            fill("·", region_count(g))
        else
            [d < 0 ? "·" : string(d) for d in distances]
        end
        print(io, _circular_cell_labels_svg(g, txts;
                                       size = size, r = r,
                                       mark_between = mark_between,
                                       bottom_leaves = bottom_leaves,
                                       top_leaves = top_leaves))
    end

    # The DUAL GRAPH overlay (optionally a
    # dual graph that goes from region 0 to every other region via a tree, and
    # writes the next reflection on the edge). Drawn LAST so it sits on top of
    # the diagram; `region_parent` comes from `circular_region_tree_parents` and is
    # filled by `_circular_tutte_morphism` — it needs the marking.
    # See render/CircularRegionTree.jl.
    region_tree && region_parent !== nothing &&
        print(io, _circular_region_tree_svg(g, region_parent, region_words;
                                            letters = region_letters,
                                            size = size, r = r,
                                            mark_between = mark_between,
                                            mark_between2 = mark_between2,
                                            bottom_leaves = bottom_leaves,
                                            top_leaves = top_leaves))

    print(io, "</svg>")
    return String(take!(io))
end

# ---- 8. Display: Base.show AND Base.showable for both MIME types --------------
#
# Without an EXPLICIT `Base.showable`, VS Code falls back to text form — see the
# comment in render/Decorated.jl. These MIME methods are purely additive; the text
# `Base.show` in circular/CircularGraph.jl stands alongside them.

# The boundary marking always belongs in the image: a bare `CircularGraph` has no
# cuts, and without a marking the leaf ring would be drawn reversed
# (`rot = 0, reflect = 1`), making the wiring unreadable — exactly the case when
# displaying `term.graph` from a reduction instead of the morphism.
#
# `display(g::CircularGraph)` therefore always draws with the CANONICAL marking
# `_circular_default_mark(g)` (black, between leaf n and 1). This fixes the ring
# orientation, so `display(g)` and `display(m)` show the same image when the
# cuts sit at the same place.
#
# WHAT THE DEBUG SWITCH CONTROLS (and what it doesn't):
#   ALWAYS            — boundary marking, slot/arm numbers, LEAF NUMBERS
#   only `debug_circular()` — region numbers
#
# Leaf numbers are always on: without them you can't tell whether a node's arms
# run the leaf ring backwards, which is the house convention. Analogous to
# render/tutte/Tutte.jl for the plain case.
_circular_default_mark(g::CircularGraph) = 0    # between leaf n and leaf 1

Base.show(io::IO, ::MIME"image/svg+xml", g::CircularGraph) =
    print(io, circular_tutte_svg(g; mark_between = _circular_default_mark(g),
                            slot_numbers = true,
                            arm_numbers = _circular_debug_on(),
                            cell_numbers = _circular_debug_on(),
                            leaf_numbers = true))
Base.show(io::IO, ::MIME"image/png", g::CircularGraph) =
    write(io, _svg_to_png(circular_tutte_svg(g; mark_between = _circular_default_mark(g),
                                        slot_numbers = true,
                                        arm_numbers = _circular_debug_on(),
                                        cell_numbers = _circular_debug_on(),
                                        leaf_numbers = true)))
Base.showable(::MIME"image/svg+xml", ::CircularGraph) = true
Base.showable(::MIME"image/png", ::CircularGraph) = true

"""
    _circular_tutte_morphism(m::CircularMorphismGraph; slot_numbers = true, cell_numbers = false) -> String

Direct analogue of `_tutte_morphism` (render/tutte/Display.jl): draws `m.graph`
with both boundary markers (black = `_circular_left_mark`, grey = `_circular_right_mark`)
plus the bottom/top orientation correction (`bottom_leaves`/`top_leaves` passed
through to `circular_tutte_svg`) — bottom lands at the bottom, top at the top,
black on the left.
"""
function _circular_tutte_morphism(m::CircularMorphismGraph; slot_numbers::Bool = true,
                    cell_numbers::Bool = false, cell_labels::Symbol = :numbers,
                    leaf_numbers::Bool = false,
                    region_tree::Bool = false,
                    kwargs...)
    # Region words AND the tree come from ONE call, or the arrow labels and the
    # cell words describe different walks. The TREE words are the ones the rules
    # trigger on (`circular_unreduced_transition`), so those are what is drawn —
    # `circular_region_words`, the edges-only variant, would show a different
    # walk and, at a `:gen2` node, a chain of edge crossings where the tree takes
    # ONE step. Computed only when actually asked for.
    tw = (cell_labels in (:words, :words_num) || region_tree) ?
         circular_region_tree_parents(m) : nothing
    words  = tw === nothing ? nothing : tw.words
    parent = tw === nothing ? nothing : tw.parent
    steplet = tw === nothing ? nothing : tw.letters
    return circular_tutte_svg(m.graph; mark_between = _circular_left_mark(m),
                  mark_between2 = _circular_right_mark(m),
                  bottom_leaves = _bottom_leaves(m), top_leaves = _top_leaves(m),
                  slot_numbers = slot_numbers, cell_numbers = cell_numbers,
                  cell_labels = cell_labels, leaf_numbers = leaf_numbers,
                  # Distances are only defined here: they are measured from the
                  # black marking, which only the morphism has. Same for the
                  # region words and the tree.
                  distances = cell_labels in (:distances, :labels_dist) ?
                              circular_region_distances(m) : nothing,
                  region_words = words, region_tree = region_tree,
                  region_parent = parent, region_letters = steplet,
                  kwargs...)
end

"""
    circular_tutte_svg(m::CircularMorphismGraph) -> String

A circular morphism is ALWAYS drawn as a morphism — identical to `display(m)`.
Analogous to `tutte_svg(::MorphismGraph)`: avoids unpacking `m.graph`, which
would lose `bottom_leaves`/`top_leaves` and thus the leaf-ring orientation.
"""
circular_tutte_svg(m::CircularMorphismGraph; kwargs...) = _circular_tutte_morphism(m; kwargs...)

# The debug switch also governs `display(m)`, so both `display_circular_tutte` and
# `display`/`show` respect `debug_circular()` consistently, showing cell numbers
# whenever debug mode is on. (`debug_labels()` is checked too, since notebooks
# also set that switch — it otherwise only controls the DECORATED renderers,
# render/Decorated.jl and render/CircularCells.jl.)
_circular_debug_on() = debug_circular() || debug_labels()
# `slot_numbers = true` is `_circular_tutte_morphism`'s default and stays that way:
# arm numbers belong in the image ALWAYS, like the marking; only leaf and
# region numbers hang off the debug switch.
Base.show(io::IO, ::MIME"image/svg+xml", m::CircularMorphismGraph) =
    print(io, _circular_tutte_morphism(m; cell_numbers = _circular_debug_on(),
                                  arm_numbers = _circular_debug_on(), leaf_numbers = true))
Base.show(io::IO, ::MIME"image/png", m::CircularMorphismGraph) =
    write(io, _svg_to_png(_circular_tutte_morphism(m; cell_numbers = _circular_debug_on(),
                                             arm_numbers = _circular_debug_on(), leaf_numbers = true)))
Base.showable(::MIME"image/svg+xml", ::CircularMorphismGraph) = true
Base.showable(::MIME"image/png", ::CircularMorphismGraph) = true

"""
    display_circular_tutte(m::CircularMorphismGraph; slot_numbers = debug_circular(),
                      cell_numbers = debug_circular(), cell_labels = :numbers,
                      leaf_numbers = debug_circular())

Draws `m` as a planar Tutte embedding with bottom/top boundary markers. In
debug mode (`debug_circular(true)` at the top of a notebook), also shows arm
numbers on each leg, region numbers in each region, and the leaf number `k`
(the index in `g.word`, as in `Leaf(k)`) just outside the boundary circle next
to each boundary leaf.

With `cell_labels = :distances`, each region shows the BFS DISTANCE to the
black marking instead of the number — the number `find_circular_d4_match` decides
on. `:none` leaves the regions blank.
"""
function display_circular_tutte(m::CircularMorphismGraph;
                           slot_numbers::Bool = debug_circular(),
                           cell_numbers::Bool = debug_circular(),
                           cell_labels::Symbol = :numbers,
                           leaf_numbers::Bool = debug_circular())
    display("image/svg+xml", _circular_tutte_morphism(m;
              slot_numbers = slot_numbers, cell_numbers = cell_numbers,
              cell_labels = cell_labels, leaf_numbers = leaf_numbers))
    return nothing
end

"""
    display_circular_tutte(g::CircularGraph; slot_numbers = debug_circular(),
                      cell_numbers = debug_circular(), cell_labels = :numbers,
                      leaf_numbers = debug_circular())

Draws `g` as a planar Tutte embedding. Same debug options as the
MorphismGraph variant, except a bare `CircularGraph` has NO marking:
`cell_labels = :distances` gives `·` everywhere here. Draw the morphism, not
`m.graph`, to see distances.
"""
function display_circular_tutte(g::CircularGraph;
                           slot_numbers::Bool = debug_circular(),
                           cell_numbers::Bool = debug_circular(),
                           cell_labels::Symbol = :numbers,
                           leaf_numbers::Bool = debug_circular())
    display("image/svg+xml", circular_tutte_svg(g;
              slot_numbers = slot_numbers, cell_numbers = cell_numbers,
              cell_labels = cell_labels, leaf_numbers = leaf_numbers))
    return nothing
end

"""
    show_circular_step(fd, m0; wiring = true)

**The standard way to display a reduction step.** Draws `fd` (a `CircularDecorated`
or `CircularGraph`) as a MORPHISM with `m0`'s cuts — i.e. **with the boundary
marking** — and below it the **wiring table** (`show_wiring`: which arm hangs
on which).

Each region shows `number|dDISTANCE` and, if not `1`, the `,label` after it.
The distance is the BFS distance to the black marking — the number
`find_circular_d4_match` and `find_circular_d4_node_match` order their candidates by. A
bare `CircularGraph` has no labels, so only number and distance are shown there.
With `distances = false`, the plain `number,label` form is used.

This puts everything needed to follow a step in ONE view: marking, slot
numbers, region labels, region distances, and the wiring.
"""
function show_circular_step(fd, m0; wiring::Bool = true, distances::Bool = true)
    g = fd isa CircularGraph ? fd : fd.graph
    # Only a `CircularDecorated` carries labels; `CircularGraph` and `CircularMorphismGraph`
    # have none (then regions show only the numbers).
    labs = hasproperty(fd, :region_labels) ? fd.region_labels : nothing
    m = CircularMorphismGraph(g, m0.cut1, m0.cut2)
    svg = _circular_tutte_morphism(m; cell_numbers = true,
                              cell_labels = distances ? :labels_dist : :labels,
                              leaf_numbers = true, slot_numbers = true,
                              region_labels = labs)
    # Shown as an image in a notebook; silently skipped in a script/REPL
    # without a display stack, so the same function works in both contexts.
    try
        display(MIME"image/svg+xml"(), svg)
    catch
        println("[SVG ", length(svg), " bytes — the diagram appears here in the notebook]")
    end
    wiring && show_wiring(g)
    if distances
        dst = circular_region_distances(m)
        println("Region distances to the mark: ",
                join(["R$i=$(d < 0 ? "·" : d)" for (i, d) in enumerate(dst)], "  "))
    end
    if labs !== nothing
        nz = findall(!isone, labs)
        println(isempty(nz) ? "all region labels = 1" :
                "non-trivial labels: " *
                join(["R$i = $(labs[i])" for i in nz], ",  "))
    end
    return nothing
end
