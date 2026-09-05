# render/Cells.jl — geometry for cell/region label anchors: ring direction, arc
# midpoints, area-centroid placement, and the cosmetic node-clearance / separation /
# boundary-pull / inside-polygon repairs. Included before CellLabels.jl.

#
# 1. `display_cell_number(g)` draws the cell number into every cell (anchor =
#    centroid of the vertices on the cell boundary, Tutte positions).
# 2. `cell_distances(m)` + `display_cell_distance(m)`: for a MorphismGraph with a
#    black mark on the left (the arc between bottom and top start), the BFS
#    distance of every cell from the mark (touching cell = 0, neighbours = 1, …).
#    Two cells are adjacent when an edge separates them (port_cell left/right).
#    Unreachable cells (e.g. free circles) get -1.

"""
    _ring_ccw(leafxy, n) -> Bool

Does the leaf ring run COUNTERCLOCKWISE in the layout, i.e. is the arc from leaf
1 to leaf 2 the short way round in the mathematically positive direction?

⚠ **Why this has to be read off, not assumed** (the same reasoning as in
[`_pull_boundary_anchors!`](@ref)):
`_halfcircle_angles` assigns the leaf angles DECREASING from π, so the ring
usually runs CLOCKWISE. A fixed `mod(θb - θa, 2pi)` then measures the LONG way
around the circle, and an arc midpoint computed that way lands on the OPPOSITE
side of the disc.
"""
function _ring_ccw(leafxy, n::Int)
    n < 2 && return true
    θ1 = atan(leafxy[1][2], leafxy[1][1])
    θ2 = atan(leafxy[2][2], leafxy[2][1])
    return mod(θ2 - θ1, 2pi) < pi
end

"""
    _gap_mid_point(leafxy, n, gg) -> Tuple{Float64,Float64}

Midpoint of the boundary arc of gap `gg` (between leaf `gg` and leaf `gg+1`) on
the unit circle, measured IN THE RUNNING DIRECTION of the ring
([`_ring_ccw`](@ref)).

**THE DIRECTION MATTERS.** A fixed `mod(θb - θa, 2pi)` is the COMPLEMENT arc on a
clockwise ring, so every arc midpoint comes out ANTIPODAL, and with it the cell
polygon and its centroid. On `id_1213212` (`light_leaf_up`, eight strips) that puts
region 1, the far-LEFT sliver, at `(+0.06, 0)` and region 8, the far-RIGHT one, at
`(−0.06, 0)` — both in the middle
and mirrored. With the running direction they land at ∓0.95, i.e. in their own
strip.
"""
function _gap_mid_point(leafxy, n::Int, gg::Int)
    pa = leafxy[gg]; pb = leafxy[mod1(gg + 1, n)]
    θa = atan(pa[2], pa[1]); θb = atan(pb[2], pb[1])
    Δ  = _ring_ccw(leafxy, n) ? mod(θb - θa, 2pi) : -mod(θa - θb, 2pi)
    θm = θa + Δ / 2
    return (cos(θm), sin(θm))
end

# ---- anchor positions (centroid of the boundary vertices of each cell) -------

# Anchor = area centroid of the cell polygon (walk vertices in order, plus arc
# midpoints between consecutive leaves). The naive vertex mean lies ON the strand
# for thin cells (e.g. a dot); the polygon centroid lies inside and touches no
# lines.
"""
    _cell_anchors_from(walks, cell_of_face, inner, leafxy, nodexy, n) -> Dict{Int,Tuple}

The anchor geometry behind `_circular_cell_anchors` (render/CircularCells.jl):
given the face walks, the face -> cell map, the inner faces, the Tutte positions
and the boundary length, place one anchor per cell at the area centroid of its
polygon.

Polygon points are the walk vertices, plus the arc midpoint of a boundary gap
between two consecutive leaves — but only when the gap belongs to THIS cell,
otherwise the centroid lands on the strand (e.g. for a chord). Falls back to the
vertex mean for a degenerate (zero-area) polygon.
"""
function _cell_anchors_from(walks, cell_of_face, inner, leafxy, nodexy, n::Int)
    # gaps per cell (only the gaps OF the cell get an arc midpoint — otherwise
    # the centroid lands on the strand, e.g. for a chord).
    gaps_of = Dict(cell.id => Set(cell.gaps) for cell in inner)
    anchors = Dict{Int, Tuple{Float64, Float64}}()
    for (fi, walk) in enumerate(walks)
        c = cell_of_face[fi]
        c == 0 && continue
        mygaps = get(gaps_of, c, Set{Int}())
        pos(v) = v[1] === :leaf ? leafxy[v[2]] : nodexy[v[2]]
        # polygon points: walk vertices; between two consecutive leaves insert
        # the arc midpoint of their boundary gap — but only if the gap belongs to
        # this cell (the other transition is the strand).
        pts = Tuple{Float64, Float64}[]
        m = length(walk)
        for i in 1:m
            v = walk[i]
            push!(pts, pos(v))
            w = walk[mod1(i + 1, m)]
            (v[1] === :leaf && w[1] === :leaf) || continue
            a, b = v[2], w[2]
            if a == b
                # a single leaf: the gap runs all the way round, arc midpoint =
                # the antipode of the leaf
                a in mygaps || continue
                p = leafxy[a]
                push!(pts, (-p[1], -p[2]))
            else
                # gap between a and b: g = a if b = a+1 (cyclically), otherwise
                # g = b (n = 2: both gaps connect the same pair, the direction in
                # the walk decides).
                gg = mod1(a + 1, n) == b ? a : b
                gg in mygaps || continue
                push!(pts, _gap_mid_point(leafxy, n, gg))
            end
        end
        isempty(pts) && continue
        # area centroid (shoelace); fallback: vertex mean.
        A = Cx = Cy = 0.0
        for i in 1:length(pts)
            (x1, y1) = pts[i]; (x2, y2) = pts[mod1(i + 1, length(pts))]
            cr = x1 * y2 - x2 * y1
            A += cr; Cx += (x1 + x2) * cr; Cy += (y1 + y2) * cr
        end
        if abs(A) > 1e-9
            anchors[c] = (Cx / (3A), Cy / (3A))
        else
            anchors[c] = (sum(p[1] for p in pts) / length(pts),
                          sum(p[2] for p in pts) / length(pts))
        end
    end
    _pull_boundary_anchors!(anchors, gaps_of, leafxy, n)
    return _push_anchors_off_nodes!(anchors, nodexy)
end


"""
    _push_anchors_off_nodes!(anchors, nodexy; clear = 0.10) -> anchors

Pushes every cell anchor away from the nodes that lie too close: in debug mode a
tuple such as `(1, 0)` would otherwise sit in the middle of a braid node and be
unreadable.

WHY THIS HAPPENS AT ALL: the anchor is the area centroid of the cell polygon. If a
cell encloses a node almost like a ring (the walk passes both sides of the same
node), its centroid falls exactly on that node — geometrically right and still
unreadable. Without debug mode it hardly showed, because the most common label `1`
is not drawn at all; the tuple always is.

The repair is purely cosmetic (only the LABEL moves, not the cell): if the anchor
is closer than `clear` to a node, it is pushed radially away from it to distance
`clear`. Iterated, so that dodging does not run straight into the next node; if it
gets stuck the last state wins — a slightly misplaced label is better than an
infinite loop.
"""
function _push_anchors_off_nodes!(anchors::Dict{Int,Tuple{Float64,Float64}},
                                  nodexy; clear::Float64 = 0.10)
    isempty(nodexy) && return anchors
    for (c, (ax, ay)) in anchors
        x, y = ax, ay
        for _ in 1:8
            moved = false
            for (_, (nx, ny)) in nodexy
                dx, dy = x - nx, y - ny
                d = hypot(dx, dy)
                d < clear || continue
                # exactly on the node: there is no dodge direction, so pick a
                # fixed one (upwards) instead of dividing by 0.
                (ux, uy) = d < 1e-9 ? (0.0, 1.0) : (dx / d, dy / d)
                x = nx + clear * ux
                y = ny + clear * uy
                moved = true
            end
            moved || break
        end
        anchors[c] = (x, y)
    end
    return anchors
end

"""
    _separate_anchors!(anchors; clear = 0.12) -> anchors

Pulls apart label anchors that (almost) coincide.

WHY: the anchor is the centroid of the polygon of WALK VERTICES. Two different
cells can have the same vertex set — with several PARALLEL edges between the same
two nodes, every bigon cell in between has the walk `[n_a, n_b]`. E.g. for
`circular_double_leaf(212321, [1,1,1,1,1,1]²)`: cells 3 and 4 (both bigons between n1
and n2) get EXACTLY the same anchor, and their labels land on top of each other and
cannot be read.

The vertices alone cannot tell the cells apart — that would need the renderer's
edge geometry. Instead the pragmatic repair: anchors closer than `clear` are
pushed apart (in a fixed direction if they coincide exactly, otherwise along the
line joining them). Purely cosmetic, the cells stay where they are.
"""
function _separate_anchors!(anchors::Dict{Int,Tuple{Float64,Float64}};
                            clear::Float64 = 0.12)
    ids = sort(collect(keys(anchors)))
    length(ids) < 2 && return anchors
    for _ in 1:12
        moved = false
        for i in 1:length(ids), j in (i + 1):length(ids)
            (ax, ay) = anchors[ids[i]]; (bx, by) = anchors[ids[j]]
            dx, dy = bx - ax, by - ay
            d = hypot(dx, dy)
            d < clear || continue
            # exactly on top of each other: there is no separating direction, so
            # pick a fixed one (horizontal) so they come apart at all.
            (ux, uy) = d < 1e-9 ? (1.0, 0.0) : (dx / d, dy / d)
            shift = (clear - d) / 2 + 1e-3
            anchors[ids[i]] = (ax - shift * ux, ay - shift * uy)
            anchors[ids[j]] = (bx + shift * ux, by + shift * uy)
            moved = true
        end
        moved || break
    end
    return anchors
end

"""
    _pull_boundary_anchors!(anchors, gaps_of, leafxy, n; frac = 0.72) -> anchors

Pulls the anchor of a BOUNDARY cell outwards, towards its own boundary gaps
(for boundary cells the label may sit closer to
the edge).

WHY: the anchor is the area centroid of the cell polygon. A boundary cell that
touches several arcs on OPPOSITE sides of the disc has its centroid almost in the
middle — e.g. for `circular_double_leaf(12321, [1,0,1,0,0]²)` the two boundary cells
#1 and #2 sit at radius 0.15 resp. 0.19, close together in the centre, even
though they occupy the whole periphery.

The repair is purely cosmetic (only the LABEL moves): the label is put at `frac`
of the way to a boundary point of its own gaps — outwards to its own edge, but not
onto the boundary line itself. Cells without gaps (interior cells) are unchanged.

WHICH boundary point (the labels don't sit at the
right spot). Two cases, and the difference is NOT cosmetic:

* The gaps of the region form ONE CONNECTED arc (cyclically; the seam may lie in
  the middle of it). The label sits in its MIDDLE — e.g. for
  `circular_double_leaf(121, [1,0,0]²)`, region 1 covers gaps 1 and 2 (arc
  midpoints at 120° and 180°), so the label sits at 150°.

* The gaps form SEVERAL separate arcs. That is not a bug but the band shape:
  `circular_double_leaf(121, [0,1,1]²)` has a region with gaps {2,4} — the band between
  the chords `Leaf5—Leaf2` and `Leaf4—Leaf3`, touching the boundary at two OPPOSITE
  places. A bisector would then lie in the middle of the disc, outside the region.
  Here we pull to the LONGEST single arc — the widest part of the region, where the
  label is most likely to fit.

Cells without gaps (interior cells) are unchanged; their centroid already lies
inside.
"""
function _pull_boundary_anchors!(anchors::Dict{Int,Tuple{Float64,Float64}},
                                 gaps_of, leafxy, n::Int; frac::Float64 = 0.72,
                                 polys = nothing)
    n == 0 && return anchors
    # ORIENTATION OF THE PICTURE (the x-coordinate is
    # sometimes swapped left and right).
    #
    # `_ring_transform` (render/tutte/Markers.jl) keeps `reflect` at 1, so the leaf
    # ring is never mirrored. With a mirrored ring THE LEAF ANGLES RUN BACKWARDS: for
    # `circular_double_leaf(121, [0,1,1]²)` leaves 1..6 sat at 143°, 90°, 36.9°, 323.1°,
    # 270°, 216.9° — decreasing instead of increasing. The direction is not fixed a
    # priori either (`_halfcircle_angles` assigns θ DECREASING from π), so it has to
    # be read off below.
    #
    # Both formulas below took a fixed `mod(…, 2pi)`, i.e. the arc COUNTERCLOCKWISE
    # from `k` to `k+1`. In the mirrored picture that is the LONG way round, and
    # the "middle" lands exactly opposite: gap 1 (between 143° and 90°, middle
    # 116.6°) came out as 296.6°. That is precisely the left/right swap.
    #
    # Repair: READ the direction off the layout instead of assuming it. The ring is
    # evenly spaced, so a single step decides.
    ccw = n < 2 ? true : begin
        θ1 = atan(leafxy[1][2], leafxy[1][1])
        θ2 = atan(leafxy[2][2], leafxy[2][1])
        mod(θ2 - θ1, 2pi) < pi          # does 1→2 run counterclockwise?
    end
    # arc from θa to θb IN THE RUNNING DIRECTION, as a signed difference.
    arc(θa, θb) = ccw ? mod(θb - θa, 2pi) : -mod(θa - θb, 2pi)
    # arc midpoint of gap k (between leaf k and k+1), as an angle.
    function gap_angle(k)
        pa = leafxy[k]; pb = leafxy[mod1(k + 1, n)]
        θa = atan(pa[2], pa[1]); θb = atan(pb[2], pb[1])
        return θa + arc(θa, θb) / 2
    end
    for (c, (ax, ay)) in anchors
        gaps = get(gaps_of, c, Set{Int}())
        isempty(gaps) && continue           # interior cell: nothing to do
        # split the gaps into cyclically connected runs.
        s = sort(collect(gaps))
        runs = Vector{Vector{Int}}()
        for k in s
            if !isempty(runs) && mod(k - last(runs[end]), n) == 1
                push!(runs[end], k)
            else
                push!(runs, [k])
            end
        end
        # seam: if the last run runs into the first, merge them — in the RUNNING
        # DIRECTION (last run first), so that `run[1]`/`run[end]` remain start and
        # end in that direction. `prepend!` would have given [1,5,6] instead of
        # [5,6,1] and put the arc middle on the opposite side.
        if length(runs) > 1 && mod(first(runs[1]) - last(runs[end]), n) == 1
            runs[end] = vcat(runs[end], runs[1])
            popfirst!(runs)
        end
        # take the longest run and aim at its middle. With exactly one run that is
        # the arc middle of the whole region.
        #
        # The run is already built in the RUNNING DIRECTION (seam included), so
        # `run[1]` → `run[end]` — do NOT take the sorted order: for {1,5,6} the run
        # is 5,6,1 and the middle sits at gap 6, whereas 1→6 would point to the
        # opposite side.
        # THE PULL MUST NOT LEAVE THE CELL (the region
        # label is not always visually inside its cell, but sometimes
        # displaced). With SEVERAL runs the pull aims at the longest
        # arc, and for a band that touches the boundary on opposite sides that
        # target can lie in a different cell — e.g. `id_1213212`, whose eight
        # strips each have two opposite gaps: without this check, their labels
        # collapse onto the two outer strips. When `polys` is passed, the runs
        # are tried longest-first and the FIRST target that lies inside the cell
        # is used; if none does, the centroid stays.
        poly = polys === nothing ? nothing : get(polys, c, nothing)
        ordered = sort(runs; by = length, rev = true)
        if poly !== nothing
            hit = nothing
            for rr in ordered
                θa0 = gap_angle(rr[1])
                θm0 = θa0 + arc(θa0, gap_angle(rr[end])) / 2
                q = (frac * cos(θm0), frac * sin(θm0))
                if _point_in_polygon(q, poly)
                    hit = q; break
                end
            end
            hit === nothing || (anchors[c] = hit)
            continue
        end
        run = runs[argmax(length.(runs))]
        θa = gap_angle(run[1])
        # here too IN THE RUNNING DIRECTION (see `arc` above) — in the mirrored
        # picture `mod(…, 2pi)` pointed to the opposite side: region 1 with gaps
        # {1,5,6} (run 5,6,1) ended up at 180° instead of 0°.
        θm = θa + arc(θa, gap_angle(run[end])) / 2
        # target point at `frac` of the radius, so the label stays INSIDE the
        # disc. Before, it was pulled from the centroid to `frac` of the WAY to the
        # boundary point — if the centroid already sits far out (band cells) the
        # label ended up at r ≈ 0.98, i.e. on the boundary line itself.
        p = (frac * cos(θm), frac * sin(θm))
        anchors[c] = p
    end
    return anchors
end

"""
    _point_in_polygon(p, poly) -> Bool

Ray casting: is `p` inside the closed polygon `poly` (a vector of `(x, y)`)?
Points exactly on an edge count as inside up to the ray's rounding — good
enough for placing a label.
"""
function _point_in_polygon(p::Tuple{Float64,Float64},
                           poly::Vector{Tuple{Float64,Float64}})
    length(poly) >= 3 || return false
    (x, y) = p
    inside = false
    j = length(poly)
    for i in 1:length(poly)
        (xi, yi) = poly[i]; (xj, yj) = poly[j]
        if (yi > y) != (yj > y)
            xc = xi + (y - yi) * (xj - xi) / (yj - yi)
            xc > x && (inside = !inside)
        end
        j = i
    end
    return inside
end

"Distance from `p` to the closest edge of `poly` (0 if `poly` is degenerate)."
function _dist_to_polygon(p::Tuple{Float64,Float64},
                          poly::Vector{Tuple{Float64,Float64}})
    length(poly) >= 2 || return 0.0
    (px, py) = p
    best = Inf
    j = length(poly)
    for i in 1:length(poly)
        (ax, ay) = poly[j]; (bx, by) = poly[i]
        dx = bx - ax; dy = by - ay
        L2 = dx * dx + dy * dy
        t = L2 < 1e-12 ? 0.0 : clamp(((px - ax) * dx + (py - ay) * dy) / L2, 0.0, 1.0)
        best = min(best, hypot(px - (ax + t * dx), py - (ay + t * dy)))
        j = i
    end
    return best
end

"""
    _anchor_inside!(anchors, polys) -> anchors

Makes sure every label sits **in its own cell** (the
region label is not always visually inside its cell, but sometimes
displaced).

The anchor is the AREA CENTROID of the cell polygon, and for a non-convex cell
— an L, a horseshoe, a cell wrapping around a node — the centroid can lie
**outside** the cell. It then reads as the label of the neighbouring cell. The
same happens after `_pull_boundary_anchors!`/`_push_anchors_off_nodes!`, which
shift anchors for cosmetic reasons without knowing the polygon.

This pass runs LAST and only touches anchors that actually fell out: it picks,
among the polygon's vertex/edge midpoints and their pairwise midpoints, the
interior candidate FARTHEST from the boundary — a cheap stand-in for the pole
of inaccessibility. Anchors already inside stay exactly where they were, so
nothing that looked right before moves.
"""
function _anchor_inside!(anchors::Dict{Int,Tuple{Float64,Float64}},
                         polys::Dict{Int,Vector{Tuple{Float64,Float64}}},
                         centroids::Dict{Int,Tuple{Float64,Float64}} =
                             Dict{Int,Tuple{Float64,Float64}}())
    for (id, poly) in polys
        haskey(anchors, id) || continue
        length(poly) >= 3 || continue
        _point_in_polygon(anchors[id], poly) && continue
        # First fall back to the UNMOVED centroid — for a cell whose centroid was
        # fine and only the cosmetic steps pushed the label out, that is exactly
        # the right place and costs no search.
        if haskey(centroids, id) && _point_in_polygon(centroids[id], poly)
            anchors[id] = centroids[id]
            continue
        end
        # candidates: the vertices pulled slightly inward, the edge midpoints,
        # and all pairwise midpoints of those.
        base = Tuple{Float64,Float64}[]
        c = (sum(p[1] for p in poly) / length(poly), sum(p[2] for p in poly) / length(poly))
        for (k, p) in enumerate(poly)
            q = poly[mod1(k + 1, length(poly))]
            push!(base, ((p[1] + q[1]) / 2, (p[2] + q[2]) / 2))
            push!(base, (p[1] + 0.25 * (c[1] - p[1]), p[2] + 0.25 * (c[2] - p[2])))
        end
        cand = copy(base)
        for a in 1:length(base), b in (a + 1):length(base)
            push!(cand, ((base[a][1] + base[b][1]) / 2, (base[a][2] + base[b][2]) / 2))
        end
        # NEAREST interior candidate, not the deepest one: the anchor that the
        # steps above produced is where the label should sit (it already
        # respects the boundary pull and the node clearance). Depth only
        # breaks ties.
        old  = anchors[id]
        best = nothing; bestk = (Inf, 0.0)
        for p in cand
            _point_in_polygon(p, poly) || continue
            k = (hypot(p[1] - old[1], p[2] - old[2]), -_dist_to_polygon(p, poly))
            k < bestk && (best = p; bestk = k)
        end
        best === nothing || (anchors[id] = best)
    end
    return anchors
end

"""
    _region_anchors(regs, walks, cell_of_face, leafxy, nodexy, n) -> Dict{Int,Tuple}

Anchor position per REGION (not per face). The renderer labels regions: a face that
touches the boundary at several separate places gets its own number at EACH of them,
instead of a single one somewhere.

The construction is the area centroid of the cell
polygon — only restricted to the gaps OF THE REGION: the walk of the corresponding
face supplies the vertices, but an arc midpoint enters the polygon only if its gap
belongs to this region. That way the centroid lands in the part of the face the
region actually occupies.

Regions without a boundary gap (interior cells: lens, free circle) use the walk
unchanged — there is nothing to restrict there.
"""
function _region_anchors(regs::Vector{Region}, walks, cell_of_face,
                         leafxy, nodexy, n::Int)
    walk_of_cell = Dict{Int, Any}()
    for (fi, walk) in enumerate(walks)
        c = cell_of_face[fi]
        c == 0 && continue
        walk_of_cell[c] = walk
    end
    anchors = Dict{Int, Tuple{Float64, Float64}}()
    gaps_of = Dict{Int, Set{Int}}()
    # The cell polygon per region, kept for the final `_anchor_inside!` pass.
    polys   = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    for R in regs
        haskey(walk_of_cell, R.cell) || continue
        walk = walk_of_cell[R.cell]
        mygaps = Set(R.gaps)
        gaps_of[R.id] = mygaps
        pos(v) = v[1] === :leaf ? leafxy[v[2]] : nodexy[v[2]]
        pts = Tuple{Float64, Float64}[]
        m = length(walk)
        for i in 1:m
            v = walk[i]
            push!(pts, pos(v))
            w = walk[mod1(i + 1, m)]
            (v[1] === :leaf && w[1] === :leaf) || continue
            a, b = v[2], w[2]
            if a == b
                a in mygaps || continue
                p = leafxy[a]
                push!(pts, (-p[1], -p[2]))
            else
                gg = mod1(a + 1, n) == b ? a : b
                gg in mygaps || continue
                push!(pts, _gap_mid_point(leafxy, n, gg))
            end
        end
        isempty(pts) && continue
        polys[R.id] = pts
        A = Cx = Cy = 0.0
        for i in 1:length(pts)
            (x1, y1) = pts[i]; (x2, y2) = pts[mod1(i + 1, length(pts))]
            cr = x1 * y2 - x2 * y1
            A += cr; Cx += (x1 + x2) * cr; Cy += (y1 + y2) * cr
        end
        if abs(A) > 1e-9
            anchors[R.id] = (Cx / (3A), Cy / (3A))
        else
            anchors[R.id] = (sum(p[1] for p in pts) / length(pts),
                             sum(p[2] for p in pts) / length(pts))
        end
    end
    # the untouched centroids, as the fallback for `_anchor_inside!`
    centroids = copy(anchors)
    _pull_boundary_anchors!(anchors, gaps_of, leafxy, n; polys = polys)
    _push_anchors_off_nodes!(anchors, nodexy)
    # separate coinciding labels (parallel edges ⇒ same vertex set ⇒ same
    # centroid). After this do NOT push away from the nodes again, or the two run
    # into each other once more.
    _separate_anchors!(anchors)
    # LAST: pull back into the own cell whatever the steps above pushed out, and
    # whatever the centroid of a non-convex cell got wrong to begin with
    # Only touches anchors that really lie outside.
    return _anchor_inside!(anchors, polys, centroids)
end
