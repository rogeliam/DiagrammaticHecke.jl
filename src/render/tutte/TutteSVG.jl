# render/tutte/TutteSVG.jl — tutte_svg(g): draws a WordGraph's Tutte embedding as an
# SVG (boundary leaves, internal nodes, edges with fan/repel/braid-angle handling,
# slot numbers, markers). Included after Tutte.jl.


"""
    tutte_svg(g; size = 600, r = 250) -> String

Draw `g` with a planar Tutte embedding (see `tutte_positions`): boundary leaves as
labelled discs on the outer circle, internal nodes at their barycentres, coloured
edges between them. A crossing-free ALTERNATIVE to `shell_svg`. Graphs WITHOUT a
boundary (barbell, …) work too; a completely empty graph draws as the empty
diagram (dashed circle only). `slot_numbers` (default `true`, as in
`circular_tutte_svg`) writes the 1-based slot number just outside each non-dot
node's rim on every leg — a debug aid for spotting wiring errors in
hand-built example graphs, where colours alone don't reveal which slot goes
where; switch off with `slot_numbers = false`.

For the full WIRING (which slot goes where) use `wiring_table(g)` — a table
beside the picture rather than inside it (printing `s→target`
on every leg made the drawing too busy).

`leaf_numbers` (default: the debug switch `_plain_debug_on()`, i.e.
`debug_labels()`/`debug_circular()`) additionally writes the LEAF NUMBER `k` (the index
in `g.word`, as in `Leaf(k)`) just outside the boundary circle next to every leaf
disc — the colour stays IN the disc. Without that number a picture does not say
which leaf is which, and that is exactly why the wiring errors in the hand-wired
examples are hard to spot without it. In the same debug
mode the boundary gets its SEAM marker: a small red ring in gap 0, i.e. at the
transition from the last to the first leaf. It uses the same marker mechanism as
the morphism (`MarkSpec` + `_marker_pt`), but is not passed as `mark_between`:
`mark_between` rotates the whole leaf ring via `_ring_transform`, and the
debug picture has the same layout as the normal one.
`edge_repel` (default `true`) makes edges AVOID foreign nodes instead
of running straight over them (render/tutte/Avoid.jl). At legs of a `:braid`
node this happens only within tight angle limits, so the slot order is
preserved. `edge_repel = false` draws every edge straight.

`dot_gap` (default `12` px = one dot diameter) is the clearance between the
edge of a dot and the edge of its neighbour (node or leaf disc).
"""
function tutte_svg(g::WordGraph; size::Int = 600, r::Real = 0.34 * size,
                   mark_between::MarkSpec = nothing,
                   mark_between2::MarkSpec = nothing,
                   bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                   top_leaves::Union{Nothing,Vector{Int}} = nothing,
                   slot_numbers::Bool = true,
                   leaf_numbers::Bool = _plain_debug_on(),
                   repel::Real = 0.02, iters::Int = 400,
                   edge_repel::Bool = true, dot_gap::Real = 2 * _DOT_R)
    # Boundary-LESS diagrams (barbell, …) are drawable as well: the Tutte layout
    # needs no pinned leaves — internal nodes seed on a small circle and spread by
    # the repulsion. A COMPLETELY empty graph is the empty diagram: just the
    # dashed outer circle — UNLESS a marker was requested (a MorphismGraph at n==0,
    # `identity_morphism()`, still has two well-defined symbolic markers at π/0 on
    # the dashed circle, see `_left_mark`/`_right_mark`) — in that case fall through to
    # the normal path, which draws no leaves/nodes/edges (all the loops below are
    # simply empty) but still places the marker dots. A free monochrome circle
    # (`Edge(c, Circle(c), Circle(c))`, no leaves or nodes otherwise) must ALSO fall
    # through, or the early return would draw only the dashed boundary and skip the
    # free-circle ring itself.
    isempty(g.word) && isempty(g.nodes) && isempty(g.edges) &&
        mark_between === nothing && mark_between2 === nothing &&
        return _empty_diagram_svg(size)
    leafxy, nodexy = circular_tutte_positions(circular(g); mark_between = mark_between, mark_between2 = mark_between2,
                                         bottom_leaves = bottom_leaves, top_leaves = top_leaves,
                                         repel = repel, iters = iters, dot_gap = dot_gap)
    n = length(g.word)
    baseθ(k) = pi/2 + 2pi * (k - 1) / max(n, 1)
    rot, reflect = _ring_transform(baseθ, n, mark_between, bottom_leaves, top_leaves)
    cx = cy = size / 2
    # angle of a marker's gap under the SAME rot/reflect used for the leaves —
    # computed directly (not re-derived from leafxy), so it works even when the
    # marker's flanking leaves coincide or the boundary is empty.
    hc_sv = _halfcircle_angles(n, bottom_leaves, top_leaves)
    # the same case distinction as in `tutte_positions` — otherwise the marker is
    # drawn somewhere other than where the leaves are
    structured_sv = !(mark_between !== nothing &&
                      !(mark_between isa Int && mark_between == 0) &&
                      (bottom_leaves === nothing || top_leaves === nothing ||
                       isempty(bottom_leaves) || isempty(top_leaves)))
    function _marker_pt(mb::MarkSpec)
        # Under HALF-CIRCLE assignment the leaf angles are not derivable from
        # `baseθ`, so the marker is determined from the actual flank angles.
        # Symbolic specs (:left/:right, n == 0) have no flanks and go through
        # `_mark_angle`.
        if structured_sv && !(mb isa Symbol)
            θh = _halfcircle_marker_angle(hc_sv, n, mb)
            θh !== nothing && return (cx + r * cos(θh), cy + r * sin(θh))
        end
        θraw = _mark_angle(baseθ, n, mb)
        θraw === nothing && return nothing
        # a bare :left/:right convention point (n == 0 only) is the FINAL angle
        # already (no ring to rotate — `rot` is always 0 there anyway; `reflect`
        # is always 1, see `_ring_transform`).
        # Every other spec (gap Int, `(:opposite,c)`, leaf pair) follows the ring's
        # rot/reflect like a normal leaf angle — `(:opposite, c)` in particular MUST
        # pick up `rot` so that it stays exactly π from gap c's own marker after
        # rotation, which is the whole point of the antipodal form. `cap_morphism(1)`
        # at n=2 is the case that shows it.
        θ = mb isa Symbol ? θraw : reflect * (θraw + rot)
        return (cx + r * cos(θ), cy + r * sin(θ))
    end
    X(p) = cx + r * p[1];  Y(p) = cy + r * p[2]
    pos(port::Leaf) = (X(leafxy[port.k]), Y(leafxy[port.k]))
    pos(port::NodePort) = (X(nodexy[port.node]), Y(nodexy[port.node]))

    # departure angles for :braid nodes only (see _braid_leg_angles docstring) —
    # in unit-disc coordinates, unaffected by the later isotropic X/Y scaling.
    θangles = _braid_leg_angles(g, leafxy, nodexy)
    braid_legdir(v, s) = (cos(θangles[v][s]), sin(θangles[v][s]))

    io = IOBuffer()
    print(io, """<svg width="$size" height="$size" viewBox="0 0 $size $size" xmlns="http://www.w3.org/2000/svg" font-family="DejaVu Sans Mono, Consolas, Liberation Mono, Courier New, monospace">""")
    print(io, """<circle cx="$cx" cy="$cy" r="$r" fill="none" stroke="#e2e7ee" stroke-dasharray="3" stroke-width="1"/>""")
    # edges — several edges between the SAME two endpoints would overdraw as one line,
    # so we FAN them: bundle edges by endpoint pair and bow each one out by a symmetric
    # offset (the middle one stays straight).
    drawable = [e for e in g.edges if !(e.a isa Circle || e.b isa Circle)]
    endkey(p::Leaf)     = (:leaf, p.k)
    endkey(p::NodePort) = (:node, p.node)
    bundle = Dict{Tuple,Vector{Int}}()             # unordered endpoint-pair → edge indices
    for (i, e) in enumerate(drawable)
        k = Tuple(sort([endkey(e.a), endkey(e.b)]))
        push!(get!(bundle, k, Int[]), i)
    end
    for (key, idxs) in bundle
        # ---- ordering + fan-side of parallel (same-endpoint-pair) edges -----------
        #
        # When a braid host has TWO parallel edges to the SAME neighbour
        # (necessarily consecutive slots, since a braid alternates s,t,s,t,… and a
        # repeated neighbour can only recur at slot+1), ordering and side-picking
        # those edges by ASCENDING COLOUR NUMBER has no relation to the node's
        # actual cyclic slot order. The geometric "side" (perpendicular offset sign)
        # must match the angular position of the tied slots' NEIGHBOURING (non-tied)
        # slots at BOTH ends, and colour order can contradict that: with slots 1,2 →
        # node1 and slots 5,6 → node6, a colour-order fan places slot2 and slot6
        # (both colour 3) angularly adjacent, breaking the required alternation.
        #
        # SO: if either endpoint is a `:braid` node, use that node's SLOT index
        # (not colour) as the tie-break order, and orient the perpendicular fan axis
        # (`px,py`) so that ascending slot number walks from the tied group's "slot
        # below" neighbour towards its "slot above" neighbour — i.e. the same
        # direction the node's OWN already-placed neighbours already use. This is
        # colour- and rotation-generic (works for any colour pair, any m, any slot
        # position of the tie) since it only ever consults slot arithmetic and the
        # neighbour's real drawn position, never a colour number.
        colour_order() = begin
            bycol = Dict{Int,Vector{Int}}()
            for i in idxs
                push!(get!(bycol, drawable[i].colour, Int[]), i)
            end
            cols = sort(collect(keys(bycol)))
            maxlen = maximum(length(q) for q in values(bycol))
            out = Int[]
            for pos in 1:maxlen, c in cols
                q = bycol[c]
                pos <= length(q) && push!(out, q[pos])
            end
            out
        end
        # locate a :braid endpoint (own slot) common to every edge in this bundle —
        # only meaningful when the bundle has >= 2 edges (a genuine parallel tie).
        braidslot(e) = begin
            if e.a isa NodePort && g.nodes[e.a.node].kind === :braid
                return (e.a.node, e.a.slot)
            elseif e.b isa NodePort && g.nodes[e.b.node].kind === :braid
                return (e.b.node, e.b.slot)
            end
            return nothing
        end
        refnode = nothing
        slotmap = Dict{Int,Int}()          # edge index -> slot at refnode
        if length(idxs) >= 2
            cands = [braidslot(drawable[i]) for i in idxs]
            if all(!isnothing, cands)
                # prefer the endpoint node shared by ALL edges' braidslot() pick;
                # since both e.a/e.b are fixed per edge, just require the SAME node
                # index across the bundle (true whenever both ends are the tied pair).
                nodes_seen = unique(c[1] for c in cands)
                if length(nodes_seen) == 1
                    refnode = nodes_seen[1]
                    for (i, c) in zip(idxs, cands)
                        slotmap[i] = c[2]
                    end
                end
            end
        end
        if refnode !== nothing
            idxs = sort(idxs, by = i -> slotmap[i])
        else
            idxs = colour_order()
        end
        # NOTE: named `nbundle`, NOT `n`. Naming it `n` would shadow the outer `n`
        # (the leaf count) with the bundle size, and `_marker_pt` (below, in the SAME
        # function scope) reads `n`: Julia `for` bodies do not get a fresh binding
        # for a name already local to the enclosing function, so the marker code
        # would see the LAST bundle's size instead of the leaf count and both
        # boundary markers would vanish.
        nbundle = length(idxs)
        # symmetric offsets centred on 0: e.g. nbundle=3 → [-1,0,1], nbundle=2 → [-0.5,0.5].
        # IMPORTANT: the perpendicular offset direction must be CANONICAL for the pair —
        # independent of each edge's stored a/b order — otherwise two parallel edges
        # stored as n1→n3 and n3→n1 flip `dx,dy` (hence `px,py`), the offsets cancel, and
        # both bows land on the SAME side and overdraw, so that n=2 looks like fewer
        # edges. Derive px,py once from the sorted endpoint pair.
        posof(k) = k[1] === :leaf ? (X(leafxy[k[2]]), Y(leafxy[k[2]])) :
                                    (X(nodexy[k[2]]), Y(nodexy[k[2]]))
        (kx1, ky1) = posof(key[1]); (kx2, ky2) = posof(key[2])
        kdx = kx2 - kx1; kdy = ky2 - ky1; kL = hypot(kdx, kdy)
        px = kL > 1e-9 ? -kdy / kL : 0.0
        py = kL > 1e-9 ?  kdx / kL : 0.0
        # orient (px,py) so ascending slot walks from "slot below the tie" to "slot
        # above the tie" at refnode, matching that node's own drawn neighbours —
        # this is what pins the fan side to the true cyclic order instead of an
        # arbitrary colour-derived one.
        #
        # Comparing the below/above neighbours' RAW direction vectors' projection
        # onto (px,py) does NOT work: it breaks whenever below/above sit roughly
        # PERPENDICULAR to the fan axis (an m=2 back-to-back braid pair whose two
        # hosts are drawn almost vertically aligned gives a nearly horizontal
        # (px,py) while below/above sit almost directly ABOVE the host, and the
        # projection's sign becomes numerically unstable there). Instead of
        # comparing abstract direction projections, actually compute
        # each tied slot's DEPARTURE angle from `refnode` under both candidate fan
        # signs, and pick whichever sign makes the tied slots' departure angles
        # fall, in ascending-slot order, on the SAME (shorter, non-wrapping) arc
        # from `below`'s angle to `above`'s angle — i.e. genuinely test the
        # resulting picture instead of an indirect proxy for it.
        if refnode !== nothing
            refnd = g.nodes[refnode]
            ns = 2 * refnd.m
            (rx, ry) = (X(nodexy[refnode]), Y(nodexy[refnode]))
            tiedslots = sort(collect(values(slotmap)))
            below = mod1(first(tiedslots) - 1, ns)
            above = mod1(last(tiedslots) + 1, ns)
            at_ref = _edges_at_node(g, refnode)
            farof = Dict{Int,Port}()
            for (_, own, far) in at_ref
                own isa NodePort || continue
                farof[own.slot] = far
            end
            angleof(far) = far isa Leaf ? atan(Y(leafxy[far.k]) - ry, X(leafxy[far.k]) - rx) :
                                           atan(Y(nodexy[far.node]) - ry, X(nodexy[far.node]) - rx)
            if haskey(farof, below) && haskey(farof, above) && below != above
                θbelow = angleof(farof[below]); θabove = angleof(farof[above])
                # angular sweep from below to above, going the way that does NOT
                # pass through any OTHER already-placed slot (there is exactly one
                # such arc since below/above are the tie's immediate neighbours).
                sweep = θabove - θbelow
                sweep > pi && (sweep -= 2pi); sweep < -pi && (sweep += 2pi)
                # candidate departure angle from refnode's own end, for a given
                # sign of (px,py) and a given tied-slot rank (1-based, ascending).
                gap = 12.0
                function depart_at(signpx, signpy, rank, ntied)
                    off = (rank - (ntied + 1) / 2)
                    mx = (kx1 + kx2) / 2 + signpx * off * gap * 2
                    my = (ky1 + ky2) / 2 + signpy * off * gap * 2
                    return atan(my - ry, mx - rx)
                end
                ntied = length(tiedslots)
                function monotonic_score(signpx, signpy)
                    # how well do the tied departures interpolate, in ascending
                    # rank, along the below→above sweep direction? Score = sum of
                    # forward progress (positive is good) minus any backward step.
                    prev = θbelow
                    score = 0.0
                    for rank in 1:ntied
                        θ = depart_at(signpx, signpy, rank, ntied)
                        d = θ - prev
                        d > pi && (d -= 2pi); d < -pi && (d += 2pi)
                        score += sweep >= 0 ? d : -d
                        prev = θ
                    end
                    return score
                end
                if monotonic_score(-px, -py) > monotonic_score(px, py)
                    px, py = -px, -py
                end
            end
        end
        for (t, i) in enumerate(idxs)
            e = drawable[i]
            (x1, y1) = pos(e.a); (x2, y2) = pos(e.b)
            col = get(_RGB, e.colour, "#666")
            # EDGE REPULSION: every node that is NOT
            # an endpoint of this edge counts as an obstacle — the curve avoids
            # them instead of running straight over them (typically over a dot).
            obst = _tutte_obstacles(g, e, nodexy, X, Y, edge_repel)
            # If an endpoint is a :braid node, steer the curve's control point at
            # THAT endpoint along its assigned slot angle
            # (`θangles`, see `_braid_leg_angles`) instead of the raw direction
            # to the target — this is what makes the DRAWN departure angle
            # equal the ASSIGNED slot angle even when two slots share a target
            # (their raw directions coincide exactly, so nothing else could
            # separate them). The endpoint itself STAYS at the node centre
            # (x1,y1)/(x2,y2) unchanged: "this segment departs node v" is read off
            # the exact coincidence of the path's `M`/end point with v's drawn
            # centre, so moving the anchor to the rim would make the leg invisible.
            # This takes precedence over the px,py fan axis (below) for any bundle
            # touching a braid node; bundles with no braid endpoint at all (two
            # trivalent nodes joined twice, say) use the fan heuristic.
            a_dir = e.a isa NodePort && g.nodes[e.a.node].kind === :braid ?
                        braid_legdir(e.a.node, e.a.slot) : nothing
            b_dir = e.b isa NodePort && g.nodes[e.b.node].kind === :braid ?
                        braid_legdir(e.b.node, e.b.slot) : nothing
            if a_dir === nothing && b_dir === nothing
                if nbundle == 1
                    # straight edge: first check whether a foreign node is in the
                    # way — if so, the line becomes a shallow Q arc instead.
                    (qx, qy, moved) = _avoid_quad(x1, y1, (x1+x2)/2, (y1+y2)/2, x2, y2, obst)
                    if moved
                        print(io, """<path d="M $x1 $y1 Q $qx $qy $x2 $y2" fill="none" stroke="$col" stroke-width="3.4"/>""")
                    else
                        print(io, """<line x1="$x1" y1="$y1" x2="$x2" y2="$y2" stroke="$col" stroke-width="3.4"/>""")
                    end
                else
                    off = (t - (nbundle + 1) / 2)             # …,-1,0,1,… (0 = straight middle)
                    gap = 12.0                          # px between adjacent parallel strands
                    mx = (x1 + x2) / 2 + px * off * gap * 2
                    my = (y1 + y2) / 2 + py * off * gap * 2
                    (mx, my, _) = _avoid_quad(x1, y1, mx, my, x2, y2, obst)
                    print(io, """<path d="M $x1 $y1 Q $mx $my $x2 $y2" fill="none" stroke="$col" stroke-width="3.4"/>""")
                end
            else
                # control point for the (single) quadratic control point that
                # `tutte_svg` draws everywhere else. A quadratic has only ONE
                # control point, so it cannot in general hit two independently
                # prescribed tangent directions at both ends EXACTLY — a
                # single point achieving angle α from x1 AND angle β from x2
                # only exists when the two forward RAYS (not full lines)
                # genuinely cross; for a back-to-back braid-braid pair with an
                # "S-shaped" tangent configuration the forward rays can
                # diverge — the LINES still meet, but only behind one of the
                # two endpoints (negative ray parameter), so exact ray
                # intersection there silently discards one end's prescribed
                # angle entirely, which can put a leg at the wrong CYCLIC
                # RANK at that end, not just an imprecise angle.
                #
                # Only the CYCLIC RANK among a node's own legs has to be
                # right, so exactness at BOTH ends is not required.
                # Solve the two rays' intersection via Cramer's rule (`sparam`
                # = signed distance from x1 along a_dir, `tparam` = signed
                # distance from x2 along b_dir):
                #   * both forward (sparam ≥ 0 AND tparam ≥ 0): genuine forward
                #     crossing — use it, both angles exact.
                #   * exactly one negative: that endpoint's exactness is
                #     unreachable from this ray pair — walk from THAT endpoint
                #     along ITS OWN prescribed direction instead (short
                #     distance `k`), which reproduces exactness at the
                #     endpoint that would otherwise have failed and preserves
                #     the correct RANK (not necessarily the exact angle) at
                #     the other endpoint too.
                #   * both negative (denom ~ 0, near-parallel prescribed
                #     directions): degenerate fallback, walk from x1.
                secL = hypot(x2 - x1, y2 - y1)
                k = 0.3 * secL
                if a_dir !== nothing && b_dir !== nothing
                    denom = a_dir[1] * (-b_dir[2]) - a_dir[2] * (-b_dir[1])
                    if abs(denom) > 1e-9
                        sparam = ((x2 - x1) * (-b_dir[2]) - (y2 - y1) * (-b_dir[1])) / denom
                        tparam = ((x2 - x1) * (-a_dir[2]) - (y2 - y1) * (-a_dir[1])) / denom
                    else
                        sparam = tparam = -1.0   # force the degenerate fallback below
                    end
                    if sparam >= 0 && tparam >= 0
                        mx = x1 + sparam * a_dir[1]; my = y1 + sparam * a_dir[2]
                    elseif sparam < 0 && tparam >= 0
                        mx = x2 + k * b_dir[1]; my = y2 + k * b_dir[2]
                    elseif tparam < 0 && sparam >= 0
                        mx = x1 + k * a_dir[1]; my = y1 + k * a_dir[2]
                    else
                        mx = x1 + k * a_dir[1]; my = y1 + k * a_dir[2]
                    end
                elseif a_dir !== nothing
                    mx = x1 + k * a_dir[1]; my = y1 + k * a_dir[2]
                else
                    mx = x2 + k * b_dir[1]; my = y2 + k * b_dir[2]
                end
                # Repulsion only within tight limits: here the control point
                # realizes the prescribed slot angles at the braid node, and only
                # their cyclic order has to be right
                # (`check_braid_alternation`) — `_avoid_quad_limited` therefore
                # allows at most ~26 degrees of rotation per end.
                (mx, my, _) = _avoid_quad_limited(x1, y1, mx, my, x2, y2, obst)
                print(io, """<path d="M $x1 $y1 Q $mx $my $x2 $y2" fill="none" stroke="$col" stroke-width="3.4"/>""")
            end
        end
    end
    # free monochrome circles (Edge(c, Circle(c), Circle(c))) — these are what makes a
    # term ZERO (R7), so each is drawn as a small solid ring; several free circles
    # (of possibly different colours) stack near the disc centre so they don't
    # overlap the nodes/leaves.
    freecircles = [e.colour for e in g.edges if e.a isa Circle && e.b isa Circle]
    for (j, c) in enumerate(freecircles)
        col = get(_RGB, c, "#666")
        fx = cx + (j - (length(freecircles) + 1) / 2) * 34.0
        fy = cy
        print(io, """<circle cx="$fx" cy="$fy" r="14" fill="none" stroke="$col" stroke-width="3.4"/>""")
    end
    # internal nodes
    for (v, nd) in enumerate(g.nodes)
        (x, y) = (X(nodexy[v]), Y(nodexy[v]))
        if nd.kind === :dot
            # SOLID, not hollow: a white fill with a thick ring reads as a full
            # disc when small and as a hole when large, so large dots would look
            # hollow. The rectangle renderer (render/Rect.jl) draws them solid too.
            col = _node_rgb(nd.colours)
            # `class="dot"` marks the dot for any post-processing of the SVG; the
            # fill colour is the node's own colour, not white.
            print(io, """<circle cx="$x" cy="$y" r="6" fill="$col" stroke="$col" stroke-width="3" class="dot"/>""")
        else
            print(io, """<circle cx="$x" cy="$y" r="$(nd.kind === :braid ? 8 : 7)" fill="$(_node_rgb(nd.colours))"/>""")
        end
    end
    # slot numbers (debug aid, as for CircularGraph — render/CircularTutte.jl): write the
    # slot number 1..degree just outside the node boundary onto the respective leg,
    # so the DRAWING says which slot is wired where (wiring errors in hand-built
    # examples can otherwise only be guessed from non-planar crossings, not seen
    # directly). Braid nodes use the already computed `θangles`
    # (smoothly fitted); dot/trivalent (degree 1/3, no ambiguity) take the raw
    # direction to the respective far end directly, with the same degree scaling as
    # below in `_node_radius`/the radius literal above (`6` dot, `7` trivalent).
    if slot_numbers
        for (v, nd) in enumerate(g.nodes)
            nd.kind === :dot && continue                  # dots have only one arm
            (x, y) = (X(nodexy[v]), Y(nodexy[v]))
            rad = nd.kind === :braid ? 8 : 7
            # slot → far end: needed for non-braid nodes to determine the leg
            # angles at all (braid nodes have the finished `θangles` for that).
            slotfar = Dict{Int,Port}()
            for (_, own, far) in _edges_at_node(g, v)
                own isa NodePort || continue
                slotfar[own.slot] = far
            end
            if nd.kind === :braid
                d = 2 * nd.m
                angles = θangles[v]
            else
                d = 3
                posof(p::Leaf) = leafxy[p.k]
                posof(p::NodePort) = nodexy[p.node]
                posof(::Circle) = nothing
                (vx, vy) = nodexy[v]                       # unit-disc coords (not pixel X/Y)
                angles = Vector{Float64}(undef, d)
                for s in 1:d
                    far = get(slotfar, s, nothing)
                    p = far === nothing ? nothing : posof(far)
                    if p === nothing
                        angles[s] = 2pi * (s - 1) / d      # unwired/self-circle fallback
                    else
                        (fx, fy) = p
                        angles[s] = atan(fy - vy, fx - vx)
                    end
                end
            end
            for s in 1:d
                (dx, dy) = (cos(angles[s]), sin(angles[s]))
                lx = x + (rad + 11) * dx
                ly = y + (rad + 11) * dy
                print(io, """<circle cx="$lx" cy="$ly" r="7" fill="#ffffff" fill-opacity="0.85" stroke="none"/>""")
                print(io, """<text x="$lx" y="$(ly + 3.5)" text-anchor="middle" font-size="10" fill="#334">$s</text>""")
            end
        end
    end
    # boundary leaves (labelled discs, same style as the circular renderer)
    #
    # LOOSE LEAF: `leaf_colour` returns `0` when NO edge
    # ends at that boundary position. Colour 0 does not exist (1, 2, 3 are allowed)
    # — the leaf is loose, the diagram broken. An inconspicuous grey "0" would hide
    # the error, so a loose leaf is drawn RED AND DASHED with a "✗", and
    # ALWAYS, not only in debug mode: a correct diagram has no loose leaf, so the
    # picture of a correct diagram does not change — and a broken one should not
    # look pretty in silence. `check_wiring` reports the same as `:leaf_degree`.
    # CENTRING THE LEAF DIGIT: `dominant-baseline="central"` is in
    # the SVG standard, but LIBRSVG IGNORES IT — and librsvg is the way into the PDF
    # (Rsvg + Cairo in the render scripts). The digit therefore sat with its
    # BASELINE on the circle centre, i.e. visibly too high (Liam's finding on the
    # title figure). `dy="0.35em"` is relative to the font size, works in browsers
    # AND librsvg, and survives scaling the `font-size` up afterwards (`boost`).
    for k in 1:length(g.word)
        c = leaf_colour(g, k)
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
    # leaf NUMBERS (debug aid): the index k from the word, radially just OUTSIDE
    # the boundary circle — the same construction as the circular version in
    # render/CircularTutte.jl. Placed outside the disc: inside it already
    # shows the COLOUR, and both in a 10 px circle would overlap.
    if leaf_numbers
        for k in 1:length(g.word)
            (x, y) = (X(leafxy[k]), Y(leafxy[k]))
            (dx, dy) = (x - cx, y - cy)
            L = hypot(dx, dy); (dx, dy) = L > 1e-9 ? (dx / L, dy / L) : (0.0, -1.0)
            lx = x + 15 * dx; ly = y + 15 * dy
            print(io, """<text x="$lx" y="$(ly + 3.5)" text-anchor="middle" font-size="10" fill="#334">$k</text>""")
        end
    end
    # SEAM marker (debug aid): the transition from the LAST to
    # the FIRST leaf, i.e. gap 0. The boundary runs counterclockwise in word order —
    # where that ring closes cannot otherwise be seen on a circle, and without this
    # point "forwards along the leaf ring" cannot be checked on the picture.
    # Drawn through the same `_marker_pt`/`MarkSpec` path as the morphism markers
    # (gap 0 = `MarkSpec` `0`), but with its own look (red ring instead of a filled
    # dot) so it cannot be confused with the black/grey bottom/top marker of a
    # morphism.
    if leaf_numbers && n >= 1
        seam = _marker_pt(0)
        seam !== nothing &&
            print(io, """<circle cx="$(seam[1])" cy="$(seam[2])" r="5.5" fill="none" stroke="#c62828" stroke-width="2"/>""")
    end
    # optional small black dot marking a side of the boundary: the GAP at
    # `mark_between`, which marks the LEFT segment of a morphism (between the start
    # of bottom and the start of top). Computed directly from the gap angle
    # under the SAME rot/reflect as the leaves — not re-derived from
    # `leafxy`, so it also works in the degenerate cases where the two flanking
    # leaves coincide or the boundary is empty (see `_marker_pt`).
    (pt = _marker_pt(mark_between)) !== nothing &&
        print(io, """<circle cx="$(pt[1])" cy="$(pt[2])" r="6.5" fill="#000"/>""")
    # second marker (grey): the RIGHT side of the boundary — between the end of
    # bottom and the end of top, same construction as the
    # left/black marker above.
    (pt2 = _marker_pt(mark_between2)) !== nothing &&
        print(io, """<circle cx="$(pt2[1])" cy="$(pt2[2])" r="6.5" fill="#999"/>""")
    print(io, "</svg>")
    return String(take!(io))
end
