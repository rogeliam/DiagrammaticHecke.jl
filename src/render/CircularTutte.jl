# render/CircularTutte.jl — planar Tutte-style drawing for CircularGraph (open-arm nodes).
#
# Kept as a separate file, purely additive on top of render/tutte/Tutte.jl: the
# WordGraph/Tutte.jl stack is untouched. Phase 1 (centre positions, reused from
# `tutte_positions`, Tutte.jl:269-306) is node-type agnostic; only
# `_node_neighbours` needed a `CircularGraph` counterpart.
#
# Unlike the plain renderer, a circular node has no fixed rotation system: several
# legs of the same node can point at the same target node, giving identical raw
# angles with nothing to disambiguate. `_circular_leg_angles` (Phase 2) solves this
# with a real ROTATION SYSTEM per node plus a cubic Bézier per leg (instead of a
# straight line): a cubic has two tangent degrees of freedom, so both slot
# angles (start AND target node) are hit exactly, so parallel edges land in
# different slots and get different angles automatically via the monotonicity
# repair.
#
# GUARANTEE: `_circular_leg_angles` only guarantees the LOCAL cyclic order at each
# node (slot order is respected when drawing) — NOT global crossing-freeness of
# the whole diagram. That is exactly the property the model requires
# (diagram/Graph.jl §PLANAR SLOT CONVENTION).
#
# Reused unchanged from the plain stack (render/tutte/*, render/Morphism.jl,
# render/Graph.jl): `_RGB`, `_RGBMIX`, `_rgbmix`, `_node_rgb`, `_ring_transform`,
# `_mark_angle`, the `_marker_pt` pattern, `_svg_to_png`, `_empty_diagram_svg`,
# `MarkSpec`.
#
# ⚠ SPLIT: the SVG EMITTER and the display glue
# (`circular_tutte_svg`, `Base.show`/`showable`, `display_circular_tutte`, `show_circular_step`)
# now live in render/CircularTutteSVG.jl, included directly after this file. What stays
# here is the LAYOUT: leaf ring, barycentre solver, departure angles, node radius.

# ---- 1. Leaf colour: `_circular_leaf_colour` lives in circular/CircularGraph.jl (also needed by
#      circular/rules/CircularBraidRules.jl §C8), not defined here.

# ---- Debug mode for display_circular_tutte ------------------------------------------

const _DEBUG_CIRCULAR = Ref(false)

"""
    debug_circular(on::Bool) -> Bool
    debug_circular() -> Bool

Global switch for `display_circular_tutte`: when `true`, arm numbers
(`slot_numbers`), node numbers (`node_numbers`) and cell numbers
(`cell_numbers`) are shown by default. Example notebooks set it once at the top.
"""
debug_circular() = _DEBUG_CIRCULAR[]
debug_circular(on::Bool) = (_DEBUG_CIRCULAR[] = on)

# ---- 2. Node neighbours -------------------------------------------------------
#
# `_node_neighbours` (render/tutte/Tutte.jl) reads only `g.edges` and the
# shared port types (`Leaf`/`NodePort`) — nothing `Node`-specific — and is
# typed `Union{WordGraph,CircularGraph}` there, so it is reused directly here.
# Consumed below by the barycentre solver under the same name.

# ---- 3. Slot → target port ----

"""
    _circular_slot_far(g::CircularGraph, v::Int) -> Dict{Int,Port}

For each wired slot of node `v`: the PORT at the other end of the edge. A
`CircularGraph`-specific counterpart to `rules/Helpers.jl`'s `_edges_at_node`
(WordGraph-typed).
"""
function _circular_slot_far(g::CircularGraph, v::Int)
    slotfar = Dict{Int,Port}()
    for e in g.edges
        # Not elseif: a self-connecting edge (both ends at v, e.g. two slots of a
        # trivalent node wired to each other) must register BOTH directions, or
        # the second slot's far port is silently missing — this is why self-loops
        # never rendered otherwise: _circular_self_loop_slots requires both slotfar[s]
        # and slotfar[t] to be set.
        if e.a isa NodePort && e.a.node == v
            slotfar[e.a.slot] = e.b
        end
        if e.b isa NodePort && e.b.node == v
            slotfar[e.b.slot] = e.a
        end
    end
    return slotfar
end

# Circular counterpart of `wiring_table` (render/tutte/Wiring.jl): a CircularNode has no
# `kind` and no fixed arm count; the arm colour sits directly in the arm sequence.
_wiring_rows(g::CircularGraph, v::Int, nd::CircularNode) =
    (_circular_slot_far(g, v), arm_count(nd), string("(", arms(nd), ")"))

_wiring_slot_colour(nd::CircularNode, s::Int) = arm_colour(nd, s)

"""
    _circular_self_loop_slots(g::CircularGraph, v::Int) -> Union{Nothing,Tuple{Int,Int}}

Two slots `(s, t)` of node `v` if they are directly connected to each other
(a self-connection, e.g. arising from `merge_at_edge` on a bigon edge — see
`_fr_needle`). `nothing` if `v` has no such loop.
"""
function _circular_self_loop_slots(g::CircularGraph, v::Int)
    slotfar = _circular_slot_far(g, v)
    for (s, far) in slotfar
        if far isa NodePort && far.node == v
            t = far.slot
            far2 = get(slotfar, t, nothing)
            if far2 isa NodePort && far2.node == v && far2.slot == s
                return (s, t)
            end
        end
    end
    return nothing
end

# ---- 4. Positions: leaf ring + barycentre solver (adapted from Tutte.jl:224-306) ----

"""
    circular_tutte_positions(g::CircularGraph; iters = 400, repel = 0.006,
                        mark_between = nothing, mark_between2 = nothing,
                        bottom_leaves = nothing, top_leaves = nothing)
        -> (leafxy, nodexy)

Coordinates for `g` in the unit disc: boundary leaves fixed on the circle in
word order, internal nodes placed by the barycentre rule plus short-range
repulsion. This is the shared solver for both diagram types — WordGraph
rendering goes through `circular(g)` + `circular_tutte_positions`. Phase 2
(`_circular_leg_angles`) determines leg angles for all node degrees via a
rotation system per node.
"""
function circular_tutte_positions(g::CircularGraph; iters::Int = 400, repel::Real = 0.02,
                             clip::Real = 0.92,
                             mark_between::MarkSpec = nothing,
                             mark_between2::MarkSpec = nothing,
                             bottom_leaves::Union{Nothing,Vector{Int}} = nothing,
                             top_leaves::Union{Nothing,Vector{Int}} = nothing,
                             dot_gap::Real = 2 * _DOT_R)
    n = length(g.word)
    baseθ(k) = pi/2 + 2pi * (k - 1) / max(n, 1)
    rot, reflect = _ring_transform(baseθ, n, mark_between, bottom_leaves, top_leaves)
    # HALF-CIRCLE assignment, same convention as the plain renderer
    # (render/tutte/Tutte.jl, `_halfcircle_angles`): the seam between the end of
    # bottom and the start of top sits on the left, bottom points down, top
    # points up. This holds even WITH a marking: `display(::CircularGraph)` sets
    # `_circular_default_mark(g) = 0` (the gap between leaf n and leaf 1), which is
    # exactly that seam. The one exception is a marking on a DIFFERENT gap,
    # which is rotated to the left.
    hc = _halfcircle_angles(n, bottom_leaves, top_leaves)
    rotate_only = mark_between !== nothing && !(mark_between isa Int && mark_between == 0) &&
                  (bottom_leaves === nothing || top_leaves === nothing ||
                   isempty(bottom_leaves) || isempty(top_leaves))
    # The angle that is ACTUALLY drawn -- the nudge must be based on this.
    θused(k) = rotate_only ? reflect * (baseθ(k) + rot) : hc[k]

    # The nudge pushes the two leaves NEXT TO a marking apart. Its direction
    # `sign(Δ)` must come from `θused` (the angles actually drawn, which run
    # DESCENDING under `_halfcircle_angles`), not from `baseθ` (ASCENDING) —
    # using `baseθ` flips the sign and pulls the leaves TOGETHER instead.
    nudge = Dict{Int,Float64}()
    function _add_nudge!(mb)
        fl = _mark_flank_leaves(n, mb)
        fl === nothing && return
        a, b = fl
        θa, θb = θused(a), θused(b)
        Δ = θb - θa
        Δ > pi && (Δ -= 2pi); Δ < -pi && (Δ += 2pi)
        # Cap at 0.6·π/n: a larger nudge would let the two seam leaves jump past
        # their marking and swap places.
        δ = min(0.12, 0.6 * pi / max(n, 1)) * sign(Δ)
        nudge[a] = get(nudge, a, 0.0) - δ
        nudge[b] = get(nudge, b, 0.0) + δ
    end
    _add_nudge!(mark_between)
    _add_nudge!(mark_between2)

    leafxy = Dict{Int,Tuple{Float64,Float64}}()
    for k in 1:n
        θ = θused(k) + get(nudge, k, 0.0)
        leafxy[k] = (cos(θ), sin(θ))
    end
    nn = length(g.nodes)
    nodexy = Dict{Int,Tuple{Float64,Float64}}()
    for v in 1:nn
        φ = 2pi * (v - 1) / max(nn, 1)
        # Initial radius: starts further out, giving the repulsion force room
        # to spread nodes evenly, same as the plain renderer.
        r0 = nn <= 1 ? 0.0 : 0.45
        nodexy[v] = (r0 * cos(φ), r0 * sin(φ))
    end
    nbrs = [_node_neighbours(g, v) for v in 1:nn]
    # Dots (degree 1) are NOT moved by the barycentre solver: a dot has only ONE
    # neighbour, so its "barycentre" is exactly that neighbour's position — the
    # dot would land on top of it instead of visibly next to it. Instead it is
    # positioned below (after the solver) along its actual departure direction.
    # As a neighbour of another node the dot still counts with its (later-set)
    # position — only the dot itself is not moved here.
    is_dot = [arm_count(nd) == 1 for nd in g.nodes]
    for _ in 1:iters
        maxmove = 0.0
        for v in 1:nn
            is_dot[v] && continue
            isempty(nbrs[v]) && continue
            sx = sy = 0.0
            for (kind, key, _) in nbrs[v]
                (x, y) = kind === :leaf ? leafxy[key] : nodexy[key]
                sx += x; sy += y
            end
            bx = sx / length(nbrs[v]); by = sy / length(nbrs[v])
            rx = ry = 0.0
            if repel > 0
                cx0, cy0 = nodexy[v]
                for u in 1:nn
                    u == v && continue
                    dx = cx0 - nodexy[u][1]; dy = cy0 - nodexy[u][2]
                    d2 = dx*dx + dy*dy
                    d2 < 1e-6 && (dx = 1e-3*(u - v); dy = 1e-3; d2 = dx*dx + dy*dy)
                    d2 < 0.8 || continue
                    f = repel / d2
                    rx += f * dx; ry += f * dy
                end
            end
            nx = 0.6 * bx + 0.4 * nodexy[v][1] + rx
            ny = 0.6 * by + 0.4 * nodexy[v][2] + ry
            rr = hypot(nx, ny)
            # `clip` = largest radius a node may take (in leaf-circle units),
            # exposed as a parameter: in very full images many nodes land on
            # exactly this radius and form a ring right in front of the leaf
            # discs. A smaller `clip` gives breathing room without bending edges
            # afterwards.
            rr > clip && (nx *= clip/rr; ny *= clip/rr)
            maxmove = max(maxmove, hypot(nx - nodexy[v][1], ny - nodexy[v][2]))
            nodexy[v] = (nx, ny)
        end
        maxmove < 1e-9 && break
    end

    # Heuristic: if a node is mostly attached to boundary leaves, pull it
    # further towards the centre so it doesn't sit right at the edge next to them.
    for v in 1:nn
        is_dot[v] && continue
        nbr = nbrs[v]
        isempty(nbr) && continue
        leaf_count = 0
        for (kind, _, _) in nbr
            kind == :leaf && (leaf_count += 1)
        end
        frac = leaf_count / length(nbr)
        if frac >= 0.6
            # factor in [0.6, 1.0]: frac=1.0 -> 0.6 (pulled in strongly),
            # frac=0.6 -> ~1.0 (little or no pull).
            factor = 1.6 - frac
            cx, cy = nodexy[v]
            nodexy[v] = (cx * factor, cy * factor)
        end
    end

    # Nodes with a self-connection (loop) need extra room for the loop (see the
    # drawing logic below) — two of their "neighbours" are self-references and
    # don't pull outward in the barycentre solver, so such nodes would otherwise
    # stick right at the leaf boundary. Extra pull to the centre, independent of
    # the leaf_count heuristic above.
    for v in 1:nn
        is_dot[v] && continue
        _circular_self_loop_slots(g, v) === nothing && continue
        cx, cy = nodexy[v]
        nodexy[v] = (cx * 0.7, cy * 0.7)
    end

    # Placing dots: once all non-dot nodes are fixed, each dot is placed along
    # the REAL slot angle of its parent node (not an artificial "away from
    # other neighbours" fan-out).
    #
    # ORDERING: `_circular_leg_angles` needs `nodexy` and is normally called only
    # AFTER `circular_tutte_positions` (in `circular_tutte_svg`). To place dots correctly
    # here, `_circular_leg_angles` is computed a second time in advance, purely for
    # this purpose (no structural change to `_circular_leg_angles` or to
    # `circular_tutte_positions`'s return value). The later call in `circular_tutte_svg`,
    # used for the actual Bézier curves, gives the same angles up to negligible
    # numerical differences, since dots are never treated as reliable
    # neighbours in the angle computation anyway (see below).
    #
    # MUTUAL DEPENDENCY, RESOLVED: `_circular_leg_angles` derives a neighbour's angle
    # from its position (`raw[s]`), including dot neighbours — but a dot is not
    # yet placed at this point. This is well-defined because `nodexy[v]` for a
    # not-yet-placed dot `v` still holds its barycentre-solver initial value
    # (dots are skipped there, `is_dot[v] && continue`), i.e. close to the
    # origin. A dot target thus lies extremely close to its own parent node,
    # and the distance-weighting mechanism in `_circular_leg_angles` (`weight(r)`,
    # saturating around 0.3, falling off quadratically for near neighbours)
    # already treats the dot neighbour as unreliable and gives it near-zero
    # weight — so the not-yet-placed dot does not meaningfully rotate its
    # parent.
    #
    # DOT DISTANCE: between a dot's edge and the
    # edge of its neighbour there should be room for exactly ONE more dot. The
    # distance therefore depends on the NEIGHBOUR (leaf disc 10 px, node disc
    # `_circular_node_radius`, enlarged to `_NODE_LABEL_R` in debug mode by the
    # node number) and is computed in pixels by `_dot_centre_distance`
    # (render/tutte/Avoid.jl). `_DOT_UNIT` converts pixels into the unit-disc
    # units this function works in: the leaf-circle radius is `0.34·600` by
    # default. A fixed `0.12` would be ~24 px, less than one dot diameter of
    # clearance.
    # If the dot is the ONLY node (the generator pictures 1-dot/2-dot/3-dot),
    # there is nothing to collide with, and it may sit far towards the centre.
    _DOT_UNIT = 1 / (0.34 * 600)
    dotdist_leaf = length(g.nodes) == 1 ? 0.72 :
                   _dot_centre_distance(_LEAF_R, dot_gap) * _DOT_UNIT
    θangles_pre = _circular_leg_angles(g, leafxy, nodexy)
    for v in 1:nn
        is_dot[v] || continue
        nbr = nbrs[v]
        isempty(nbr) && continue
        (kind, key, _) = only(nbr)
        if kind === :leaf
            (nx, ny) = leafxy[key]
            r = hypot(nx, ny)
            # Direction from the leaf TOWARDS the disc centre (not away from
            # it!): the leaf already sits on the unit circle, so the dot must be
            # placed INWARD, or it would land outside the drawing disc.
            (dx, dy) = r > 1e-9 ? (-nx / r, -ny / r) : (1.0, 0.0)
            nodexy[v] = (nx + dotdist_leaf * dx, ny + dotdist_leaf * dy)
        else
            (nx, ny) = nodexy[key]
            # Find dot `v`'s slot at the parent node `key`, to read its actual
            # departure angle from the pre-computed rotation.
            parent_slotfar = _circular_slot_far(g, key)
            slot = only(s for (s, far) in parent_slotfar
                        if far isa NodePort && far.node == v)
            θ = θangles_pre[key][slot]
            (dx, dy) = (cos(θ), sin(θ))
            prad = max(_circular_node_radius(g.nodes[key]), _NODE_LABEL_R)
            dd = length(g.nodes) == 1 ? 0.72 : _dot_centre_distance(prad, dot_gap) * _DOT_UNIT
            nodexy[v] = (nx + dd * dx, ny + dd * dy)
        end
    end

    return leafxy, nodexy
end

# ---- 5. Phase 2: departure angles θ(v,s), the core piece ----------------------

"""
    _circular_leg_angles(g::CircularGraph, leafxy, nodexy) -> Dict{Int,Vector{Float64}}

For each internal node `v` of degree `d`: a vector `θ` of length `d` with the
departure angle of each slot `1..d`, constructed so the angles run STRICTLY
MONOTONICALLY around the full circle in slot order (the planar cyclic
convention, diagram/Graph.jl §PLANAR SLOT CONVENTION) and sum to exactly `2π`.

Algorithm (four steps):
1. `raw[s]` = direction from `nodexy[v]` to the centre of the target port wired
   at slot `s`; unwired or target essentially coincident (`hypot < 1e-9`) ⇒ `NaN`.
2. Best uniform phase `φ` (least-squares fit to the even distribution
   `unif[s] = φ + 2π(s-1)/d`) over all non-`NaN` slots — the rotation system,
   depending only on degree, never on node kind.
3. Blend between the raw angle and the uniform distribution with `blend = 0.55`
   (`NaN` slots take `unif[s]` directly).
4. Monotonicity repair: unfold forward, clamp each of the `d` cyclic gaps to at
   least `0.35·(2π/d)`, then renormalize so the sum is exactly `2π`.

This GUARANTEES ONLY the LOCAL cyclic order at `v` (step 4 forces strictly
ascending angles in slot order) — NOT global crossing-freeness of the whole
diagram (two different nodes can still have crossing legs if their centres
happen to lie that way). That is exactly the property the model actually needs
(diagram/Graph.jl §PLANAR SLOT CONVENTION: slot order IS the cyclic order,
nothing more is required).
"""
function _circular_leg_angles(g::CircularGraph, leafxy, nodexy)
    nn = length(g.nodes)
    result = Dict{Int,Vector{Float64}}()
    posof(p::Leaf) = leafxy[p.k]
    posof(p::NodePort) = nodexy[p.node]
    for v in 1:nn
        d = arm_count(g.nodes[v])
        (vx, vy) = nodexy[v]
        slotfar = _circular_slot_far(g, v)
        raw = fill(NaN, d)
        # Distance to the target port per slot — a very near neighbour (e.g. a
        # dot sitting almost on the circular node itself) gives a numerically
        # UNRELIABLE direction (small position shifts swing the angle a lot) and
        # should influence the node's rotation less than a distant, stable
        # boundary-leaf neighbour.
        dist = fill(NaN, d)
        for s in 1:d
            far = get(slotfar, s, nothing)
            (far === nothing || far isa Circle) && continue
            (fx, fy) = posof(far)
            dx = fx - vx; dy = fy - vy
            r = hypot(dx, dy)
            r < 1e-9 && continue
            raw[s] = atan(dy, dx)
            dist[s] = r
        end
        # Weight per slot: saturates around 0.3 distance units (boundary leaves
        # are typically ~1.0 from the centre), falls off quadratically to near 0
        # for very near neighbours.
        weight(r) = isnan(r) ? 0.0 : clamp(r / 0.3, 0.0, 1.0)^2
        # Step 2: best uniform phase φ over the known slots, WEIGHTED by distance.
        sinsum = cossum = 0.0
        any_known = false
        for s in 1:d
            isnan(raw[s]) && continue
            w = weight(dist[s])
            w <= 0 && continue
            any_known = true
            u = 2pi * (s - 1) / d
            sinsum += w * sin(raw[s] - u)
            cossum += w * cos(raw[s] - u)
        end
        φ = any_known ? atan(sinsum, cossum) : 0.0
        unif = [φ + 2pi * (s - 1) / d for s in 1:d]
        # Step 3: blend.
        # reduce blending towards the uniform distribution so the raw directions
        # from neighbouring ports dominate more — this reduces visible rotation
        # of the node's slot system relative to the leaves. Per-slot blend is
        # additionally weakened (towards the uniform angle) for very near
        # neighbours, since their raw direction is unreliable.
        #
        # blend0 itself grows with the degree: at low degree (3) the raw
        # neighbour directions should still dominate (e.g. a trivalent node
        # orienting itself towards its actual targets), but at high degree
        # (6, e.g. a braid node) unevenly-spaced neighbours (a near dot-arm
        # pulling the barycentre off-centre) otherwise distort the legs far from
        # their even 60°-apart rest positions — so lean much more towards
        # uniform spacing there.
        blend0 = clamp(0.20 + 0.15 * (d - 3), 0.20, 0.80)
        θ = Vector{Float64}(undef, d)
        for s in 1:d
            if isnan(raw[s])
                θ[s] = unif[s]
            else
                w = weight(dist[s])
                blend = 1.0 - w * (1.0 - blend0)          # w=1 -> blend0; w=0 -> 1 (fully uniform)
                Δ = raw[s] - unif[s]
                Δ = mod(Δ + pi, 2pi) - pi                 # wrap into (-π, π]
                θ[s] = unif[s] + (1 - blend) * Δ
            end
        end
        # Step 4: monotonicity repair — unfold forward, clamp gaps, then
        # renormalize to sum = 2π. Use a smaller minsector to allow tighter
        # spacing when raw angles suggest it.
        gaps = Vector{Float64}(undef, d)
        for s in 1:d
            s2 = mod1(s + 1, d)
            Δ = θ[s2] - θ[s]
            (s == d) && (Δ += 2pi)                        # the last jump closes the circle
            Δ = mod(Δ, 2pi)
            Δ == 0 && (Δ = 2pi)                            # edge case: exactly 0 ⇒ full revolution
            gaps[s] = Δ
        end
        minsector = 0.25 * (2pi / d)
        gaps = max.(gaps, minsector)
        gaps .*= 2pi / sum(gaps)                          # renormalize: sum exactly 2π
        θfixed = Vector{Float64}(undef, d)
        θfixed[1] = θ[1]
        for s in 2:d
            θfixed[s] = θfixed[s - 1] + gaps[s - 1]
        end
        result[v] = θfixed
    end
    return result
end

# ---- 6. Radius and fill of a circular node ------------------------------------------

"""
    _circular_node_radius(nd::CircularNode) -> Float64

Draw radius for `nd`: arm count 1 (a dot) → `6`; otherwise growing with arm
count, `clamp(6.0 + 0.55·(deg-2), 7, 13)` — arm count 6 (the plain m=3 braid) →
`8.2`, matching the braid radius `r=8` in render/tutte/Tutte.jl.
"""
function _circular_node_radius(nd::CircularNode)
    d = arm_count(nd)
    d == 1 && return 6.0
    return clamp(6.0 + 0.55 * (d - 2), 7.0, 13.0)
end

"""
    _circular_node_fill(nd::CircularNode) -> String

Colour rule: `_rgb_weighted(arms(nd))` — the arm-count-weighted mean of the arm
colours in RGB space; for `n` arms of colour 1 and `m` of colour 3 that's
`(n·colour1 + m·colour3)/(n+m)`. Monochrome nodes are unaffected (the mean of
equal values is the value itself). `_rgb_weighted` makes no assumption about
the number of distinct colours, unlike the earlier fixed-mix `_node_rgb`
approach it replaced.
"""
_circular_node_fill(nd::CircularNode) = _rgb_weighted(arms(nd))

# ---- Cell anchors for the cell-number overlay ----------------------------------
#
# `_circular_cell_anchors` and `_circular_cell_labels_svg` live in render/CircularCells.jl, not
# here; `_push_anchors_off_nodes!` is the shared helper from render/Cells.jl.
