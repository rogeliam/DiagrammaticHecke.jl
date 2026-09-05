# render/tutte/Avoid.jl — two small geometry helpers used by both Tutte
# renderers (render/tutte/Tutte.jl and render/CircularTutteSVG.jl).
#
# What they are for:
#   1. A dot must not visually merge with its node or leaf — the gap is exactly ONE
#      dot diameter between the dot's edge and the node's/leaf's edge
#      (`_dot_centre_distance`).
#   2. Node repulsion alone leaves edges running straight through foreign nodes and
#      dots. `_avoid_quad`/`_avoid_cubic` push the control points of a Bezier curve
#      sideways until the curve clears every foreign node.
#
# All sizes are in PIXELS (drawing coordinates), matching the fixed radii `6`
# (dot), `10` (leaf disc), and `_circular_node_radius`.

const _DOT_R  = 6.0     # drawn radius of a dot
const _LEAF_R = 10.0    # drawn radius of a boundary-leaf disc
# In debug mode the node number enlarges the disc to `8 + 2.5·digits` (see
# `node_numbers` in CircularTutteSVG.jl) — for the dot gap we conservatively
# use the larger value, otherwise the gap would again be too small in the
# default (=debug) image.
const _NODE_LABEL_R = 10.5

"""
    _dot_centre_distance(parent_r, gap = 2 * _DOT_R) -> Float64

Centre-to-centre distance of a dot from its neighbour: the neighbour's edge
plus `gap` clearance plus the dot's own radius. The default `gap = 2·_DOT_R` is
exactly **one dot diameter** — so another dot would still fit between the dot
and the node/leaf. Adjustable via `dot_gap` in
`circular_tutte_positions`/`tutte_svg`/`circular_tutte_svg`.
"""
_dot_centre_distance(parent_r::Real, gap::Real = 2 * _DOT_R) = parent_r + gap + _DOT_R

"""
    _avoid_shift(sample, obstacles; clearance = _DOT_R) -> (dx, dy)

Core of the edge repulsion. `sample(t)` gives the curve point at `t ∈ [0,1]`,
`obstacles` is a list `(x, y, radius)` of nodes that are NOT an endpoint of the
edge. Returns the push the curve would need at its midpoint so that it keeps
`radius + clearance` clearance around every node (0,0 if nothing is violated).
"""
function _avoid_shift(sample::Function, obstacles; clearance::Real = _DOT_R)
    sx = sy = 0.0
    isempty(obstacles) && return (0.0, 0.0)
    for (ox, oy, orad) in obstacles
        req = orad + clearance
        best = Inf; bx = by = 0.0
        for t in 0.05:0.05:0.95
            (px, py) = sample(t)
            d = hypot(px - ox, py - oy)
            if d < best
                best = d; bx = px; by = py
            end
        end
        best >= req && continue
        (dx, dy) = (bx - ox, by - oy)
        L = hypot(dx, dy)
        if L < 1e-6                      # curve passes through the centre
            (dx, dy) = (0.0, 1.0); L = 1.0
        end
        need = req - best
        sx += need * dx / L; sy += need * dy / L
    end
    return (sx, sy)
end

"""
    _avoid_quad(x1, y1, mx, my, x2, y2, obstacles; clearance = _DOT_R) -> (mx, my, moved)

Pushes the control point `(mx,my)` of a quadratic Bezier curve away from the
`obstacles` (iteratively, by at most `0.8·chord`). `moved` says whether any push
happened at all — an edge that would otherwise be drawn as a straight line must
then be emitted as a `Q` path instead of a `line`.
"""
function _avoid_quad(x1, y1, mx, my, x2, y2, obstacles; clearance::Real = _DOT_R)
    cap = 0.8 * hypot(x2 - x1, y2 - y1)
    tx = ty = 0.0                                   # accumulated shift
    for _ in 1:12
        cx = mx + tx; cy = my + ty
        b(t) = ((1-t)^2 * x1 + 2*(1-t)*t*cx + t^2 * x2,
                (1-t)^2 * y1 + 2*(1-t)*t*cy + t^2 * y2)
        (sx, sy) = _avoid_shift(b, obstacles; clearance = clearance)
        hypot(sx, sy) < 0.25 && break
        # The curve's apex moves only HALF as far as the control point
        # displacement — hence factor 2.
        tx += 2 * sx; ty += 2 * sy
        L = hypot(tx, ty)
        L > cap && (tx *= cap / L; ty *= cap / L)
    end
    return (mx + tx, my + ty, hypot(tx, ty) > 1e-9)
end

"""
    _avoid_cubic(x1, y1, c1x, c1y, c2x, c2y, x2, y2, obstacles; clearance = _DOT_R)
        -> (c1x, c1y, c2x, c2y, moved)

Like `_avoid_quad`, but for the cubic Beziers of the circular renderer. BOTH
control points are shifted by the same vector, so that the prescribed departure
directions at the nodes (the slot angles) are preserved.
"""
function _avoid_cubic(x1, y1, c1x, c1y, c2x, c2y, x2, y2, obstacles;
                      clearance::Real = _DOT_R)
    cap = 0.8 * hypot(x2 - x1, y2 - y1)
    tx = ty = 0.0
    for _ in 1:12
        a1x = c1x + tx; a1y = c1y + ty; a2x = c2x + tx; a2y = c2y + ty
        b(t) = ((1-t)^3*x1 + 3*(1-t)^2*t*a1x + 3*(1-t)*t^2*a2x + t^3*x2,
                (1-t)^3*y1 + 3*(1-t)^2*t*a1y + 3*(1-t)*t^2*a2y + t^3*y2)
        (sx, sy) = _avoid_shift(b, obstacles; clearance = clearance)
        hypot(sx, sy) < 0.25 && break
        # Apex of a cubic with equally-shifted control points: 3/4 of the
        # displacement, hence factor 4/3.
        tx += (4/3) * sx; ty += (4/3) * sy
        L = hypot(tx, ty)
        L > cap && (tx *= cap / L; ty *= cap / L)
    end
    return (c1x + tx, c1y + ty, c2x + tx, c2y + ty, hypot(tx, ty) > 1e-9)
end

"""
    _tutte_obstacles(g, e, nodexy, X, Y, on) -> Vector{Tuple{Float64,Float64,Float64}}

Obstacle list for edge `e`: all nodes of `g` except its own endpoints, with
their DRAWN radius (dot `6`, braid `8`, otherwise `7`; in debug mode the node
additionally carries its number, hence `_NODE_LABEL_R` as a lower bound for
non-dots). `X`/`Y` convert unit-circle coordinates to pixel coordinates.
`on = false` gives the empty list, so `edge_repel = false` draws every edge
straight.
"""
function _tutte_obstacles(g, e, nodexy, X, Y, on::Bool)
    out = Tuple{Float64,Float64,Float64}[]
    on || return out
    skip = Set{Int}()
    e.a isa NodePort && push!(skip, e.a.node)
    e.b isa NodePort && push!(skip, e.b.node)
    for (v, nd) in enumerate(g.nodes)
        v in skip && continue
        rad = nd.kind === :dot ? _DOT_R :
              max(nd.kind === :braid ? 8.0 : 7.0, _NODE_LABEL_R)
        push!(out, (X(nodexy[v]), Y(nodexy[v]), rad))
    end
    return out
end

"""
    _circular_obstacles(g, e, nodexy, X, Y, on) -> Vector{Tuple{Float64,Float64,Float64}}

Same for the circular renderer; the radius comes from `_circular_node_radius`
(resp. `_DOT_R` for single-arm nodes).
"""
function _circular_obstacles(g, e, nodexy, X, Y, on::Bool)
    out = Tuple{Float64,Float64,Float64}[]
    on || return out
    skip = Set{Int}()
    e.a isa NodePort && push!(skip, e.a.node)
    e.b isa NodePort && push!(skip, e.b.node)
    for v in 1:length(g.nodes)
        v in skip && continue
        nd = g.nodes[v]
        rad = arm_count(nd) == 1 ? _DOT_R : max(_circular_node_radius(nd), _NODE_LABEL_R)
        push!(out, (X(nodexy[v]), Y(nodexy[v]), rad))
    end
    return out
end

"""
    _avoid_quad_limited(x1, y1, mx, my, x2, y2, obstacles;
                        clearance = _DOT_R, max_angle = 0.20) -> (mx, my, moved)

Like `_avoid_quad`, but the control point may rotate the DEPARTURE ANGLES at
both ends by at most `max_angle` (radians, default ~26°). This is the variant
used for legs at `:braid` nodes (render/tutte/Tutte.jl): there the control
point realizes the assigned slot angles, and only the cyclic ORDER of the legs
must remain correct — a small rotation is harmless, a large one could flip the
order.
"""
function _avoid_quad_limited(x1, y1, mx, my, x2, y2, obstacles;
                             clearance::Real = _DOT_R, max_angle::Real = 0.45)
    (nx, ny, moved) = _avoid_quad(x1, y1, mx, my, x2, y2, obstacles; clearance = clearance)
    moved || return (mx, my, false)
    ang(ax, ay, bx, by) = atan(by - ay, bx - ax)
    a1 = ang(x1, y1, mx, my); a2 = ang(x2, y2, mx, my)
    wrap(d) = atan(sin(d), cos(d))
    λ = 1.0
    for _ in 1:12                      # bisection on the push distance
        cx = mx + λ * (nx - mx); cy = my + λ * (ny - my)
        d1 = abs(wrap(ang(x1, y1, cx, cy) - a1))
        d2 = abs(wrap(ang(x2, y2, cx, cy) - a2))
        max(d1, d2) <= max_angle && return (cx, cy, λ > 1e-3)
        λ /= 2
    end
    return (mx, my, false)
end
