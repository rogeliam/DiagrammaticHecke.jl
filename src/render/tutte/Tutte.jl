# render/tutte/Tutte.jl — node-neighbour lookup and the braid-node departure-angle
# rotation system (_braid_leg_angles) used by the Tutte SVG drawer. Included before
# TutteSVG.jl.

#
# An ALTERNATIVE to the concentric shell renderer (`shell_svg`, render/Morphism.jl).
# The shell renderer shows the peel structure but can overlap; Tutte gives a clean
# crossing-free layout for a planar diagram.
#
# Method (Tutte 1963): pin the boundary leaves on a circle in word order (the outer
# face), then place every internal node at the BARYCENTRE of its neighbours. For a
# 3-connected planar graph this is provably crossing-free; for our diagrams it is a
# clean, stable layout even when 3-connectivity isn't guaranteed. Solved by iterating
# the barycentre condition to a fixpoint (no external solver needed).
#
# reuses _RGB / _node_rgb from render/Morphism.jl.

# neighbours of node v: list of (kind, key) where kind ∈ (:leaf,:node), plus the colour.
#
# Takes `Union{WordGraph,CircularGraph}`: the body only ever reads `g.edges` (a
# `Vector{Edge}` on both types) and the shared `Leaf`/`NodePort` port types —
# nothing `Node`-specific — so one implementation serves both.
# `circular_tutte_positions` calls it directly.
function _node_neighbours(g::Union{WordGraph,CircularGraph}, v::Int)
    nb = Tuple{Symbol,Int,Int}[]           # (:leaf|:node, key, colour)
    for e in g.edges
        for (p, q) in ((e.a, e.b), (e.b, e.a))
            p isa NodePort && p.node == v || continue
            if q isa Leaf
                push!(nb, (:leaf, q.k, e.colour))
            elseif q isa NodePort
                push!(nb, (:node, q.node, e.colour))
            end
        end
    end
    return nb
end

# ---- braid-node departure angles -----------------------------------
#
# WHY THIS EXISTS. There is no per-node rotation system: a leg's departure angle is
# just `atan` of the raw direction to its target point. That
# is fine almost everywhere, but at a `:braid` node two DIFFERENT slots can point
# to the SAME neighbour node (necessarily consecutive slots, since a braid
# alternates colours s,t,s,t,… and a repeated neighbour can only recur at
# slot+1) — their raw angles then coincide EXACTLY, and sorting by angle (which
# is what `diagram/Graph.jl`'s §PLANAR SLOT CONVENTION requires the drawing to
# respect: "the slot indices 1,2,…,2m ARE the planar cyclic order") has nothing
# to distinguish them by, so the tie can and does land in the wrong order —
# e.g. a braid node with slots [3,2,3,2,3,2] can draw in angle order
# [3,4,6,5,1,2], breaking alternation at the 6/5 tie. The existing parallel-
# bundle/fan heuristic just above only reorders ALREADY-TIED parallel edges by
# picking a bow side — it cannot help a leg whose raw direction is simply wrong
# relative to its slot when the tie is between two angles that are identical
# rather than merely close.
#
# THE FIX (the same as `_circular_leg_angles` in render/CircularTutte.jl, which solves the
# identical problem for `CircularGraph`'s wider-degree nodes — see that module's
# header for the full derivation). Give every `:braid` node its own genuine
# rotation system: a departure angle per slot that (a) tracks the raw direction
# to each slot's target when that's reliable, (b) falls back to even spacing
# when it isn't (unwired slot, or a target sitting on top of the node), and (c)
# is repaired to be STRICTLY MONOTONIC in slot order and sum to exactly 2π, so
# sorting the drawn legs by angle reproduces the model's slot order by
# construction. Scoped to `:braid` nodes only — `:dot` (1 slot) and
# `:trivalent` (3 slots, generically well separated after the barycentre solve)
# do not show this failure mode and stay on the plain-`atan` path.
#
# GUARANTEE, STATED HONESTLY (same caveat as `_circular_leg_angles`): this fixes the
# LOCAL cyclic order at each braid node (what the model / the test require) —
# NOT global crossing-freedom of the whole diagram.
"""
    _braid_leg_angles(g::WordGraph, leafxy, nodexy) -> Dict{Int,Vector{Float64}}

For every `:braid` node `v` (degree `d = 2m`): a vector `θ` of length `d` with
the departure angle of each slot `1..d`, constructed so the angles are strictly
increasing in slot order and sum to exactly `2π`. Non-`:braid` nodes are absent
from the returned `Dict`. Coordinates are in the same unit-disc system as
`leafxy`/`nodexy` (pre `X`/`Y` scaling); the angle is unchanged by that later
isotropic scaling.

Algorithm (identical shape to `_circular_leg_angles`, CircularTutte.jl):
1. `raw[s]` = direction from `nodexy[v]` to the slot's far port's centre;
   `NaN` if the slot is unwired, its far end is a `Circle`, or the far point is
   numerically on top of `v` (`hypot < 1e-9`).
2. Best uniform phase `φ` (least-squares fit of `unif[s] = φ + 2π(s-1)/d`) over
   the known slots.
3. Blend `raw` toward `unif` (`NaN` slots take `unif` directly).
4. Monotonicity repair: unfold forward, clamp each cyclic gap to a minimum
   sector, renormalise so the gaps sum to exactly `2π`.
"""
function _braid_leg_angles(g::WordGraph, leafxy, nodexy)
    result = Dict{Int,Vector{Float64}}()
    posof(p::Leaf) = leafxy[p.k]
    posof(p::NodePort) = nodexy[p.node]
    for v in 1:length(g.nodes)
        nd = g.nodes[v]
        nd.kind === :braid || continue
        d = 2 * nd.m
        (vx, vy) = nodexy[v]
        at = _edges_at_node(g, v)
        slotfar = Dict{Int,Port}()
        for (_, own, far) in at
            own isa NodePort || continue
            slotfar[own.slot] = far
        end
        raw = fill(NaN, d)
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
        weight(r) = isnan(r) ? 0.0 : clamp(r / 0.3, 0.0, 1.0)^2
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
        gaps = Vector{Float64}(undef, d)
        for s in 1:d
            s2 = mod1(s + 1, d)
            Δ = θ[s2] - θ[s]
            (s == d) && (Δ += 2pi)                        # the last hop closes the circle
            Δ = mod(Δ, 2pi)
            Δ == 0 && (Δ = 2pi)                            # edge case: exactly 0 ⇒ full turn
            gaps[s] = Δ
        end
        minsector = 0.25 * (2pi / d)
        gaps = max.(gaps, minsector)
        gaps .*= 2pi / sum(gaps)                          # renormalise: sum exactly 2π
        θfixed = Vector{Float64}(undef, d)
        θfixed[1] = θ[1]
        for s in 2:d
            θfixed[s] = θfixed[s - 1] + gaps[s - 1]
        end
        result[v] = θfixed
    end
    return result
end
