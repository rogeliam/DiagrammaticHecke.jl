# RENDER — boundary markers (MarkSpec), used by the circular renderer too.
# render/tutte/Markers.jl  —  boundary-marker / angle geometry for the Tutte renderer.
#
# Contains: MarkSpec + the gap-angle helpers (_gap_angle, _nearest_gap, _mark_angle,
# _mark_flank_leaves, _ring_transform, _halfcircle_marker_angle, _halfcircle_angles).
# Used by tutte_positions/tutte_svg (render/tutte/Tutte.jl) and by
# render/tutte/Display.jl, render/CircularTutte.jl. Included BEFORE
# render/tutte/Tutte.jl since MarkSpec is used as a
# type annotation there (needs to exist at parse time).
#
# ---- boundary markers: GAPS, not leaf pairs --------------------
#
# A boundary marker is conceptually a GAP c ∈ 0..n-1 (the gap sits between leaf c
# and leaf c+1, cyclically — MorphismGraph.jl's cut convention) — never a pair of
# leaves. A gap has a well-defined angle: the bisector of its two flanking leaves
# on the SHORT arc; since `baseθ` is linear in k, that bisector is exactly
# `baseθ(c) + π/n` (the point exactly between leaf c and leaf c+1), and this stays
# correct even across the seam (k=n → k=1) because of that linearity.
#
# A LEAF-PAIR tuple `(a, b)` — the two leaves flanking the marked arc, marker at
# their bisector — is accepted as well and normalised to an angle by `_mark_angle`;
# render/Cells.jl and render/Decorated.jl use that form.
#
# Two more spec shapes cover the degenerate cases where the gap alone cannot carry
# two visually distinguishable marker dots:
#   * a bare `Symbol` (`:left`/`:right`) is the CONVENTION marker for `n == 0`
#     (no leaves at all exist, so no gap/pair can be formed): `:left` → π,
#     `:right` → 0, in the SAME rotated frame as everything else (`rot` is
#     always 0 when n == 0 anyway, so this is also the raw angle).
#   * `(:opposite, c)` marks the OTHER side of gap `c` (angle `_gap_angle(c) + π`,
#     i.e. exactly π away from `c`'s own marker in the SAME rotated frame) —
#     used when `cut1 == cut2` (n > 0): both cuts sit on the SAME single gap, so
#     `_right_mark` can't just repeat that gap (the two dots would coincide) but
#     also can't use an independent absolute angle (that would NOT stay π away
#     after rotation, see the n=2 `cap_morphism(1)` regression this fixed).
#     `(:opposite, c)` is guaranteed exactly antipodal to gap `c` under any `rot`.
const MarkSpec = Union{Nothing,Int,Tuple{Int,Int},Symbol,Tuple{Symbol,Int}}

# angle of gap c: the gap sits between leaf c and leaf c+1 (MorphismGraph.jl's cut
# convention — leaves are 1-indexed 1..n, c ranges 0..n-1, c=0 is the gap between
# leaf n and leaf 1). Its angle is the bisector `baseθ(c) + π/n` using leaf c's
# angle where c is taken as a 1-indexed leaf (mod1, so c=0 ↦ leaf n, matching
# "gap 0 is between leaf n and leaf 1").
function _gap_angle(baseθ::Function, n::Int, c::Int)
    n <= 0 && return 0.0
    return baseθ(mod1(c, n)) + pi / n
end

# the gap (0..n-1) whose own angle is closest to `target`. It snaps the
# antipodal-of-a-gap marker (`(:opposite, c)`) into a genuine gap instead of the
# raw `_gap_angle(c) + π` angle, which can land EXACTLY on a leaf (on
# `dot_morphism(1)` the single gap of n==1 is π away from leaf 1, so its raw
# antipode is 2π away, i.e. leaf 1 itself; the same collision occurs for every ODD
# n, since n odd makes `_gap_angle(c) + π` fall exactly on a leaf angle for
# every gap c). Gap angles THEMSELVES never coincide with a leaf angle by
# construction (a gap's bisector always sits strictly between two distinct
# leaves, even for n==1 where that
# "between" is the whole rest of the circle) — so snapping the antipode onto the
# nearest gap guarantees the grey marker always lands in an actual leaf gap.
# For n >= 2 (odd n in particular) this nearest gap is ALWAYS a gap other than
# `c` itself, so black and grey stay visually distinct. n == 1 is special:
# there is only ONE gap in existence, so
# the "nearest gap" to any target is gap 0 itself — the caller handles n == 1
# separately (see `_mark_angle`) instead of calling this helper for it.
function _nearest_gap(baseθ::Function, n::Int, target::Real)
    n <= 0 && return 0
    best_c, best_d = 0, Inf
    for c in 0:(n - 1)
        d = abs(mod(_gap_angle(baseθ, n, c) - target + pi, 2pi) - pi)
        if d < best_d
            best_c, best_d = c, d
        end
    end
    return best_c
end

# normalise a marker spec (gap Int, OLD leaf-pair Tuple, symbolic convention, the
# antipodal-of-a-gap form, or nothing) to an angle, or `nothing` if it doesn't
# resolve (out of range).
function _mark_angle(baseθ::Function, n::Int, mb::MarkSpec)
    mb === nothing && return nothing
    if mb isa Symbol
        return mb === :left ? pi : 0.0
    end
    if mb isa Tuple{Symbol,Int}
        _, c = mb
        n <= 0 && return nothing
        target = _gap_angle(baseθ, n, c) + pi
        if n == 1
            # n == 1: exactly one leaf, one gap — the raw antipode of that gap
            # IS the leaf's own angle (gap and leaf sit on a diameter through
            # each other), so there is no second gap to snap to. Nudge the
            # antipode by a fixed margin (short of a full quarter turn) so the
            # grey marker stays inside the single free gap, near-opposite the
            # black marker but never ON the leaf. `dot_morphism(1)` is the case
            # that needs it.
            return target - pi / 4
        end
        # n >= 2: snap to the nearest ACTUAL gap instead of using the raw
        # antipodal angle, which coincides exactly with a leaf for every odd n
        # (see `_nearest_gap`'s docstring).
        return _gap_angle(baseθ, n, _nearest_gap(baseθ, n, target))
    end
    if mb isa Int
        n <= 0 && return nothing
        return _gap_angle(baseθ, n, mb)
    end
    # OLD form: leaf-pair tuple (a, b) — bisector of the short arc from a to b.
    a, b = mb
    (1 <= a <= n && 1 <= b <= n) || return nothing
    θa, θb = baseθ(a), baseθ(b)
    Δ = θb - θa
    Δ > pi && (Δ -= 2pi); Δ < -pi && (Δ += 2pi)
    return θa + Δ / 2
end

# the two leaves flanking a marker spec, for the cosmetic nudge — nothing if the
# spec doesn't resolve to two distinct real leaves (e.g. a gap on n==0, the
# degenerate single-gap case where "flanking leaves" would coincide, or a
# symbolic convention marker which has no real flanking leaves).
function _mark_flank_leaves(n::Int, mb::MarkSpec)
    (mb === nothing || mb isa Symbol || mb isa Tuple{Symbol,Int}) && return nothing
    if mb isa Int
        (n > 0 && 0 <= mb <= n - 1) || return nothing
        a = mod1(mb, n); b = mod1(mb + 1, n)
        a == b && return nothing
        return (a, b)
    end
    a, b = mb
    (1 <= a <= n && 1 <= b <= n && a != b) || return nothing
    return (a, b)
end

# shared rot/reflect computation (used by both `tutte_positions`, for the actual
# leaf coordinates, and `tutte_svg`, to place the marker dots at the SAME gap
# angles without re-deriving them from already-placed leaf coordinates — which
# breaks for degenerate boundaries where the two flanking leaves of a marker
# coincide or don't exist). n == 0 ⇒ rot = 0, reflect = 1 (nothing to rotate).
#
# WHY A RIGID MOTION IS NOT ENOUGH. `rot`/`reflect` describe a rotation plus a
# mirror, so they preserve the ANGLE between the two boundary gaps, which is fixed
# at `2π·min(|bottom|,|top|)/n` on the short way. That is enough to put bottom's
# mean y below top's — the criterion is the SIGN of a y-difference, and `reflect`
# has full freedom over it — but it CANNOT also put the left marker at x<0 and the
# right marker at x>0 in every case: when the smaller arc spans < 90° of the circle
# (top_frac = |top|/n < 1/4), both gaps lie within 90° of each other, so EVERY
# rotation angle puts them both left or both right of x=0. On
# `tensor(cap_morphism(1), merge_morphism(1))` (n=5, |bottom|=4, |top|=1) none of
# the four combinations of "rotate cut1 or cut2 to π" × "reflect or not" satisfies
# xl<0 ∧ xr>0.
#
# `_halfcircle_angles` therefore PRESCRIBES the leaf positions instead of deriving
# them from a rigid motion.
#
# House convention: the marking between the start of bottom and top sits on the
# left. The bottom word runs ccw W-SW-S-…, and the top word runs cw W-NW-N…; with
# no bottom/top word given, the whole word is bottom and top is empty, i.e. read
# from the bottom upwards.
#
# The three requirements (marker left, bottom below, top above) then hold
# CONSTRUCTIVELY, including the corner case |bottom|=4, |top|=1 where every rigid
# motion puts both markers left.
#
# `reflect` is consequently always 1: a mirrored leaf ring makes the wiring read the
# wrong way round.
function _ring_transform(baseθ::Function, n::Int, mark_between::MarkSpec,
                          bottom_leaves::Union{Nothing,Vector{Int}},
                          top_leaves::Union{Nothing,Vector{Int}})
    rot = 0.0
    θm = _mark_angle(baseθ, n, mark_between)
    θm !== nothing && (rot = pi - θm)
    # `reflect` is ALWAYS 1 — see `_halfcircle_angles` and the block above. It is
    # still returned so that the four call sites (here, `tutte_svg`, and the two
    # counterparts in CircularTutte.jl) can compute `reflect * (…)` uniformly.
    return rot, 1
end

"""
    _halfcircle_marker_angle(hc, n, mb) -> Union{Nothing,Float64}

The angle of a boundary marker under HALF-CIRCLE assignment: the angular midpoint
of the short arc between its two flanking leaves, read off the already assigned
angles `hc`. Under uniform full-circle assignment `_mark_angle` does this via
`baseθ`; with leaves spread over half circles the position derived from `baseθ` does
not match the real ones. `nothing` if the marker has no two flanking leaves.
"""
function _halfcircle_marker_angle(hc::Dict{Int,Float64}, n::Int, mb::MarkSpec)
    fl = _mark_flank_leaves(n, mb)
    fl === nothing && return nothing
    a, b = fl
    (haskey(hc, a) && haskey(hc, b)) || return nothing
    θa, θb = hc[a], hc[b]
    d = mod(θb - θa + pi, 2pi) - pi        # short arc from a to b
    return θa + d / 2
end

"""
    _halfcircle_angles(n, bottom_leaves, top_leaves) -> Dict{Int,Float64}

The leaf angles under the house convention: the marking is on the left, always
fixed; then the bottom word ccw first, then the top word reversed.

* the **marker lies LEFT** (angle π) — fixed, in every case;
* from there runs **ONE continuous ccw ring**: first the **bottom word** in its
  order, then the **top word reversed**;
* all `n` leaves are spread **uniformly** over the FULL circle.

With no bottom/top given, **the whole word is bottom and top is empty** — so it is
always read from the bottom upwards.

Because the ring is continuous the spacing is always equal and the size ratio does
not matter: going from length 5 to 2, the bottom word simply reaches into the
north-east. bottom thus sits **rather below** and top **rather above**, without
forcing the arcs onto half circles.

SVG convention: y grows DOWNWARDS, so "below" in the picture is `sin θ > 0`; ccw in
the picture therefore means **decreasing** θ from π.
"""
function _halfcircle_angles(n::Int, bottom_leaves::Union{Nothing,Vector{Int}},
                            top_leaves::Union{Nothing,Vector{Int}})
    n <= 0 && return Dict{Int,Float64}()
    bl = bottom_leaves === nothing ? collect(1:n) : bottom_leaves
    tl = top_leaves === nothing ? Int[] : top_leaves
    # If the specification is incomplete, read everything as bottom rather than
    # leave leaves without a position.
    (isempty(bl) && isempty(tl)) && (bl = collect(1:n))
    # Ring order from the marker: ONE continuous ccw ring, i.e. bottom LEAVES
    # forwards, then top LEAVES forwards. `top_leaves` is already the ccw boundary
    # arc (`_top_leaves`, cut2+1 … cut1); the top WORD is its reversal. "top
    # reversed" means the WORD reversed, which is exactly the leaf list FORWARDS.
    # (Not `reverse(tl)`: that draws the top word ccw instead of cw, whereas a
    # morphism w→v must read ccw as w·v^{-1}, e.g. id_[2,1] as 2,1,1,2.)
    order = vcat(bl, tl)
    # Append leaves in neither list — otherwise they would have no position (can
    # happen with an incomplete bottom/top specification).
    seen = Set(order)
    for k in 1:n
        k in seen || push!(order, k)
    end
    θ = Dict{Int,Float64}()
    m = length(order)
    # DECREASING from the seam at π: sin θ > 0 is the lower half of the picture
    # (SVG y grows downwards), so θ running π → 0 goes W → SW → S → SE → E, i.e.
    # counter-clockwise in the picture. (Increasing θ from π would run UPWARDS and
    # put the bottom word in the upper half.)
    for (i, k) in enumerate(order)
        θ[k] = pi - 2pi * (i - 0.5) / max(m, 1)
    end
    return θ
end
