# RENDER — display glue for the Tutte layout (WordGraph).
# render/tutte/Display.jl  —  unified display_* helpers + MorphismGraph marker helpers.
#
# Contains: display_rect/display_shell/display_tutte, _empty_diagram_svg,
# _left_mark/_right_mark, _tutte_morphism, tutte_svg(::MorphismGraph), and the
# Base.show MIME methods that make Tutte the default display.
#
# ---- unified display_* helpers ---------------------------
#
# Clearer than `display(rect(f))`: call `display_rect(f)` and it renders + shows in one
# step. One per layout. Each shows the SVG (browser-Jupyter) and PNG (VS Code) so both
# front-ends work, matching the WordGraph convention.

# accept a MorphismGraph or a plain WordGraph where a graph layout is asked for
_asgraph(x::WordGraph) = x
_asgraph(x::MorphismGraph) = x.graph

"Show `f` as a bottom→top RECTANGLE (string diagram). `f` is a MorphismGraph."
function display_rect(f::MorphismGraph)
    display("image/svg+xml", morphism_rect_svg(f))
    return nothing
end

"Show `g` as CONCENTRIC SHELLS (peel-based). `g` is a WordGraph or MorphismGraph."
function display_shell(g)
    display("image/svg+xml", shell_svg(_asgraph(g)))
    return nothing
end

"Show `g` with a planar TUTTE embedding. `g` is a WordGraph or MorphismGraph.
A MorphismGraph also gets the black left-side marker."
function display_tutte(g::MorphismGraph)
    display("image/svg+xml", _tutte_morphism(g))
    return nothing
end
function display_tutte(g)
    display("image/svg+xml", tutte_svg(_asgraph(g)))
    return nothing
end

# ---- Tutte is the DEFAULT display ------------------------
#
# `display(g)` draws the planar Tutte layout by default. These definitions
# come after render/Graph.jl / morphism/* so they WIN. The flat, clumped
# single-circle renderer is still reachable by name as `diagram_svg`/`diagram_png`.
#
# For a MorphismGraph we also mark the LEFT side of the boundary with a small black dot:
# the arc between the START OF BOTTOM and the START OF TOP. The bottom word is read
# forward, so its start is the FIRST bottom leaf. The top word is read REVERSED, so the
# first LETTER of the top word is the LAST top leaf. The arc between those two leaves is
# the left edge of the morphism (reading counter-clockwise bottom then top-reversed).
#
"SVG of the EMPTY diagram (no boundary, no nodes): just the dashed outer circle."
function _empty_diagram_svg(size::Int; r::Real = 0.34 * size)
    c = size / 2
    return """<svg viewBox="0 0 $size $size" xmlns="http://www.w3.org/2000/svg" font-family="DejaVu Sans Mono, Consolas, Liberation Mono, Courier New, monospace">
<circle cx="$c" cy="$c" r="$r" fill="none" stroke="#e2e7ee" stroke-dasharray="3" stroke-width="1"/></svg>"""
end

"""
    _left_mark(m::MorphismGraph) -> MarkSpec

The LEFT boundary marker of `m`: the GAP where bottom starts / top ends — this is
exactly the gap `m.cut1` sits at (MorphismGraph.jl's cut convention: the gap
between leaf `cut` and leaf `cut+1`). Markers are GAPS, not leaf pairs, so this
never returns `nothing` for a non-empty boundary.

Degenerate cases (n = number of leaves):
  * `n == 0` (`identity_morphism()`, empty boundary): no gap exists at all — the
    symbolic convention marker `:left` (angle π, fixed pre-rotation).
  * `cut1 == cut2 == 0` (bottom = whole boundary, top = ε) or
    `cut1 == cut2 == _EMPTY_BOTTOM_CUT` (bottom = ε, top = whole boundary): both
    cuts coincide on ONE gap. The single gap is genuinely where bottom "starts"
    (reading the whole/only arc) in both readings, so gap `mod(m.cut1, n)`
    (`_EMPTY_BOTTOM_CUT` normalises to the same single gap via `mod`) is correct
    and well-defined; `_right_mark` uses the SAME gap's other side (angle 0
    instead of π under the `:left`/`:right` convention) so the two dots don't
    coincide — see `_right_mark`.
"""
function _left_mark(m::MorphismGraph)
    n = length(m.graph.word)
    n == 0 && return :left
    return mod(m.cut1, n)
end

"""
    _right_mark(m::MorphismGraph) -> MarkSpec

The RIGHT boundary marker of `m`: the gap where bottom ends / top starts — gap
`m.cut2`. Symmetric counterpart of `_left_mark`, likewise gap-based. `n == 0`: the
symbolic `:right` convention marker (angle 0).

`cut1 == cut2` (either degenerate whole/ε reading, `n > 0`): `_left_mark` and
`_right_mark` land on the SAME gap (there is only one, geometrically) — that gap
alone cannot carry two visually distinct dots. Instead `_right_mark` returns the
ANTIPODAL-of-a-gap spec `(:opposite, gap)` (angle `gap`'s own angle `+ π`, exactly
opposite it in the SAME rotated frame — NOT an independent absolute angle: an
independent `0` would drift away from being opposite gap `0`'s marker once `rot ≠
0` rotates gap `0` off of `π`, which is exactly the `n=2` `cap_morphism(1)`
regression this form fixes). The boundary here is entirely one word (bottom or
top), so the "start" and "end" of that one arc are the SAME point topologically,
but the drawing still needs two distinguishable dots — the gap itself (black,
"here bottom/top starts") and its antipode (grey, "here it ends") mark the two
sides of that single gap. This is a drawing convention for a single-point
degeneracy, not a claim that the gap has two locations.
"""
function _right_mark(m::MorphismGraph)
    n = length(m.graph.word)
    n == 0 && return :right
    m.cut1 == m.cut2 && return (:opposite, mod(m.cut1, n))
    return mod(m.cut2, n)
end

# Bottom/top leaf lists, used ONLY to decide the orientation (mirror) fix —
# kept separate from the markers themselves so a caller that omits them still
# gets the rotate-to-marker behaviour without the orientation fix, as
# render/Cells.jl and render/Decorated.jl do, passing `mark_between` alone.
_tutte_morphism(m::MorphismGraph; kwargs...) =
    tutte_svg(m.graph; mark_between = _left_mark(m),
                       mark_between2 = _right_mark(m),
                       bottom_leaves = _bottom_leaves(m),
                       top_leaves = _top_leaves(m), kwargs...)

"""
    tutte_svg(m::MorphismGraph) -> String

A morphism is ALWAYS drawn as a morphism — identical to `display(m)`.

This method exists so that `tutte_svg(m)` need not be reconstructed by hand as
`tutte_svg(m.graph; mark_between = _left_mark(m))`: such a hand-rolled call is
missing `bottom_leaves`/`top_leaves`, which ORIENT the leaf ring (bottom below,
top above). Without them `_halfcircle_angles` reads the whole word as bottom, and
the image looks twisted even though the wiring itself is still correct.
"""
tutte_svg(m::MorphismGraph; kwargs...) = _tutte_morphism(m; kwargs...)

# THE DEBUG SWITCH ALSO APPLIES TO `display(g)` — the same construction as the
# circular version (render/CircularTutte.jl, `_circular_debug_on`). Without this pass-through
# only an explicit `tutte_svg(g; leaf_numbers = true)` call would consult the
# switch, but `display` would not — and in the notebooks everything goes through
# `display`.
# ⚡ LEAF NUMBERS ALWAYS BELONG IN THE IMAGE.
#
# They do NOT depend on the debug switch (`_plain_debug_on()`): without leaf numbers
# one cannot check whether the arms of a node run backwards around the leaf ring,
# which is the house convention ("arms 2,3,4,5,6 cannot be connected to leaves
# 1,2,3,4,5 — it has to be arms 6,5,4,3,2"). The switch does control the
# CELL/region numbers, which quickly clutter the image.
Base.show(io::IO, ::MIME"image/svg+xml", g::WordGraph) =
    print(io, tutte_svg(g; leaf_numbers = true))
Base.show(io::IO, ::MIME"image/png", g::WordGraph) =
    write(io, _svg_to_png(tutte_svg(g; leaf_numbers = true)))
Base.show(io::IO, ::MIME"image/svg+xml", m::MorphismGraph) =
    print(io, _tutte_morphism(m; leaf_numbers = true))
Base.show(io::IO, ::MIME"image/png", m::MorphismGraph) =
    write(io, _svg_to_png(_tutte_morphism(m; leaf_numbers = true)))
