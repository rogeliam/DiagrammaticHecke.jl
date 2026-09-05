# RENDER — SVG/PNG output for WordGraph (`_svg_to_png` is used by the circular
# renderer too).
# render/Graph.jl  —  draw a WordGraph DIRECTLY as its diagram (not as a path)
#
# A rule is a WordGraph: a boundary word, ONE (or few) internal vertices, and edges.
# The GRAPH itself is drawn — boundary leaves on a circle, the central
# vertex (braid / trivalent) inside with strands from the leaves, dots as degree-1
# stubs, self-loops as arcs joining two slots of the SAME vertex. NOT a multi-ring
# reduction path, which would produce spurious extra rings.
#
# reuse _RGB / _rgbmix / _node_rgb / _empty_svg from render/Morphism.jl

# colour a vertex glyph
_vertex_fill(nd::Node) = _node_rgb(nd.colours)

# incidence: edges at node ni, as (edge_index, this_port, other_port)
function _edges_at(g::WordGraph, ni::Int)
    hits = Tuple{Int,Port,Port}[]
    for (ei, e) in enumerate(g.edges)
        if e.a isa NodePort && e.a.node == ni
            push!(hits, (ei, e.a, e.b))
        elseif e.b isa NodePort && e.b.node == ni
            push!(hits, (ei, e.b, e.a))
        end
    end
    return hits
end

"""
    diagram_svg(g::WordGraph; size = 320, r = 118, label = "") -> String

Draw the WordGraph `g` as its planar diagram: boundary leaves evenly on a circle in
word order (coloured by their letter), each internal vertex placed inside, strands
drawn from leaves/ports to vertices. A `:dot` is a small white-ringed disc at the end
of its single strand; a `:trivalent` / `:braid` is a solid disc where strands meet.
A self-loop (an edge whose two ports are on the SAME vertex) is drawn as an arc.
"""
function diagram_svg(g::WordGraph; size::Int = 320, r::Real = 118,
                     label::AbstractString = "")
    lets = letters(g.word)
    n = length(lets)
    cx = cy = size / 2
    # leaf positions: leaf 1 at top, clockwise
    leafθ(k) = -pi/2 + 2pi * (k - 1) / max(n, 1)
    lx(k) = cx + r * cos(leafθ(k));  ly(k) = cy + r * sin(leafθ(k))

    # vertex positions: a vertex sits at the circular mean of the boundary leaves it
    # touches, pulled inward; a vertex touching no boundary sits at the centre.
    vpos = Dict{Int,Tuple{Float64,Float64}}()
    # vertices touching NO boundary (e.g. barbell dots): spread them on a small circle
    # around the centre so they don't overlap.
    noboundary = [ni for ni in 1:length(g.nodes)
                  if all(!(oth isa Leaf) for (_, _, oth) in _edges_at(g, ni))]
    for (ni, _) in enumerate(g.nodes)
        θs = Float64[]
        for (_, _, oth) in _edges_at(g, ni)
            oth isa Leaf && push!(θs, leafθ(oth.k))
        end
        if isempty(θs)
            j = findfirst(==(ni), noboundary)
            if length(noboundary) == 1
                vpos[ni] = (cx, cy)
            else
                φ = 2pi * (j - 1) / length(noboundary)
                vpos[ni] = (cx + 0.22r * cos(φ), cy + 0.22r * sin(φ))
            end
        else
            θ̄ = _mean_angle(θs)
            depth = r * (0.46 - 0.04 * min(length(θs), 4))
            vpos[ni] = (cx + depth * cos(θ̄), cy + depth * sin(θ̄))
        end
    end
    # a dot vertex should sit a bit BEYOND its (single) neighbour, toward the centre,
    # so its stub strand is visible. Nudge dots inward from their leaf/vertex anchor.
    for (ni, nd) in enumerate(g.nodes)
        nd.kind === :dot || continue
        at = _edges_at(g, ni)
        isempty(at) && continue
        (_, _, oth) = at[1]
        (ox, oy) = oth isa Leaf ? (lx(oth.k), ly(oth.k)) : vpos[oth.node]
        # Place the dot on the segment from its neighbour toward the centre, keeping a
        # CLEAR gap so it never sits on top of a vertex. Push a fixed distance beyond the
        # neighbour (not a fraction of it), so a dot capping a near-boundary trivalent is
        # still pulled well inward.
        dx = cx - ox; dy = cy - oy; L = hypot(dx, dy) + 1e-9
        gap = min(0.42r, 0.55L)                      # how far along toward the centre
        vpos[ni] = (ox + gap * dx / L, oy + gap * dy / L)
    end

    port_xy(p::Leaf) = (lx(p.k), ly(p.k))
    port_xy(p::NodePort) = vpos[p.node]
    port_xy(::Circle) = (cx, cy)                    # unused: circle edges drawn specially

    io = IOBuffer()
    # Explicit width/height (not just viewBox) so the diagram renders at a compact
    # fixed size in a notebook cell instead of stretching to fill the container.
    # `disp` is the on-screen size in px; the coordinate system stays `size`.
    disp = 300
    print(io, """<svg width="$disp" height="$(round(Int, disp*(size+24)/size))" """ *
              """viewBox="0 0 $size $(size+24)" xmlns="http://www.w3.org/2000/svg" font-family="DejaVu Sans Mono, Consolas, Liberation Mono, Courier New, monospace">""")
    print(io, """<circle cx="$cx" cy="$cy" r="$r" fill="none" stroke="#e2e7ee" stroke-dasharray="3" stroke-width="1"/>""")

    # count parallel edges between the same node pair, so we can bow them apart.
    pairkey(e) = e.a isa NodePort && e.b isa NodePort ?
        (min(e.a.node, e.b.node), max(e.a.node, e.b.node)) : nothing
    paircount = Dict{Any,Int}(); pairseen = Dict{Any,Int}()
    for e in g.edges
        k = pairkey(e); k === nothing && continue
        paircount[k] = get(paircount, k, 0) + 1
    end

    # free monochrome circles (Circle-marker edges): draw as small standalone loops
    # near the centre, spread if several. These evaluate to 0 (R7) but we still show
    # them so a reduction step is inspectable.
    circ_edges = [e for e in g.edges if e.a isa Circle || e.b isa Circle]
    for (ci, e) in enumerate(circ_edges)
        col = get(_RGB, e.colour, "#666")
        off = length(circ_edges) == 1 ? 0.0 : 0.24r
        φ = 2pi * (ci - 1) / max(length(circ_edges), 1)
        ox = cx + off * cos(φ); oy = cy + off * sin(φ)
        print(io, """<circle cx="$ox" cy="$oy" r="14" fill="none" stroke="$col" stroke-width="2.4"/>""")
    end

    # edges first (under the vertices)
    for e in g.edges
        (e.a isa Circle || e.b isa Circle) && continue   # handled above
        col = get(_RGB, e.colour, "#666")
        selfloop = e.a isa NodePort && e.b isa NodePort && e.a.node == e.b.node
        (x1, y1) = port_xy(e.a); (x2, y2) = port_xy(e.b)
        if selfloop
            # a loop at one vertex: bow it out toward the centre with a wide C
            (vx, vy) = vpos[e.a.node]
            nx = vx - cx; ny = vy - cy; L = hypot(nx, ny) + 1e-9
            # perpendicular direction to fan the loop
            perpx = -ny / L; perpy = nx / L
            c1x = vx + 26perpx - 22*nx/L; c1y = vy + 26perpy - 22*ny/L
            c2x = vx - 26perpx - 22*nx/L; c2y = vy - 26perpy - 22*ny/L
            print(io, """<path d="M $vx $vy C $c1x $c1y $c2x $c2y $vx $vy" fill="none" stroke="$col" stroke-width="2.4"/>""")
        else
            k = pairkey(e)
            # bow offset: spread parallel edges symmetrically
            off = 0.0
            if k !== nothing && paircount[k] > 1
                idx = get(pairseen, k, 0); pairseen[k] = idx + 1
                off = (idx - (paircount[k]-1)/2) * 26.0
            end
            mx = (x1+x2)/2; my = (y1+y2)/2
            dx = x2-x1; dy = y2-y1; L = hypot(dx,dy)+1e-9
            perpx = -dy/L; perpy = dx/L
            qx = mx + 0.10*(cx-mx) + off*perpx
            qy = my + 0.10*(cy-my) + off*perpy
            print(io, """<path d="M $x1 $y1 Q $qx $qy $x2 $y2" fill="none" stroke="$col" stroke-width="2.4"/>""")
        end
    end

    # boundary leaves: coloured ring + letter. Colour and letter come from the EDGE
    # attached to that leaf (the truth), NOT from `lets[k]` — the boundary word is stored
    # in min-rotation normal form, so `lets[k]` need not match the edge at Leaf(k).
    leafcol = Dict{Int,Int}()
    for e in g.edges
        e.a isa Leaf && (leafcol[e.a.k] = e.colour)
        e.b isa Leaf && (leafcol[e.b.k] = e.colour)
    end
    for k in 1:n
        lc  = get(leafcol, k, lets[k])           # edge colour if wired, else the word letter
        col = get(_RGB, lc, "#666")
        print(io, """<circle cx="$(lx(k))" cy="$(ly(k))" r="9" fill="white" stroke="$col" stroke-width="2"/>""")
        print(io, """<text x="$(lx(k))" y="$(ly(k)+4)" text-anchor="middle" font-size="11" fill="$col">$lc</text>""")
    end

    # vertices
    for (ni, nd) in enumerate(g.nodes)
        (x, y) = vpos[ni]
        if nd.kind === :dot
            col = get(_RGB, nd.colours[1], "#666")
            print(io, """<circle cx="$x" cy="$y" r="5.5" fill="$col" stroke="white" stroke-width="1.5"/>""")
        else
            print(io, """<circle cx="$x" cy="$y" r="7" fill="$(_vertex_fill(nd))"/>""")
        end
    end

    if !isempty(label)
        print(io, """<text x="$cx" y="$(size+16)" font-size="13" fill="#666" text-anchor="middle">$label</text>""")
    end
    print(io, "</svg>")
    return String(take!(io))
end

# (the named example graphs moved to diagram/Fixtures.jl)

"""
    save_diagram(path, g; kwargs...)

Write `diagram_svg(g; kwargs...)` to `path` (creating parent dirs).
"""
function save_diagram(path::AbstractString, g::WordGraph; kwargs...)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, diagram_svg(g; kwargs...))
    end
    return path
end

# ---------------------------------------------------------------------------
# Rich display: let IJulia / Jupyter render a WordGraph as its SVG diagram.
# ---------------------------------------------------------------------------
#
# Jupyter shows a value by querying `show(io, ::MIME, x)` for the MIME types it can
# render, richest first. Declaring `image/svg+xml` makes a bare `g` in a cell draw
# the planar diagram instead of the one-line text summary. The plain-text
# `show(io, g)` (diagram/Graph.jl) stays as the REPL/`print` fallback, and
# `show_wordgraph(g)` is the structural console dump.

# NOTE: the rich-display `show(::MIME, ::WordGraph)` methods live in
# render/tutte/Display.jl — `display(g)` draws the planar Tutte layout. `diagram_svg`
# and `diagram_png` here are the flat single-circle renderer.
Base.showable(::MIME"image/svg+xml", ::WordGraph) = true

# ---- PNG fallback (VS Code notebooks skip image/svg+xml) --------------------
#
# VS Code's notebook renderer prefers image/png and does not reliably display
# image/svg+xml, so we ALSO offer a PNG: rasterise the same `diagram_svg` string
# with Rsvg → Cairo. Jupyter picks the richest MIME each front-end supports, so
# browser-Jupyter still gets the crisp SVG and VS Code gets the PNG.

"""
    diagram_png(g::WordGraph; scale = 1) -> Vector{UInt8}

Rasterise `diagram_svg(g)` to PNG bytes at `scale`× the SVG's nominal size. The
default `scale = 1` gives a compact ~320 px image for a notebook cell (this is the
size VS Code shows, since it renders the PNG, not the SVG); raise it for a sharper
or larger picture, e.g. `diagram_png(g; scale = 2)`. Also handy for
`write("g.png", diagram_png(g))`.
"""
# Rasterise any SVG string to PNG bytes via Rsvg → Cairo (shared by the circular and
# rectangle renderers).
function _svg_to_png(svg::AbstractString, scale::Real = 1)
    handle = Rsvg.handle_new_from_data(String(svg))
    dim = Rsvg.handle_get_dimensions(handle)
    w = ceil(Int, dim.width * scale)
    h = ceil(Int, dim.height * scale)
    surface = Cairo.CairoImageSurface(w, h, Cairo.FORMAT_ARGB32)
    ctx = Cairo.CairoContext(surface)
    Cairo.scale(ctx, scale, scale)
    Rsvg.handle_render_cairo(ctx, handle)
    buf = IOBuffer()
    Cairo.write_to_png(surface, buf)
    return take!(buf)
end

diagram_png(g::WordGraph; scale::Real = 1) = _svg_to_png(diagram_svg(g), scale)

# PNG rich-display method also lives in render/tutte/Display.jl (Tutte default).
Base.showable(::MIME"image/png", ::WordGraph) = true

# ---- running notebook-style code in a plain script ---------------------------

"""
    notebook_display_sink!()

Register a display that ACCEPTS the rich MIMEs a notebook front-end would render
(`image/svg+xml`, `image/png`, `text/html`, `text/latex`), so notebook-style code
can be run by `julia script.jl`.

Without it, any `display("image/svg+xml", …)` — `display_cell_number`,
`display_circular_decorated`, … — dies with

    MethodError: no method matching show(::IOStream, ::MIME"image/svg+xml", ::String)

because a bare script has no MIME-capable display. That is an artefact of script
mode, not a fault in the code being run.

The sink does not merely swallow the call: it asserts the payload is non-empty,
so a renderer that silently produced nothing still fails. That makes it usable as
a check — run a notebook's cells as a script and the drawing code is really
exercised, just not shown.

Inside Jupyter you do not need this: use the project's own kernel, installed with

    using IJulia; installkernel("DiagrammaticHecke", "--project=" * pwd())
"""
function notebook_display_sink!()
    pushdisplay(_NotebookSink())
    return nothing
end

struct _NotebookSink <: AbstractDisplay end

const _SINK_MIMES = ("image/svg+xml", "image/png", "text/html", "text/latex")

for M in _SINK_MIMES
    @eval Base.displayable(::_NotebookSink, ::MIME{Symbol($M)}) = true
end

# A payload that is ALREADY a string (`display("image/svg+xml", svg)`) is taken
# as it is; anything else is RENDERED. Rendering is the point: the sink exists so
# that a script really exercises the drawing code, and `string(x)` on a diagram
# would only produce its text form and assert on that instead.
function Base.display(::_NotebookSink, m::MIME, x)
    m in MIME.(_SINK_MIMES) || throw(MethodError(display, (m, x)))
    payload = x isa AbstractString ? x : repr(m, x)
    @assert !isempty(payload) "empty $(m) payload"
    return nothing
end

# `display(x)` without a MIME picks the first rich form the object supports, the
# same choice a notebook front-end makes. Falling through to `text/plain` here
# would silently skip the renderer.
function Base.display(d::_NotebookSink, x)
    for M in _SINK_MIMES
        m = MIME(M)
        showable(m, x) && return display(d, m, x)
    end
    show(devnull, "text/plain", x)
    return nothing
end
