# render/Rect.jl  —  draw a MorphismGraph as a RECTANGLE (bottom → top), not a circle.
#
# The circular renderer (render/Graph.jl) draws the boundary on a disc and hides the
# bottom/top cuts. For a morphism the conventional string-diagram picture is a
# rectangle: the DOMAIN word on the bottom edge (left→right), the CODOMAIN word on the
# top edge (left→right), internal vertices stacked in between, time flowing upward.
# A morphism with a top and a bottom is printed with its top and bottom as top and
# bottom.
#
# Layout is SIMPLE (not a planarity engine):
#   * bottom leaves evenly on the bottom edge, in `_bottom_leaves` order;
#   * top leaves evenly on the top edge, in `top`-word order (reverse of _top_leaves);
#   * each internal vertex gets a HEIGHT from its edge-distance to the bottom leaves
#     (so stacked braids rise in order) and an X from the mean of what it touches;
#   * edges are drawn as gentle vertical Béziers.
# Edge crossings are NOT avoided in general — for the braid-join diagrams we build
# (identities through, one braid vertex per layer) it reads cleanly; complex diagrams
# may cross. Documented limitation.
#
# reuses _RGB / _node_rgb from render/Morphism.jl.


"""
    morphism_rect_svg(m::MorphismGraph; size = 320, label = "", lw = 1.0, marks = 1.0,
                      frame = true, rails = false) -> String

Draw `m` as a rectangle: bottom word on the bottom edge, top word on the top edge,
internal vertices stacked in between (time upward). Simple layout — see file header.

* `lw` scales every stroke width, `marks` scales every marker radius and font size.
  Both default to `1.0` (the on-screen size). Raise them for a beamer slide, where
  the default hairlines vanish: `morphism_rect_svg(m; lw = 2, marks = 1.6)`.
  Scaling `size` alone does NOT help — the coordinates grow, the strokes do not.
* `frame = false` drops the dashed rectangle and the `bot`/`top` labels (for a
  slide, where the frame is noise).
* `rails = true` draws just two dashed HORIZONTAL lines, on the bottom and the top
  edge, and nothing else — the domain/codomain boundary without the box or the
  labels. The BOTTOM rail is black, the TOP one grey, so the two can be told apart.
  Independent of `frame`; on a slide use `frame = false, rails = true`.
"""
function morphism_rect_svg(m::MorphismGraph; size::Int = 320, label::AbstractString = "",
                           lw::Real = 1.0, marks::Real = 1.0, frame::Bool = true,
                           rails::Bool = false, dotgap::Real = 0.22,
                           dotstair::Real = 0.0, fan::Real = 20.0,
                           caprise::Real = 1.0, dotscale::Real = 1.0)
    g = m.graph
    bl = _bottom_leaves(m)                 # bottom leaf indices, left→right
    tl = reverse(_top_leaves(m))           # top leaf indices so the TOP WORD reads L→R
    botset = Set(bl); topset = Set(tl)

    _w(x) = round(x * lw, digits = 2)      # stroke width
    _r(x) = round(x * marks, digits = 2)   # radius / font size

    W = size; H = size
    padx = 34.0; pady = 30.0
    ybot = H - pady; ytop = pady
    # COMMON GRID for both edges (with 3 leaves below and 2 above,
    # spreading each row over the full width alone dragged the top letters out to the
    # far corners). Both rows use the same column spacing, set by the wider row, and
    # the narrower row is CENTRED in it — so a strand that just goes straight up
    # stays roughly vertical.
    ncol = max(length(bl), length(tl))
    xat(i, n) = ncol <= 1 ? W / 2 :
        padx + (W - 2padx) * ((i - 1) + (ncol - n) / 2) / (ncol - 1)

    lxb = Dict(k => xat(i, length(bl)) for (i, k) in enumerate(bl))
    lxt = Dict(k => xat(i, length(tl)) for (i, k) in enumerate(tl))

    nn = length(g.nodes)
    # HEIGHT: one distinct layer per node, in node-INDEX order. For a graph built by
    # `braid_top`/`Zamo` (one vertex appended per step, bottom→top), the node index IS
    # the build order = the vertical layer, so no two vertices share a height. A BFS
    # heuristic would put several nodes on one level and let them overlap.
    vy = Dict{Int,Float64}()
    for ni in 1:nn
        frac = nn == 1 ? 0.5 : (ni) / (nn + 1)       # in (0,1), evenly spaced, off the edges
        vy[ni] = ybot + (ytop - ybot) * frac
    end
    # A :dot TERMINATES its strand, so it does not need a layer of its own: drawn at
    # the generic layer height its stub runs far up the picture and crosses whatever
    # else is there (on the light leaves of 121 — the 2-dot cut
    # across the 1-strand). Pull every dot that hangs off a BOTTOM leaf down close to
    # that leaf, and every dot hanging off a TOP leaf up close to it; a dot between
    # two inner nodes keeps its computed height.
    # `dotgap` controls HOW close: at 0.22 the dot sat right on the leaf marker and
    # hid it. Larger = farther from the rail.
    # `dotstair`: a STAIRCASE instead of equal height. The light-leaf construction
    # places the dots STEP BY STEP; with `dotstair > 0` the first dot sits lowest,
    # each further one `dotstair` higher — the picture shows its own build order.
    # The rank is the node index, which IS the build order (see HEIGHT above).
    # `dotstair = 0` (default): all dots sit at equal height.
    botdots = [ni for (ni, nd) in enumerate(g.nodes) if nd.kind === :dot &&
               any(p_ isa NodePort && p_.node == ni && q_ isa Leaf && haskey(lxb, q_.k)
                   for e in g.edges for (p_, q_) in ((e.a, e.b), (e.b, e.a)))]
    topdots = [ni for (ni, nd) in enumerate(g.nodes) if nd.kind === :dot &&
               any(p_ isa NodePort && p_.node == ni && q_ isa Leaf && haskey(lxt, q_.k)
                   for e in g.edges for (p_, q_) in ((e.a, e.b), (e.b, e.a)))]
    brank = Dict(ni => i for (i, ni) in enumerate(botdots))
    trank = Dict(ni => i for (i, ni) in enumerate(topdots))
    for (ni, nd) in enumerate(g.nodes)
        nd.kind === :dot || continue
        for e in g.edges, (p_, q_) in ((e.a, e.b), (e.b, e.a))
            p_ isa NodePort && p_.node == ni || continue
            if q_ isa Leaf
                if haskey(lxb, q_.k)
                    gap = dotgap + dotstair * (get(brank, ni, 1) - 1)
                    vy[ni] = ybot + gap * (ytop - ybot)
                end
                if haskey(lxt, q_.k)
                    gap = dotgap + dotstair * (get(trank, ni, 1) - 1)
                    vy[ni] = ytop + gap * (ybot - ytop)
                end
            end
        end
    end
    # X: leaf-touching nodes sit at the mean x of their leaves; nodes that touch no leaf
    # get their x by PROPAGATION from already-placed neighbours (so inner braids don't
    # all pile onto W/2). Iterate to a fixpoint.
    vx = Dict{Int,Float64}()
    leafx = Dict{Int,Float64}()
    for ni in 1:nn
        xs = Float64[]
        for e in g.edges, (p, q) in ((e.a, e.b), (e.b, e.a))
            p isa NodePort && p.node == ni || continue
            q isa Leaf || continue
            haskey(lxb, q.k) && push!(xs, lxb[q.k])
            haskey(lxt, q.k) && push!(xs, lxt[q.k])
        end
        isempty(xs) || (leafx[ni] = sum(xs) / length(xs))
    end
    for ni in 1:nn; vx[ni] = get(leafx, ni, W / 2); end
    for _ in 1:(nn + 2)                              # propagate x from neighbours
        for ni in 1:nn
            haskey(leafx, ni) && continue             # fixed by its leaves
            xs = Float64[]
            for e in g.edges, (p, q) in ((e.a, e.b), (e.b, e.a))
                p isa NodePort && p.node == ni || continue
                q isa NodePort && push!(xs, vx[q.node])
                q isa Leaf && (haskey(lxb, q.k) ? push!(xs, lxb[q.k]) :
                               haskey(lxt, q.k) && push!(xs, lxt[q.k]))
            end
            isempty(xs) || (vx[ni] = sum(xs) / length(xs))
        end
    end

    port_xy(p::Leaf) = haskey(lxb, p.k) ? (lxb[p.k], ybot) : (lxt[p.k], ytop)
    port_xy(p::NodePort) = (vx[p.node], vy[p.node])
    port_xy(::Circle) = (W / 2, H / 2)

    io = IOBuffer()
    disp = 300
    exH = H + 24
    print(io, """<svg width="$disp" height="$(round(Int, disp*exH/W))" """ *
              """viewBox="0 0 $W $exH" xmlns="http://www.w3.org/2000/svg" font-family="monospace">""")
    # the rectangle frame (dashed): bottom and top edges emphasised
    if frame
        print(io, """<rect x="$padx" y="$ytop" width="$(W-2padx)" height="$(ybot-ytop)" fill="none" stroke="#e2e7ee" stroke-dasharray="3" stroke-width="1"/>""")
        print(io, """<text x="$(padx-13)" y="$(ybot+4)" font-size="$(_r(10))" fill="#9aa4b2" text-anchor="end">bot</text>""")
        print(io, """<text x="$(padx-13)" y="$(ytop+4)" font-size="$(_r(10))" fill="#9aa4b2" text-anchor="end">top</text>""")
    end
    # just the two boundary rails: a dashed line below the domain
    # word and one above the codomain word, no box and no labels. Scales with `lw`
    # so it stays visible next to thick strands on a beamer slide.
    if rails
        rx1 = padx - 6; rx2 = W - padx + 6
        # BOTTOM = black, TOP = grey. The two rails are what
        # disappears when the rectangle is closed up into a circle, so they must be
        # TELLABLE APART: on the "planar graphs" slide the point is exactly that the
        # domain and the codomain stop being distinguished.
        for (y, col) in ((ybot, "#1a1a1a"), (ytop, "#9aa4b2"))
            print(io, """<line x1="$rx1" y1="$y" x2="$rx2" y2="$y" stroke="$col" """ *
                      """stroke-dasharray="$(_w(4)) $(_w(3))" stroke-width="$(_w(1.4))"/>""")
        end
    end

    # free circles (R7 markers) — small loops in the middle band
    circ = [e for e in g.edges if e.a isa Circle || e.b isa Circle]
    for (ci, e) in enumerate(circ)
        col = get(_RGB, e.colour, "#666")
        ox = W/2 + (length(circ) == 1 ? 0.0 : 40*(ci - (length(circ)+1)/2)); oy = H/2
        print(io, """<circle cx="$ox" cy="$oy" r="13" fill="none" stroke="$col" stroke-width="2.4"/>""")
    end

    # edges (under vertices): vertical-ish cubic Béziers. Several edges between the SAME
    # two vertices would overdraw as one curve (port_xy maps to the node, not the slot),
    # so we FAN a bundle: spread the Bézier control points by a symmetric horizontal
    # offset (the middle strand stays straight).
    drawable = [e for e in g.edges if !(e.a isa Circle || e.b isa Circle)]
    endkey(p::Leaf)     = (:leaf, p.k)
    endkey(p::NodePort) = (:node, p.node)
    bundle = Dict{Tuple,Vector{Int}}()
    for (i, e) in enumerate(drawable)
        k = Tuple(sort([endkey(e.a), endkey(e.b)]))
        push!(get!(bundle, k, Int[]), i)
    end
    for (_, idxs) in bundle
        nb = length(idxs)
        for (t, i) in enumerate(idxs)
            e = drawable[i]
            col = get(_RGB, e.colour, "#666")
            (x1, y1) = port_xy(e.a); (x2, y2) = port_xy(e.b)
            # `fan`: how far apart a BUNDLE of parallel edges between the same two
            # vertices is spread. The three strands between the two braid vertices of
            # `r5_c`, for example, sit almost on top of each other at a fixed 20px.
            off = nb == 1 ? 0.0 : (t - (nb + 1) / 2) * float(fan)   # px, 0 = straight middle
            if e.a isa NodePort && e.b isa NodePort && e.a.node == e.b.node
                # self-loop at a vertex: small bump to the side
                print(io, """<path d="M$x1,$y1 C$(x1+26),$(y1-22) $(x1+26),$(y1+22) $x2,$y2" stroke="$col" stroke-width="$(_w(2.4))" fill="none"/>""")
            elseif e.a isa Leaf && e.b isa Leaf && y1 == y2
                # CAP / CUP: both ends on the SAME edge of the rectangle. The generic
                # branch below puts both control points at `my == y1`, i.e. ON that edge
                # — the arc came out as a straight line running along the boundary and
                # through whatever leaf sat between the two ends (e.g.
                # `light_leaf([1,2,1], [1,0,1])`, whose cap joins leaves 1 and 3).
                # Arch it into the rectangle instead: bottom leaves bulge UP, top leaves
                # bulge DOWN, by an amount set by how far apart the two feet are.
                inward = (y1 ≈ ybot) ? -1.0 : 1.0
                # Rise: the cap has to clear whatever sits between its two feet. Since
                # dots hang LOW (right next to their leaf), the rise must be large
                # enough to arch over them rather than graze past, while staying
                # below the node band. `caprise`: a multiplier on that arch, for
                # cases where the default still grazes a dot.
                rise = inward * caprise * min(0.52 * abs(x2 - x1) + 16.0,
                                              0.58 * (ybot - ytop))
                print(io, """<path d="M$x1,$y1 C$x1,$(y1 + rise) $x2,$(y2 + rise) $x2,$y2" stroke="$col" stroke-width="$(_w(2.4))" fill="none"/>""")
            else
                my = (y1 + y2) / 2
                print(io, """<path d="M$x1,$y1 C$(x1+off),$my $(x2+off),$my $x2,$y2" stroke="$col" stroke-width="$(_w(2.4))" fill="none"/>""")
            end
        end
    end

    # leaves: same marker as the circular renderer — a white disc with a coloured ring
    # and the letter inside (render/Graph.jl:165-166).
    for k in vcat(bl, tl)
        c = leaf_colour(g, k); col = get(_RGB, c, "#666")
        (x, y) = k in botset ? (lxb[k], ybot) : (lxt[k], ytop)
        print(io, """<circle cx="$x" cy="$y" r="$(_r(9))" fill="white" stroke="$col" stroke-width="$(_w(2))"/>""")
        print(io, """<text x="$x" y="$(y + _r(4))" text-anchor="middle" font-size="$(_r(11))" fill="$col">$c</text>""")
    end

    # vertices. The DOT is filled, not hollow (on a slide the hollow
    # ring reads as a hole in the strand) and noticeably bigger than the other markers —
    # it is the one generator you must be able to spot from the back row.
    for (ni, nd) in enumerate(g.nodes)
        (x, y) = (vx[ni], vy[ni])
        if nd.kind === :dot
            # `dotscale`: without it the dot scales with `marks` just like the leaf
            # discs — at the eight light leaves (marks = 3.0) that made it as big as a
            # leaf and touching the merge node, so it needs its own, smaller scale.
            print(io, """<circle cx="$x" cy="$y" r="$(_r(6.5 * dotscale))" fill="$(_node_rgb(nd.colours))" stroke="$(_node_rgb(nd.colours))" stroke-width="$(_w(2.2 * dotscale))"/>""")
        else
            fill = _node_rgb(nd.colours)
            rad = nd.kind === :braid ? 7.5 : 6.0
            print(io, """<circle cx="$x" cy="$y" r="$(_r(rad))" fill="$fill"/>""")
        end
    end

    if !isempty(label)
        print(io, """<text x="$(W/2)" y="$(exH-7)" font-size="$(_r(11))" fill="#556" text-anchor="middle">$label</text>""")
    end
    print(io, "</svg>")
    return String(take!(io))
end

"PNG of the rectangle rendering (via Rsvg/Cairo, same path as diagram_png)."
morphism_rect_png(m::MorphismGraph; scale::Real = 1) =
    _svg_to_png(morphism_rect_svg(m), scale)

# ---- display wrapper: `rect(m)` renders as a rectangle in a notebook ---------
#
# `display(m)` keeps drawing the circle (the default). Wrap in `rect(m)` to get the
# rectangle view instead: `display(rect(g))`. Both MIME types offered (SVG in
# browser-Jupyter, PNG in VS Code) — same convention as WordGraph.

"Wrapper so `display(rect(m))` draws the MorphismGraph `m` as a rectangle (bottom→top)."
struct Rect
    m::MorphismGraph
end
rect(m::MorphismGraph) = Rect(m)

Base.show(io::IO, ::MIME"image/svg+xml", r::Rect) = print(io, morphism_rect_svg(r.m))
Base.show(io::IO, ::MIME"image/png", r::Rect) = write(io, morphism_rect_png(r.m))
Base.showable(::MIME"image/svg+xml", ::Rect) = true
Base.showable(::MIME"image/png", ::Rect) = true
