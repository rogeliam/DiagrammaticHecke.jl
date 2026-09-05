# RENDER — wiring table / debug labels; also read by render/CircularTutte.jl.
# render/tutte/Wiring.jl  —  debug labelling + the wiring table/console dump.
#
# Contains: the global debug-label switch (debug_labels/_DEBUG_LABELS,
# _plain_debug_on) plus show_wiring/wiring_table and their row helpers
# (_wiring_rows, _wiring_slot_colour, _wiring_edge_colour, _port_label).
#
# ---- Debug mode ---------------------------------
#
# TEMPORARY extra labelling for cross-checking hand-built diagrams by hand — NOT
# the normal display. It affects the CELL labels: instead of just the polynomial,
# each cell shows the TUPLE `(label, distance)`, exactly the pair that
# `find_d4_match` (rules/DecoratedRules.jl) decides on.
#
# A global switch instead of a passed-through keyword, because notebooks set it
# ONCE at the top and every subsequent display (including `display(g)` via the
# Base.show MIME methods, which take no keywords) should see it.
#
# IN THE IMAGE, each leg shows only the ARM NUMBER: the full wiring `s→target`
# there makes the image too cluttered. The wiring is available as a TABLE next to
# the image instead — `wiring_table`, below.

const _DEBUG_LABELS = Ref(false)

"""
    debug_labels(on::Bool) -> Bool
    debug_labels() -> Bool

Toggle debug labelling globally. When on,
`decorated_svg`/`circular_decorated_svg` write the TUPLE `(label, distance)`
into every cell instead of just the polynomial. The gallery/example notebooks
set it to `true` once at the top.

The leg labelling in the image is NOT affected by this — it always shows just
the arm number (`slot_numbers`). Anyone who needs the wiring uses `wiring_table`.
"""
debug_labels() = _DEBUG_LABELS[]
debug_labels(on::Bool) = (_DEBUG_LABELS[] = on)

# THE DEBUG SWITCH FOR THE PLAIN RENDERER, mirroring `_circular_debug_on`
# (render/CircularTutte.jl): both switches apply to both renderers, so a
# notebook doesn't get two different images depending on which one it set at
# the top. What depends on it is documented at `tutte_svg`: LEAF NUMBERS and the
# SEAM marking. The arm numbers (`slot_numbers`) do NOT depend on
# it — they always belong in the image.
# ⚠️ `debug_circular` is only defined in render/CircularTutte.jl, which is included
# AFTER this file. That's fine: the call happens at runtime, by which point the
# module is fully loaded.
_plain_debug_on() = debug_labels() || debug_circular()

"""
    _port_label(p) -> String

Short notation for an edge endpoint, as used when cross-checking the wiring by
hand: `n.s` for `NodePort(n, s)`, `Ls` for `Leaf(s)`,
`∘` for a free circle. Used by `wiring_table`.
"""
_port_label(p::NodePort) = string(p.node, ".", p.slot)
_port_label(p::Leaf)     = string("L", p.k)
_port_label(::Circle)    = "∘"

"""
    show_wiring(g)
    show_wiring(m::MorphismGraph)

The wiring of `g` as a TEXT table on `stdout` — the console counterpart of
`wiring_table`. One row per node per arm, with arm number, arm colour, and the
opposite endpoint.

For a `MorphismGraph` an extra `bottom`/`top` column is added: for a boundary leaf
it records which side of the morphism it lies on.
"""
function show_wiring(g; io::IO = stdout, side = nothing)
    for (v, nd) in enumerate(g.nodes)
        slotfar, d, kindtxt = _wiring_rows(g, v, nd)
        println(io, "Node ", v, " ", kindtxt)
        println(io, "   Arm | Colour | connected to", side === nothing ? "" : " | Side")
        for s in 1:d
            far  = get(slotfar, s, nothing)
            scol = _wiring_slot_colour(nd, s)
            ecol = far === nothing ? nothing : _wiring_edge_colour(g, v, s)
            coltxt = (ecol === nothing || ecol == scol) ? string(scol) :
                     string(scol, "/", ecol, " MISMATCH")
            sidetxt = (side === nothing || !(far isa Leaf)) ? "" :
                      string(" | ", side(far.k))
            println(io, "   ", lpad(s, 3), " | ", lpad(coltxt, 5), " | ",
                    rpad(far === nothing ? "—" : _port_label(far), 6), sidetxt)
        end
    end
    return nothing
end

show_wiring(m::MorphismGraph; io::IO = stdout) =
    show_wiring(m.graph; io = io,
                side = let bl = Set(_bottom_leaves(m)); k -> k in bl ? "bottom" : "top" end)

"""
    wiring_table(g) -> String

The wiring of `g` as an HTML table: one row per node per arm, with ARM NUMBER,
ARM COLOUR, and the OPPOSITE ENDPOINT (`_port_label`: `n.s`, `Lk`, `∘`). Together
with the arm numbers in the image, this lets you read off which arm of which
colour goes where, without tracing curves.

The "colour" column names TWO numbers when they disagree: the SLOT colour (what
the node prescribes for that arm) and the EDGE colour (what the edge actually
attached there carries). If they don't match, the wiring is faulty and the row is
highlighted in red. An unwired slot gets `—`.

Works the same way for `WordGraph` and `CircularGraph` (the arm colour comes from
`_slot_colour` resp. `arm_colour`).
"""
function wiring_table(g)
    io = IOBuffer()
    print(io, """<table style="border-collapse:collapse;font-family:monospace;font-size:13px">""")
    print(io, """<tr><th style="text-align:left;padding:3px 10px;border-bottom:1px solid #ccc">Node</th>""" *
              """<th style="text-align:left;padding:3px 10px;border-bottom:1px solid #ccc">Arm</th>""" *
              """<th style="text-align:left;padding:3px 10px;border-bottom:1px solid #ccc">Colour</th>""" *
              """<th style="text-align:left;padding:3px 10px;border-bottom:1px solid #ccc">connected to</th></tr>""")
    for (v, nd) in enumerate(g.nodes)
        slotfar, d, kindtxt = _wiring_rows(g, v, nd)
        for s in 1:d
            far  = get(slotfar, s, nothing)
            scol = _wiring_slot_colour(nd, s)
            ecol = far === nothing ? nothing : _wiring_edge_colour(g, v, s)
            bad  = ecol !== nothing && scol !== nothing && ecol != scol
            coltxt = (ecol === nothing || ecol == scol) ? string(scol) :
                     string("slot ", scol, " / edge ", ecol)
            bg = bad ? " background:#ffe6e6;" : ""
            print(io, """<tr style="$(bg)">""")
            print(io, """<td style="padding:3px 10px">$(s == 1 ? string(v, " ", kindtxt) : "")</td>""")
            print(io, """<td style="padding:3px 10px">$s</td>""")
            print(io, """<td style="padding:3px 10px">$coltxt$(bad ? " ⚠" : "")</td>""")
            print(io, """<td style="padding:3px 10px">$(far === nothing ? "—" : _port_label(far))</td></tr>""")
        end
    end
    print(io, "</table>")
    return String(take!(io))
end

# Per node: slot→opposite endpoint, arm count, and a short type description.
# Separate implementations for WordGraph (Node with `kind`) and CircularGraph
# (CircularNode, arm sequence only).
function _wiring_rows(g::WordGraph, v::Int, nd::Node)
    slotfar = Dict{Int,Port}()
    for (_, own, far) in _edges_at_node(g, v)
        own isa NodePort && (slotfar[own.slot] = far)
    end
    d = nd.kind === :dot ? 1 : (nd.kind === :braid ? 2 * nd.m : 3)
    return slotfar, d, string("(", nd.kind, " ", nd.colours, ")")
end

_wiring_slot_colour(nd::Node, s::Int) = nd.kind === :braid ? _slot_colour(nd, s) : nd.colours[1]

# The colour of the edge attached to slot `s` of node `v`.
function _wiring_edge_colour(g, v::Int, s::Int)
    for e in g.edges
        for p in (e.a, e.b)
            p isa NodePort && p.node == v && p.slot == s && return e.colour
        end
    end
    return nothing
end
