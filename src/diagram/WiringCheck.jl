# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/WiringCheck.jl — Euler test, check_wiring, audit mode: euler,
# check_wiring/check_wiring_table/show_wiring_check, is_planar_embedding,
# slot_colour_audit and their helpers.
#
# ---- Euler test and convention test -----------------------------------------
#
# Both live here ONCE (the tracer is here anyway), generic over `deg` like the face
# walk itself; `_circular_is_planar_embedding` is a private alias.

"""
    _euler_generic(word, nodes, edges, deg) -> Int

V − E + F of the combinatorial map. `F` is the number of φ cycles of the face walk
(inner faces PLUS the outer face), i.e. `length(t.cycles)`; free circles (Circle
ports) are excluded from `E` because they contribute no cycle in the tracer either —
they add +1 to `E` and +1 to `F` and cancel anyway. Number for number the same as
`V = |nodes|+|word|, E = |edges|+|word|, F = face_count+1`.
"""
function _euler_generic(word::CircularWord, nodes::AbstractVector,
                        edges::Vector{Edge}, deg)
    t = _trace_generic(word, nodes, edges, deg)
    V = length(word) + length(nodes)
    E = count(e -> !(e.a isa Circle || e.b isa Circle), edges) + length(word)
    return V - E + length(t.cycles)
end

"""
    euler(g) -> Int

The Euler characteristic `V − E + F` of the EMBEDDING of `g` (`WordGraph`,
`MorphismGraph`, `CircularGraph`, `CircularMorphismGraph`). It measures the rotation system,
not the graph: `2` means sphere/disc (planar embedding), less means positive genus —
the wiring then contradicts the house convention (orientation block at the top of
diagram/Faces.jl).

The house test; `test/faces.jl` computes it over all fixtures.
"""
euler(g::WordGraph) = _euler_generic(g.word, g.nodes, g.edges, _node_degree)
euler(m::MorphismGraph) = euler(m.graph)

"""
    is_planar_embedding(g) -> Bool

`euler(g) == 2`, i.e. the embedding is planar and the wiring compatible with the
house convention.
"""
is_planar_embedding(g) = euler(g) == 2

# ---- check_wiring: the house convention computed directly -------------------
#
# The convention: the BOUNDARY runs counter-clockwise (leaves in word order), the
# ARMS at a node run clockwise (slot order). That is a statement about the
# ORIENTATION of the rotation system, and orientation is global — the face walk is
# self-consistent under ANY wiring (it pairs darts and follows φ cycles regardless),
# only the genus of the surface it runs on changes. So there is NO complete local
# test; the complete test is `euler(g) == 2`.
#
# What there is are two places where the convention has local content, because a
# SECOND, independent order is available to compare the slot order against — the
# leaf ring:
#
#   (1) LEAF↔NODE. If a node hangs on several leaves, its slots (cw) run BACKWARDS
#       along the leaf ring (ccw). That is the direct translation of "cw at one end,
#       ccw at the other" for node↔boundary connections. The cyclic comparison has
#       content only from THREE leaves on — with two, every cyclic sequence is both
#       forwards and backwards.
#
#   (2) NODE↔NODE. Here the second order is missing; a single node has no external
#       reference. The honest substitute is a MEASUREMENT rather than a formula:
#       reverse the arm sequence of ONE node (slot s ↦ d+1−s, i.e. cw ↦ ccw at that
#       node) and if `euler` rises, that node was attached to its neighbours
#       mirrored — co-oriented instead of counter-oriented at both ends.
#
# An empty return value therefore means: neither of the two locally visible
# violations. Whether the embedding really is planar is what `euler` says; if it is
# not, that appears additionally as the violation `:euler`.

"""
    WiringViolation

A violation of the house convention as reported by `check_wiring`. `kind` ∈
`(:leaf_degree, :slot_colour, :leaf_ring, :node_ring, :reversed_node, :euler)`,
`node` is the node concerned (`0` = global), `detail` the measured reason in plain
text.
"""
struct WiringViolation
    kind::Symbol
    node::Int
    detail::String
end

# Direction of a sequence of PAIRWISE DISTINCT leaf indices on the ring of length n:
# +1 = forwards (ccw, word direction), −1 = backwards, 0 = undecidable (fewer than
# three leaves) or not monotone at all (more than one turn around).
function _ring_direction(ks::Vector{Int}, n::Int)
    m = length(ks)
    (m >= 3 && n >= 3) || return 0
    step(a, b) = mod(b - a, n)
    fwd = sum(step(ks[i], ks[mod1(i + 1, m)]) for i in 1:m)
    fwd == n && return 1
    bwd = sum(step(ks[mod1(i + 1, m)], ks[i]) for i in 1:m)
    bwd == n && return -1
    return 0
end

# slot → opposite end, per node.
function _slot_targets(nodes::AbstractVector, edges::Vector{Edge})
    out = [Dict{Int,Port}() for _ in 1:length(nodes)]
    for e in edges, (p, q) in ((e.a, e.b), (e.b, e.a))
        p isa NodePort && (out[p.node][p.slot] = q)
    end
    return out
end

# The colour that slot `slot` of the node carries by the node's own definition — for
# BOTH worlds, since `_check_wiring_generic` is used by `WordGraph` and `CircularGraph`.
# `nothing` if the slot does not exist or the node type fixes no slot colour. See
# check (0b) in `_check_wiring_generic`.
_slot_colour_generic(nd::Node, slot::Int) =
    1 <= slot <= _node_degree(nd) ? _slot_colour(nd, slot) : nothing
# `CircularNode` (circular/CircularGraph.jl) is included AFTER this file, so the type is not known
# here — hence the dispatch on the field `arms` rather than on the type.
function _slot_colour_generic(nd, slot::Int)
    hasproperty(nd, :arms) || return nothing
    a = getproperty(nd, :arms)
    return 1 <= slot <= length(a) ? a[slot] : nothing
end

# Check (0b) alone — slot colour against edge colour, without the leaf ring and
# without Euler: kept separate from `_check_wiring_generic` because the audit mode
# below runs it on EVERY constructed diagram, where the Euler computation would be
# far too expensive.
function _slot_colour_violations(nodes::AbstractVector, edges::Vector{Edge})
    out = WiringViolation[]
    for e in edges, p in (e.a, e.b)
        p isa NodePort || continue
        1 <= p.node <= length(nodes) || continue
        c = _slot_colour_generic(nodes[p.node], p.slot)
        c === nothing && continue
        c == e.colour && continue
        push!(out, WiringViolation(:slot_colour, p.node,
            "slot $(p.slot) carries colour $c by the node's definition, but the " *
            "edge ending there has colour $(e.colour) — the node's colour order " *
            "does not match the wiring (for `:braid`: `colours` swapped)"))
    end
    return out
end

# ---- audit mode: check EVERY constructed diagram against (0b) ----------------
#
# "Slot colour ≠ edge colour" arises in hand-wired FIXTURES, and those are scattered
# over tests, scripts and notebooks — many not as a literal `Node(...)` but via
# helper constructors. The audit hooks into the 3-argument constructors of
# `WordGraph`/`CircularGraph`, where everything passes through. The call site is
# recorded from the stack so the culprit can be found.
const _SC_AUDIT = Ref(false)
const _SC_AUDIT_LOG = Vector{Tuple{String,Int,String,Vector{WiringViolation}}}()

"""
    slot_colour_audit(on::Bool = true)

Turns on recording of `:slot_colour` violations for ALL `WordGraph`/`CircularGraph`s
constructed from now on. Switching it ON clears the log; read it with
[`slot_colour_audit_log`](@ref).

Costs one pass over the edges per diagram — meant for test runs, not the production
path.
"""
function slot_colour_audit(on::Bool = true)
    # Clear ONLY when switching on. Otherwise switching off at the end of a run
    # erases exactly the log one wants to read afterwards — that mistake can
    # report a whole run as "clean" although a fixture with known violations
    # exists.
    on && empty!(_SC_AUDIT_LOG)
    _SC_AUDIT[] = on
    return on
end

"""
    slot_colour_audit_log() -> Vector{Tuple{String,Int,String,Vector{WiringViolation}}}

The audit log: one entry `(file, line, signature, violations)` for the FIRST call
site outside `src/`, i.e. the fixture or example.

The `signature` (boundary word + node list) is needed because inside a `@testset` the
line points at the testset header, not at the construction — only the signature says
WHICH fixture is meant.
"""
slot_colour_audit_log() = _SC_AUDIT_LOG

# First stack frame outside the package `src/` — that is the fixture or example. In
# addition: if other `src/` code came BEFORE it (other than the constructor and the
# audit itself), a RULE built the diagram and it is an intermediate state, not a
# hand-wired fixture. Without that distinction the rules' intermediate states end up
# filed under the test file that called them, and the report is useless.
const _SC_OWN = ("_sc_audit", "_sc_caller", "WordGraph", "CircularGraph")

function _sc_caller()
    internal = false
    for fr in stacktrace()
        f = String(fr.file)
        if occursin(joinpath("DiagrammaticHecke.jl", "src"), f)
            String(fr.func) in _SC_OWN || (internal = true)
            continue
        end
        occursin(r"^\./|^[.]?(boot|essentials|client)\.jl", f) && continue
        return (f, fr.line, internal)
    end
    return ("?", 0, internal)
end

function _sc_audit(word::CircularWord, nodes::AbstractVector, edges::Vector{Edge})
    _SC_AUDIT[] || return nothing
    vs = _slot_colour_violations(nodes, edges)
    if !isempty(vs)
        f, l, internal = _sc_caller()
        sig = string(internal ? "[intermediate state of a rule] " : "[fixture] ",
                     "word=", join(letters(word)), "  nodes=[",
                     join((sprint(show, nd) for nd in nodes), ", "), "]")
        push!(_SC_AUDIT_LOG, (f, l, sig, vs))
    end
    return nothing
end

# Reverse the arm sequence of node `v`: slot s ↦ d+1−s. Only the PORTS are
# renumbered, `nodes` and `word` stay untouched — this is about the orientation in
# the rotation system only, not about colours.
function _reverse_node_slots(nodes::AbstractVector, edges::Vector{Edge},
                             v::Int, deg)
    d = deg(nodes[v])
    f(p::Port) = p isa NodePort && p.node == v ? NodePort(v, d + 1 - p.slot) : p
    return Edge[Edge(e.colour, f(e.a), f(e.b)) for e in edges]
end

function _check_wiring_generic(word::CircularWord, nodes::AbstractVector,
                               edges::Vector{Edge}, deg)
    out = WiringViolation[]
    n = length(word)
    targets = _slot_targets(nodes, edges)

    # (0) EXACTLY ONE edge ends at EVERY boundary position.
    #
    # `is_wired` (circular/CircularGraph.jl) checks the node slots and colour fidelity but
    # EXPRESSLY not the boundary. A diagram can therefore have boundary positions
    # where no edge ends at all — and nobody notices, because `leaf_colour` returns
    # `0` for want of an edge and the renderer dutifully paints a "0". But colour 0
    # does not exist (1, 2, 3 are allowed); the leaf is LOOSE. Two edges at the same
    # boundary position are the same error from the other side: the boundary is a
    # ring of leaves with ONE strand each.
    leafdeg = zeros(Int, n)
    for e in edges, p in (e.a, e.b)
        p isa Leaf && 1 <= p.k <= n && (leafdeg[p.k] += 1)
    end
    for k in 1:n
        leafdeg[k] == 1 && continue
        push!(out, WiringViolation(:leaf_degree, 0,
            leafdeg[k] == 0 ?
                "leaf $k is LOOSE — no edge ends there (the renderer shows it as " *
                "colour 0, which does not exist)" :
                "$(leafdeg[k]) edges end at leaf $k instead of one"))
    end

    # (0b) SLOT COLOUR = EDGE COLOUR at every node port.
    #
    # A `:braid` node fixes, through the order of its `colours`, which colour the odd
    # and which the even slots carry (`_slot_colour`, diagram/Graph.jl §PLANAR SLOT
    # CONVENTION). If that says `[1,3]` but the edges at the odd slots carry 3, the
    # wiring is self-contradictory — and nothing else reports it: `is_planar_braid`
    # checks only the cyclicity of the slots, `euler` only the embedding, `is_wired`
    # only slot occupancy and edge-colour fidelity port against port.
    #
    # The consequence is severe because the RULES determine their roles via
    # `_slot_colour`: they then grab the wrong connections, and the resulting leaf
    # colours can come out as NO rotation of the boundary word — the word has
    # effectively changed, which no diagram operation may do.
    append!(out, _slot_colour_violations(nodes, edges))

    # (1) leaf↔node: the slots run backwards along the leaf ring.
    for (v, nd) in enumerate(nodes)
        ks = Int[]
        for s in 1:deg(nd)
            q = get(targets[v], s, nothing)
            q isa Leaf && push!(ks, q.k)
        end
        dir = _ring_direction(ks, n)
        dir == 1 && push!(out, WiringViolation(:leaf_ring, v,
            "leaves in slot order $(ks) run FORWARDS along the leaf ring; the " *
            "convention requires backwards (i.e. $(reverse(ks)))"))
        dir == 0 && length(ks) >= 3 && push!(out, WiringViolation(:leaf_ring, v,
            "leaves in slot order $(ks) are not monotone on the leaf ring at all " *
            "— the arms cross"))
    end

    # (1b) node↔node DIRECTLY, where it is locally decidable.
    #
    # For leaves this is check (1), but only from three leaves at the node on. For
    # node↔node there was otherwise only the indirect counter-test (2) via `euler`,
    # and two nodes compensating each other keep `euler` at 2 and stay invisible.
    #
    # It is directly decidable as soon as two nodes are joined by AT LEAST THREE
    # edges: walking the slots of `u` clockwise, the opposite slots at `v` must run
    # counter-clockwise — the two cyclic sequences are reverses of each other. With
    # exactly TWO connections (a bigon) this is combinatorially empty: two cyclic
    # sequences of length 2 are always reverses of each other. There `euler` remains
    # the only instrument — so this check supplements (2), it does not replace it.
    #
    # The case ≥ 3 is exactly the "there and back" rules (two m=3 braids over three
    # strands, R2-3/R5) and the braid–trivalent relation (R2-4/R9).
    conn = Dict{Tuple{Int,Int}, Vector{Tuple{Int,Int}}}()   # (u,v) → [(slot_u, slot_v)]
    for e in edges
        (e.a isa NodePort && e.b isa NodePort) || continue
        u, v = e.a.node, e.b.node
        u == v && continue                                   # self-loop: no orientation
        key, pair = u < v ? ((u, v), (e.a.slot, e.b.slot)) :
                            ((v, u), (e.b.slot, e.a.slot))
        push!(get!(conn, key, Tuple{Int,Int}[]), pair)
    end
    for ((u, v), ps) in conn
        length(ps) >= 3 || continue
        sort!(ps; by = first)                       # slots of u clockwise
        vs_slots = [q for (_, q) in ps]
        dv = deg(nodes[v])
        # Same direction measurement as on the leaf ring, but on the slot ring of
        # `v`: +1 = co-oriented (wrong), −1 = counter-oriented (right).
        dir = _ring_direction(vs_slots, dv)
        dir == 1 && push!(out, WiringViolation(:node_ring, v,
            "nodes n$u and n$v are joined by $(length(ps)) edges; the slots " *
            "$(first.(ps)) at n$u correspond to the slots $(vs_slots) at n$v, and " *
            "those run CO-ORIENTED. The convention requires counter-oriented " *
            "(i.e. $(reverse(vs_slots)))"))
    end

    # (2) node↔node: measured counter-test by reversing the arm sequence.
    χ = _euler_generic(word, nodes, edges, deg)
    for v in 1:length(nodes)
        deg(nodes[v]) >= 2 || continue
        e2 = _reverse_node_slots(nodes, edges, v, deg)
        χ2 = _euler_generic(word, nodes, e2, deg)
        χ2 > χ && push!(out, WiringViolation(:reversed_node, v,
            "reversing the arm sequence lifts V−E+F from $χ to $χ2 — this node " *
            "attaches to its neighbours co-oriented instead of counter-oriented"))
    end

    χ != 2 && push!(out, WiringViolation(:euler, 0,
        "V−E+F = $χ instead of 2: the embedding is not planar"))
    return out
end

"""
    check_wiring(g) -> Vector{WiringViolation}

Checks the house convention (boundary ccw, arms cw) on `g` — for `WordGraph`,
`MorphismGraph`, `CircularGraph` and `CircularMorphismGraph`. An empty list means clean.

Four kinds are reported (derivation in the comment block above this function):

* `:leaf_degree` — not exactly one edge ends at some boundary position. Not an
  orientation finding but a completeness one: `is_wired` does not check the
  boundary, and a loose leaf otherwise only shows up as "colour 0".
* `:leaf_ring` — a node hangs on ≥ 3 leaves and its slots do not run backwards
  along the leaf ring.
* `:reversed_node` — reversing this node's arm sequence RAISES `euler`; measured,
  not guessed. The node hangs mirrored.
* `:euler` — global: `V−E+F ≠ 2`.

The table form is `check_wiring_table(g)` (HTML, in the style of `wiring_table`)
and `show_wiring_check(g)` (console).
"""
check_wiring(g::WordGraph) = _check_wiring_generic(g.word, g.nodes, g.edges, _node_degree)
check_wiring(m::MorphismGraph) = check_wiring(m.graph)

const _WIRING_KIND_TXT = Dict(
    :leaf_degree   => "loose leaf",
    :leaf_ring     => "leaf↔node",
    :reversed_node => "node↔node",
    :euler         => "global")

"""
    check_wiring_table(g) -> String

The violations from `check_wiring(g)` as an HTML table, in the style of
`wiring_table` (render/tutte/Wiring.jl): one row per violation with kind, node and
the measured reason. If there is nothing to report, one green row says so.
"""
function check_wiring_table(g)
    vs = check_wiring(g)
    io = IOBuffer()
    print(io, """<table style="border-collapse:collapse;font-family:monospace;font-size:13px">""")
    print(io, """<tr><th style="text-align:left;padding:3px 10px;border-bottom:1px solid #ccc">kind</th>""" *
              """<th style="text-align:left;padding:3px 10px;border-bottom:1px solid #ccc">node</th>""" *
              """<th style="text-align:left;padding:3px 10px;border-bottom:1px solid #ccc">finding</th></tr>""")
    if isempty(vs)
        print(io, """<tr style=" background:#e8f7e8;"><td style="padding:3px 10px" colspan="3">""" *
                  """no violations — V−E+F = $(euler(g))</td></tr>""")
    else
        for v in vs
            print(io, """<tr style=" background:#ffe6e6;">""")
            print(io, """<td style="padding:3px 10px">$(get(_WIRING_KIND_TXT, v.kind, v.kind)) ⚠</td>""")
            print(io, """<td style="padding:3px 10px">$(v.node == 0 ? "—" : string(v.node))</td>""")
            print(io, """<td style="padding:3px 10px">$(v.detail)</td></tr>""")
        end
    end
    print(io, "</table>")
    return String(take!(io))
end

"""
    show_wiring_check(g; io = stdout)

Console form of `check_wiring_table` — one line per violation.
"""
function show_wiring_check(g; io::IO = stdout)
    vs = check_wiring(g)
    if isempty(vs)
        println(io, "check_wiring: no violations  (V−E+F = ", euler(g), ")")
        return nothing
    end
    println(io, "check_wiring: ", length(vs), " violation(s)")
    for v in vs
        println(io, "  ", rpad(get(_WIRING_KIND_TXT, v.kind, string(v.kind)), 14),
                v.node == 0 ? "     " : rpad(string("n", v.node), 5), v.detail)
    end
    return nothing
end
