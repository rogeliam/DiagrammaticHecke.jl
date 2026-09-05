# circular/rules/CircularDecoratedRules.jl — D4 ("dot + connection") on CircularDecorated/CircularGraph.
#
# The D4 mechanism on the circular world: the relation
# `dot_A = dot_B + dot_C - BC_decorated_trivalent`, with the matching criterion
# "neighbour cell closer to the marking, edge of the same colour", on
# `CircularGraph`/`CircularDecorated`/`CircularMorphismGraph`.
#
# "DOT" IN THE CIRCULAR PICTURE. In the plain picture a dot is its own node kind
# (`nd.kind === :dot`, degree 1). The circular picture has no fixed node kind for
# it — "1 arm = hollow dot": a `CircularNode` with
# `arm_count == 1` plays exactly the same role, regardless of its `kind`
# (`:mono` for colour 2, `:mixed` for colour 1/3). `find_circular_d4_match`/
# `apply_circular_d4` therefore match on `arm_count(nd) == 1` instead of
# `nd.kind === :dot`.
#
# "TRIVALENT" IN THE CIRCULAR PICTURE. The new node that `apply_d4` builds for the
# BC term is not its own degree-3 node kind here, but the special case
# `circular_node([i,i,i])` of the general circular constructor (3 arms, same colour `i`)
# — `circular_node` recognizes this automatically as `:mono` (colour 2) or
# `:mixed` (colour 1/3).
#
# ⚠ The DRIVER — `reduce_to_circular_leave`
# with its growth guard, its finishing step and `circular_fusion_step` — lives in
# circular/rules/CircularLeaveDriver.jl, included directly after this file. What stays here
# is the D4 relation itself (`apply_circular_d4`, `find_circular_d4_match`, the cyclic end
# order) plus the scalar extraction.
#

"""
    CIRCULAR_ZAMO_ENABLED

Measurement switch (Zamo from the pairing): when `false`,
`reduce_to_circular_leave` skips BOTH Zamo variants (whole-term step
`circular_zamo_step` and region/direction step `circular_zamo_direction_step`) —
reduction then runs entirely over the remaining rules ("Zamo-free"). Default
is `true`; normal operation is unaffected. Toggle with
`CIRCULAR_ZAMO_ENABLED[] = false` and reset afterwards (always bracket
measurements in `try … finally`).
"""
const CIRCULAR_ZAMO_ENABLED = Ref(true)

"""
    CircularDecoratedMorphism

A `CircularMorphismGraph` (with cuts/marking) plus one polynomial label per inner
cell — the guided circular variant of the decorated diagram, analogous to
`DecoratedMorphism`. A thin pairing, no new degree of freedom.

The third, optional argument is the label of the OUTER region (`outer_label`,
default `1`) — the same field as in
`CircularDecorated`. **It MUST be passed along at every re-pairing**: the driver
rebuilds the morphism per term and step from `(m, dd.region_labels)`, and
without this third argument the outer label would silently drop out at each
of the ~15 call sites.
"""
struct CircularDecoratedMorphism
    m::CircularMorphismGraph
    region_labels::Vector{SoergelPoly}
    outer_label::SoergelPoly

    function CircularDecoratedMorphism(m::CircularMorphismGraph, labels::Vector{SoergelPoly},
                                  outer::SoergelPoly = one(SoergelPoly))
        length(labels) == region_count(m.graph) || throw(ArgumentError(
            "expected $(region_count(m.graph)) labels (one per region), got $(length(labels))"))
        return new(m, labels, outer)
    end
end

"The underlying decorated CircularGraph (cuts forgotten)."
_circular_decorated(fdm::CircularDecoratedMorphism) =
    CircularDecorated(fdm.m.graph, fdm.region_labels, fdm.outer_label)

"The label of the OUTER region (see `CircularDecorated`)."
outer_label(fdm::CircularDecoratedMorphism) = fdm.outer_label

"""
    with_marks(dd::CircularDecorated, m::CircularMorphismGraph) -> CircularDecoratedMorphism
    with_marks(g::CircularGraph, m::CircularMorphismGraph)      -> CircularDecoratedMorphism

Re-attach the BOUNDARY MARKING to a reduction RESULT — for notebooks and
measurements that draw intermediate or final states.

The marking is NOT a property of the diagram but of the cut position: it lives
as `cut1`/`cut2` on the [`CircularMorphismGraph`](@ref). A `CircularDecorated` (let alone
a bare `CircularGraph`) doesn't know it — `display`ing one directly gives a picture
WITHOUT the marking, which is missing information rather than a renderer bug.
`reduce_to_circular_leave`/`reduce_circular` return exactly such mark-less objects.

The boundary word and cut positions are unchanged by reduction (every rule is a
relation between diagrams with the SAME boundary), so it is safe to reuse the
cuts of the starting position `m` — that's exactly what this function does.
"""
with_marks(dd::CircularDecorated, m::CircularMorphismGraph) =
    CircularDecoratedMorphism(CircularMorphismGraph(dd.graph, m.cut1, m.cut2),
                         dd.region_labels, dd.outer_label)
with_marks(g::CircularGraph, m::CircularMorphismGraph) = with_marks(circular_decorated(g), m)

"""
    circular_decorated_svg(fdm::CircularDecoratedMorphism; size = 600, r = 0.34·size) -> String

Like `circular_decorated_svg(::CircularDecorated)` (render/CircularCells.jl), but with
everything that follows from `fdm.m` — analogous to
`decorated_svg(::DecoratedMorphism)`:

- BOTH boundary markers, `_circular_left_mark(fdm.m)` (black) and
  `_circular_right_mark(fdm.m)` (grey) — the grey one matters beyond decoration:
  `circular_tutte_positions` calls `_add_nudge!(mark_between2)` on it, which affects
  leaf placement;
- `bottom_leaves`/`top_leaves`, i.e. the ORIENTATION of the leaf ring (they
  determine leaf-angle assignment, `_halfcircle_angles`) — without them the
  decorated circular path would draw nodes mirrored relative to
  `circular_tutte_svg(fdm.m)`;
- the region distances `circular_region_distances(fdm.m)`, always (not only in
  debug mode): they pick the one region in which a face's polynomial is
  displayed (`_label_region_of_face`, render/Cells.jl). Debug mode
  (`debug_labels()`) only gates whether each region additionally shows the
  tuple `(number, polynomial, distance)`.
"""
circular_decorated_svg(fdm::CircularDecoratedMorphism; size::Int = 600, r::Real = 0.34 * size) =
    circular_decorated_svg(_circular_decorated(fdm); size = size, r = r,
                      mark_between = _circular_left_mark(fdm.m),
                      mark_between2 = _circular_right_mark(fdm.m),
                      bottom_leaves = _bottom_leaves(fdm.m),
                      top_leaves = _top_leaves(fdm.m),
                      distances = circular_region_distances(fdm.m),
                      region_based = true)

"""
    display_circular_decorated(fdm::CircularDecoratedMorphism)

Draws the decorated circular morphism (Tutte) with cell labels AND the left
boundary marker. Analogous to `display_decorated(::DecoratedMorphism)`.
"""
display_circular_decorated(fdm::CircularDecoratedMorphism) = (display(fdm); nothing)

# MIME methods for BOTH types. Without them `display(fdm)` falls back to text
# form, and VS Code shows no image at all — it skips `image/svg+xml` and needs
# PNG (see the comment in render/Decorated.jl:25-30, which documents this
# exact issue). `show_circular_combo` draws its terms through this
# path so that markers, orientation, and distances all come from ONE source.
Base.show(io::IO, ::MIME"image/svg+xml", fdm::CircularDecoratedMorphism) =
    print(io, circular_decorated_svg(fdm))
Base.show(io::IO, ::MIME"image/png", fdm::CircularDecoratedMorphism) =
    write(io, _svg_to_png(circular_decorated_svg(fdm)))
Base.showable(::MIME"image/svg+xml", ::CircularDecoratedMorphism) = true
Base.showable(::MIME"image/png", ::CircularDecoratedMorphism) = true

"""
    circular_outer_region_distance(m::CircularMorphismGraph) -> Int

Distance of the OUTER region from the black marking — the circular analogue of
[`outer_cell_distance`](@ref). `0` whenever there is a
marking, `-1` when there is none (no marking, no distance).

Careful with "no marking": a `CircularMorphismGraph` always HAS one as long as its
boundary word is non-empty — `cut1 = cut2 = 0` is the canonical gap between leaf
`n` and leaf 1, not the absence of a mark (`_circular_left_mark`). Only an empty
boundary word gives `-1` here, exactly as in [`circular_region_distances`](@ref). The
genuinely unmarked object is a bare `CircularDecorated`, which has no cuts and hence
never reaches this function.

**Why a constant ("strict wall").** The
boundary circle is a wall one cannot walk through, so the outer region is not
reachable from the inside at all and a BFS would be the wrong picture. But the
marking itself sits in a boundary GAP, and the outer region lies directly in
front of that gap: it touches the marking without crossing anything. Hence its
distance is `0` by definition, not by search.

**Consequence.** [`circular_extract_scalars`](@ref) always pulls a non-constant
`outer_label` into the coefficient once the diagram carries a marking: a label
on the outside IS a global scalar. Nothing gets stuck out there any more, and
no rule ever has to move a label across the boundary — the fusion stays strictly
inside ([`circular_fusion_step`](@ref)).

**The inner distances keep their own rule.** The inner distances ([`circular_region_distances`](@ref))
keep the wall: boundary regions of a component that the marking cannot reach
through edges get `-1`, and their labels travel through the component
chaining ([`circular_extract_floating_components`](@ref)).
"""
function circular_outer_region_distance(m::CircularMorphismGraph)
    g = m.graph
    n = length(g.word)
    (n == 0 || _circular_left_mark(m) === nothing) && return -1
    return 0
end

"""
    circular_extract_scalars(fdm::CircularDecoratedMorphism) -> Tuple{SoergelPoly, CircularDecoratedMorphism}

Pulls scalar
factors out of the cell labels — (a) any constant ≠ 1 in an arbitrary region,
and (b) any label (even non-constant) of a region at marking-distance 0.
NECESSARY here (unlike the nine base rules in CIRCULAR_RULES, which know no cell
labels at distance 0): `reduce_circular`/`reduce_circular_full` are NOT marking-aware (no
`extract_scalars` analogue in CircularDriver.jl — a bare `CircularDecorated` has no cuts),
so they would e.g. never pull an `alpha_i` label on the D4 trivalent term's own
(distance-0) cell into the coefficient. `reduce_to_circular_leave` therefore calls
this function explicitly before every `reduce_circular` round, exactly where
`apply_d4`'s plain analogue (`reduce_decorated`) handles it implicitly via
`extract_scalars`.

**Distances over REGIONS**, matching `find_circular_d4_match`'s convention: a face
made of several regions would otherwise get one distance per the face notion
but several per the region notion, letting a label get pulled in wrongly (or
not). Rule, storage, and extraction all use the same region-based notion.

**The OUTER region counts too.** For `outer_label` the same
rule applies, with [`circular_outer_region_distance`](@ref) as the distance: a
constant always moves into the factor, a genuine polynomial only at distance 0
— a scalar may only be pulled out of the diagram if it sits
in the region at distance 0.

**This distance is ALWAYS 0 once there is a marking** — the outer region
touches the marking's boundary gap. An outer label in the marked driver
therefore always becomes part of the coefficient; only an UNMARKED
`CircularDecorated` (distance `-1`) leaves it in place.
"""
function circular_extract_scalars(fdm::CircularDecoratedMorphism)
    dist = circular_region_distances(fdm.m)
    factor = one(SoergelPoly)
    labels = copy(fdm.region_labels)
    for c in 1:length(labels)
        f = labels[c]
        isone(f) && continue
        if degree(f) == 0 || dist[c] == 0
            factor = factor * f
            labels[c] = one(SoergelPoly)
        end
    end
    outer = fdm.outer_label
    if !isone(outer) && (degree(outer) == 0 || circular_outer_region_distance(fdm.m) == 0)
        factor = factor * outer
        outer = one(SoergelPoly)
    end
    return factor, CircularDecoratedMorphism(fdm.m, labels, outer)
end

# ---- D4 — the "dot + connection" relation, on CircularGraph ------------------------
#
# CYCLIC ORDER OF THE THREE STRAND ENDS.
#
# The restriction "all three strand ends must be boundary leaves" is NOT a
# real condition of the rule, only an ARTIFACT of using `sort(by = p -> p.k)`
# for the cyclic order: for boundary leaves that shortcut works directly. A
# `NodePort` has no `p.k`, but its position in the cyclic order is still
# well-defined — it just has to be read off GEOMETRICALLY instead.
#
# WHERE THE ORDER COMES FROM. The three strand ends (`dot_strand_end`, `endA`,
# `endB`) all lie on the boundary of EXACTLY TWO cells: the dot's own cell
# (`dot_cell` — the capping edge has the same cell on both sides) and the "BC"
# cell (the other side of `ei`). Removing the capping edge AND `ei` (exactly
# what term 1 does) merges these two cells into ONE — the interior of the
# new 3-armed node plus the area around it. The cyclic order sought is
# exactly the order in which
# `dot_strand_end`/`endA`/`endB` follow each other along the boundary of THIS
# merged region — for boundary leaves that is just the leaf order (they sit on
# the diagram's own boundary); for a `NodePort` the position follows from the
# neighbouring-node geometry, with no extra case needed.
#
# HOW THIS IS COMPUTED: `_trace_generic` (diagram/Faces.jl) keeps, for every
# cell, its FACE CYCLE as an ordered dart list (already the basis of
# `_circular_face_walks`, circular/CircularFaces.jl). Taking the face cycles of `dot_cell` and
# `bc_cell` (or just one, if both coincide — happens for small diagrams) and
# removing the four darts of the two vanishing edges (capping edge there/back,
# `ei` there/back), the cycles break into at most four open chains whose ends
# are exclusively `dot_strand_end`/`endA`/`endB` (special case: since the dot
# node always has exactly 1 arm, both darts of the capping edge sit right next
# to each other in the same cycle — the "chain" between them is then empty,
# see `_circular_d4_cyclic_ends`). Concatenating these chains in the order where one
# chain's end matches the next one's start gives the combined boundary — and
# in it, `dot_strand_end`/`endA`/`endB` already appear in the cyclic order
# sought (up to a global rotation, irrelevant for `circular_node([i,i,i])` as a node
# type, but significant for the CONCRETE wiring).
#
# WHY THIS THROWS RATHER THAN GUESSES: a wrong order produces a non-planar
# node, which does NOT show up as an exception but as a silently wrong cell
# structure. `_circular_d4_cyclic_ends` therefore
# throws (e.g. on unexpected chain topology) instead of guessing an order.
#
# CHORD-CASE EXCEPTION: if `ei` is a CHORD of the dot cell — both sides of `ei` are
# the same cell (`dot_cell == bc_cell_old`, e.g. when the diagram has only one
# inner cell) — removing the four darts does not merge two cells; instead the
# ONE face cycle splits into three chains that pair up as a 2-cycle
# {dse <-> X} plus a fixed point {Y -> Y} rather than a 3-cycle. This case is
# still planar and handled: the new trivalent node goes into the pocket
# containing the dot (bounded by the dse->X chain, the capping edge, and `ei`)
# and splits that pocket into 2 regions — valid, because the cyclic order of
# the three ends around the pocket boundary is still well-defined: (dse, X,
# Y). `_circular_d4_cyclic_ends` recognizes this signature (only possible when
# fi_dot == fi_bc) and returns exactly this order; any OTHER unexpected
# topology still throws.

"""
    _circular_d4_cyclic_ends(g, dot_node, capport, dse, endA, endB) -> Vector{Port}

Determines the cyclic order of the three strand ends `dse` (`dot_strand_end`),
`endA`, `endB` for the new 3-armed circular node that `apply_circular_d4`'s term 1
builds — the general replacement for `sort([dse, endA, endB]; by = p -> p.k)`,
which also works for `NodePort` ends (see the section header above for the
geometric derivation). `capport` is the dot's own `NodePort` (slot 1).

Returns the three ports in ONE cyclic rotation, starting at `dse` (the
concrete rotation is irrelevant for `circular_node([i,i,i])`, but the relative
order fixes the wiring). Throws if the chain topology matches neither the
expected 3-cycle form nor the chord signature (2-cycle {dse <-> X} plus fixed
point {Y -> Y} at `dot_cell == bc_cell_old`, see the section header) — no
fallback to guessing.
"""
function _circular_d4_cyclic_ends(g::CircularGraph, dot_node::Int, capport::NodePort,
                             dse::Port, endA::Port, endB::Port,
                             dot_cell::Int, bc_cell_old::Int)
    t = _trace_generic(g.word, g.nodes, g.edges, _circular_degree)
    fi_dot = findfirst(==(dot_cell), t.cell_of_face)
    fi_bc  = findfirst(==(bc_cell_old), t.cell_of_face)
    (fi_dot === nothing || fi_bc === nothing) && error(
        "_circular_d4_cyclic_ends: dot_cell/bc_cell not found in the face cycle " *
        "— unexpected")

    # the four darts of the two vanishing edges (capping edge, `ei` —
    # referenced via their end ports, not via an edge index this function
    # doesn't know)
    cap_d_out = t.port_dart[capport]        # (node) -> dse
    cap_d_in  = t.darts[cap_d_out].rev      # dse -> (node)
    ei_d_a    = t.port_dart[endA]           # endA -> endB
    ei_d_b    = t.darts[ei_d_a].rev         # endB -> endA
    vanishing = Set([cap_d_out, cap_d_in, ei_d_a, ei_d_b])

    # each vanishing dart fixes which strand end it OPENS a new chain segment
    # towards (`after_port`, dart points IN the direction of that end) and
    # which end it CLOSES a chain segment coming FROM (`before_port`, dart
    # arrives FROM that end).
    before_port = Dict(cap_d_in => dse, ei_d_a => endA, ei_d_b => endB)
    after_port  = Dict(cap_d_out => dse, ei_d_a => endB, ei_d_b => endA)

    # split a face cycle at the vanishing darts into open chains (start-end
    # pairs of the three strand ends); if two vanishing darts sit right next
    # to each other (dot with only 1 arm: capping edge there and back with no
    # real segment in between), an EMPTY chain results (start==end).
    function chains_of(cyc)
        vpos = [k for k in 1:length(cyc) if cyc[k] in vanishing]
        out = Tuple{Port,Port}[]
        for idx in 1:length(vpos)
            d_here = cyc[vpos[idx]]
            haskey(after_port, d_here) || continue   # this dart opens no segment
            start_port = after_port[d_here]
            d_next = cyc[vpos[mod1(idx + 1, length(vpos))]]
            haskey(before_port, d_next) || error(
                "_circular_d4_cyclic_ends: unexpected chain topology on the face cycle")
            push!(out, (start_port, before_port[d_next]))
        end
        return out
    end

    # FACES FROM THE DARTS, not from the cell numbers: if a cell consists of
    # SEVERAL faces, `findfirst(==(cell), t.cell_of_face)` would pick an
    # arbitrary one, whereas the four vanishing darts determine their face
    # uniquely. This covers 1 face (chord case), 2 (normal case), and 3
    # uniformly.
    face_of = Dict{Int,Int}()
    for (fi, cyc) in enumerate(t.cycles), d in cyc
        face_of[d] = fi
    end
    faces = unique(face_of[d] for d in (cap_d_out, cap_d_in, ei_d_a, ei_d_b))
    chains = vcat((chains_of(t.cycles[fi]) for fi in faces)...)
    by_start = Dict{Port,Port}()
    for (s, e) in chains
        haskey(by_start, s) && error(
            "_circular_d4_cyclic_ends: two chains start at the same strand end — " *
            "chain topology not unambiguous")
        by_start[s] = e
    end
    length(by_start) == 3 || error(
        "_circular_d4_cyclic_ends: expected exactly 3 chain segments (one per strand end), " *
        "got $(length(by_start))")

    # chain together starting at `dse`, until the cycle closes back at `dse`.
    #
    # CHORD CASE (verified): if `ei` is a CHORD of the dot cell
    # (dot_cell == bc_cell_old, both ei darts in the same face cycle), the
    # chaining does not split into a 3-cycle but into a 2-cycle {dse <-> X}
    # plus a fixed point {Y -> Y} ({X, Y} = {endA, endB}). This is NOT a
    # reason to throw here: the new trivalent node goes into the POCKET
    # containing the dot (bounded by the dse->X chain + capping edge + ei);
    # this pocket is split into 2 regions by the trivalent connection, which
    # is perfectly fine planarly. The cyclic order of the ends around the
    # pocket boundary is thus well-defined: (dse, X, Y) — the dse->X chain is
    # the pocket wall, Y sits on the pocket boundary too as the second `ei`
    # end. Proceed in this order like in the 3-cycle case (same direction
    # convention below); ALL other unexpected topologies still throw.
    order = Port[]
    X = get(by_start, dse, nothing)
    chord = fi_dot == fi_bc && X !== nothing && X != dse &&
            get(by_start, X, nothing) == dse
    if chord
        ys = [p for p in keys(by_start) if p != dse && p != X]
        chord = length(ys) == 1 && by_start[only(ys)] == only(ys)
        chord && (order = Port[dse, X, only(ys)])
    end
    if !chord
        cur = dse
        for _ in 1:3
            push!(order, cur)
            haskey(by_start, cur) || error(
                "_circular_d4_cyclic_ends: chain concatenation breaks off — unexpected topology")
            cur = by_start[cur]
        end
        cur == dse || error(
            "_circular_d4_cyclic_ends: chain concatenation does not close back at " *
            "dot_strand_end after 3 steps — unexpected topology")
    end
    # `order` runs in the direction "chain follows chain" — the OPPOSITE
    # direction of the clockwise slot order sought (measured against
    # `sort(by = p -> p.k)` on all three boundary-leaf cases). Only reverse
    # the direction, do NOT shift the starting point: `order[1]` is already
    # `dse` and must stay that way — a bare `reverse(order)` would push `dse`
    # to the end, which is WRONG: the concrete rotation matters, since it
    # fixes which cell ends up carrying `alpha_i`, not just the cyclic
    # class. So: keep the first element (`dse`), reverse the REST.
    result = vcat(order[1], reverse(order[2:end]))

    # Fix the CONCRETE rotation (not just the cyclic CLASS): `result` is
    # already one of the three possible rotations of the cyclic order sought
    # — all three are equally VALID/planar for the resulting
    # `circular_node([i,i,i])` (rotation-symmetric as a node type), but give
    # DIFFERENT concrete cell ids for the new "BC" cell (cell ids are pure
    # numbering artefacts, see diagram/Faces.jl, but `apply_circular_d4` still has
    # to make ONE concrete choice). So that the pure boundary-leaf case
    # reproduces `sort(by = p -> p.k)` LITERALLY, not just up to rotation
    # (regression anchor, see test/circulardecoratedrules.jl), we rotate here so
    # that — if at least one end is a boundary leaf — the leaf with the
    # SMALLEST leaf index comes first (exactly what `sort` does anyway). If
    # all three ends are `NodePort`s there is no numeric comparison basis;
    # the geometrically determined rotation is kept as-is (each of the three
    # is equally valid, see above).
    leaf_idx = findall(p -> p isa Leaf, result)
    isempty(leaf_idx) && return result
    start = leaf_idx[argmin([result[j].k for j in leaf_idx])]
    return vcat(result[start:end], result[1:start-1])
end

"""
    apply_circular_d4(fd::CircularDecorated, dot_node::Int, ei::Int) -> Union{Nothing, CircularComboR}

The D4 rule ("dot + connection", see `apply_d4`) on
`CircularDecorated`. `dot_node` must be a circular node with `arm_count == 1` (the circular
dot); `ei` an edge of the SAME colour as the dot that is NEITHER the dot's own
capping edge NOR touches the dot node itself on either end (the "connection of
the other two"). Replaces the dot + this edge by three terms (analogous to
`dot_A = dot_B + dot_C - BC_decorated_trivalent`):
1. a new 3-armed circular node (`circular_node([i,i,i])`) at the same three strand ends,
   with `alpha_i` on the cell that `ei` bounds in `g`
   (coefficient `-1`);
2. a circular dot (`circular_node([i])`) on one end of `ei`, the other end connected to
   the existing dot strand by a fresh edge (coefficient `+1`);
3. symmetrically with the other end of `ei` (coefficient `+1`).
`nothing` if `dot_node` is not a degree-1 circular dot, `ei` is its own capping
edge, or `ei` touches the dot node itself.

GENERAL CASE: the three strand ends may be arbitrary `Port`s (boundary leaves
AND/OR `NodePort`s at other nodes) — restricting to boundary leaves is not a
real condition of the rule, only an artifact of sorting by `p.k` for the
cyclic slot order, see [`_circular_d4_cyclic_ends`](@ref) and the section header of
this file. Throws (rather than guessing) when the chain topology on the face
cycle does not match the expected form.
"""
function apply_circular_d4(fd::CircularDecorated, dot_node::Int, ei::Int)
    g = fd.graph
    1 <= dot_node <= length(g.nodes) || return nothing
    nd = g.nodes[dot_node]
    arm_count(nd) == 1 || return nothing
    at = _circular_edges_at_node(g, dot_node)
    length(at) == 1 || return nothing
    (cap_ei, _capport, dot_strand_end) = at[1]   # dot_strand_end = the other end of the capping edge
    i = arm_colour(nd, 1)
    ei == cap_ei && return nothing
    1 <= ei <= length(g.edges) || return nothing
    e = g.edges[ei]
    e.colour == i || return nothing
    # `ei` must not touch the dot node itself (otherwise it's the capping
    # edge or a degenerate case)
    _touches_node(p::NodePort, ni::Int) = p.node == ni
    _touches_node(::Port, ::Int) = false
    (_touches_node(e.a, dot_node) || _touches_node(e.b, dot_node)) && return nothing

    endA, endB = e.a, e.b   # the two ends of the connecting edge ("B" and "C")

    # the "BC" cell: the side of `ei` that is NOT the dot's own cell — needed
    # both for the cyclic order (below) and for the label transfer (term 1),
    # hence computed here already.
    l, r = g.cells.port_cell[e.a]
    dot_cell, _ = g.cells.port_cell[NodePort(dot_node, 1)]
    bc_cell_old = l == dot_cell ? r : l
    bc_cell_old >= 1 || error(
        "apply_circular_d4: the BC cell lies in the outer face — not supported")

    # CYCLIC ORDER of the three strand ends (general case, see the section
    # header above): `_circular_d4_cyclic_ends` generalizes this GEOMETRICALLY over
    # the face cycles of `dot_cell`/`bc_cell_old` and reproduces EXACTLY (not
    # just up to rotation) `sort(by = p -> p.k)` when all three ends are
    # boundary leaves (regression-guarded in test/circulardecoratedrules.jl).
    #
    # Some edge choices are geometrically inadmissible (e.g. the second bigon
    # edge, if the other edge would first simplify to a needle = 0). In these
    # cases `_circular_d4_cyclic_ends` throws — D4 is then not a valid move here,
    # and this returns `nothing` instead of propagating an error.
    ends_sorted = try
        _circular_d4_cyclic_ends(g, dot_node, NodePort(dot_node, 1), dot_strand_end,
                            endA, endB, dot_cell, bc_cell_old)
    catch e
        e isa ErrorException || rethrow()
        return nothing
    end

    # The cyclic order can collapse on inadmissible topologies (e.g.
    # [Leaf(3), Leaf(3), Leaf(3)] on the second bigon edge). The new trivalent
    # node would then not be planar; D4 is not a move here.
    length(unique(ends_sorted)) == 3 || return nothing

    # ---- Term 1: new 3-armed circular node at {dot_strand_end, endA, endB},
    #      alpha_i on the cell that `ei` bounds in `g` -----------
    gplus1 = CircularGraph(g.word, vcat(g.nodes, [circular_node([i, i, i])]), g.edges)
    tid = length(gplus1.nodes)
    keep1 = [ee for (j, ee) in enumerate(gplus1.edges) if j != cap_ei && j != ei]
    for (slot, p) in enumerate(ends_sorted)
        push!(keep1, Edge(i, p, NodePort(tid, slot)))
    end
    g1 = _circular_delete_nodes(gplus1, Set([dot_node]), keep1)

    # ---- LABELS: over REGIONS, not over faces --------------------------------
    #
    # `find_circular_d4_match` decides over regions; a face made of several regions
    # would otherwise let a single `alpha_i` sit on several image-separated
    # areas, since face and region distances can disagree exactly where a dot
    # disappears.
    #
    # The BC REGION is the Dot region's neighbour beyond `ei` — the same one
    # the rule matches on in the first place.
    R_dot = circular_region_of_dot(g, dot_node)
    adj, _ = circular_region_adjacency(g)
    bc_hit = [R2 for (R2, ej) in adj[R_dot] if ej == ei]
    isempty(bc_hit) && error(
        "apply_circular_d4: `ei` does not separate the dot region $R_dot from any " *
        "neighbour region — unexpected (find_circular_d4_match requires exactly that)")
    R_bc = first(bc_hit)

    rn1 = regions(g1)
    bc_gaps = regions(g)[R_bc].gaps
    if !isempty(bc_gaps)
        bc_new = findfirst(S -> bc_gaps[1] in S.gaps, rn1)
        bc_new === nothing && error(
            "apply_circular_d4: BC gap $(bc_gaps[1]) is in no region of the new graph")
    else
        # BC REGION WITHOUT BOUNDARY GAPS (an error class): not findable via
        # the boundary gap, but findable GEOMETRICALLY:
        # after inserting the new trivalent node, the BC region is the sector
        # BETWEEN the two `ei` spokes (the arms to `endA`/`endB`) — bounded by
        # the `endA`/`endB` arms, their neighbouring arms at the carrier nodes,
        # and the new trivalent node. That is exactly the trivalent sector
        # "clockwise after `slotA`" when the next slot clockwise is `slotB`
        # (otherwise the other way round); the third sector (with `dse`)
        # belongs to the split dot region. Agrees with the gap-based path
        # above on the boundary-leaf case (regression-checked via the
        # `_circular_d4_bc_new` comparison).
        slotA = findfirst(==(endA), ends_sorted)
        slotB = findfirst(==(endB), ends_sorted)
        (slotA === nothing || slotB === nothing) && error(
            "apply_circular_d4: endA/endB not in ends_sorted — unexpected")
        tid1 = length(g1.nodes)   # the trivalent node is the last node of g1
                                  # (only dot_node died, tid placed at the end)
        sr1 = circular_node_sector_regions(g1)[tid1]
        # CYCLIC successor of `slotA` among the THREE slots (1→2, 2→3, 3→1).
        # `mod1(slotA, 3) + 1` is the same for slots 1 and 2 but yields 4 for
        # slot 3 — a `BoundsError` on the 3-element `ends_sorted` (the minimal
        # figure is the bigon fixture in test/circulardecoratedrules.jl).
        # `mod1(slotA + 1, 3)` wraps correctly.
        bc_new = ends_sorted[mod1(slotA + 1, 3)] == endB ? sr1[slotA] : sr1[slotB]
        (bc_new < 1 || bc_new > length(rn1)) && error(
            "apply_circular_d4: BC sector region $bc_new invalid — unexpected")
    end
    # The new 3-armed node SPLITS the dot region into two. Only the BC side
    # carries `alpha_i`; the other half has no predecessor and gets `1` from
    # the transfer.
    labels1 = _circular_transfer_labels(g, g1, fd.region_labels;
                                   overrides = Dict(bc_new => alpha(i)))
    # The outer label is carried through UNCHANGED: D4
    # only rebuilds the interior, the outer region stays the outer region.
    # Without this third argument it would silently drop out here.
    term1 = -one(SoergelPoly) *
            CircularCombo{SoergelPoly}(CircularDecorated(g1, labels1, fd.outer_label))

    # ---- Term 2/3: dot on endA resp. endB, the other two strand ends
    #      (dot_strand_end and the respective other end of `ei`) connected ----
    isone(fd.region_labels[R_dot]) || error(
        "apply_circular_d4: the dot's own region carries a non-trivial label " *
        "— unexpected (the dot should be a pure bimodule capping)")
    function dotted_term(keepend::Port, otherend::Port)
        gplus = CircularGraph(g.word, vcat(g.nodes, [circular_node([i])]), g.edges)
        did = length(gplus.nodes)
        keep = [ee for (j, ee) in enumerate(gplus.edges) if j != cap_ei && j != ei]
        push!(keep, Edge(i, keepend, NodePort(did, 1)))
        push!(keep, Edge(i, dot_strand_end, otherend))
        g2 = _circular_delete_nodes(gplus, Set([dot_node]), keep)
        # Over regions no special case is needed: the dot region merges with
        # its neighbour (the rule "the non-trivial label wins", see
        # `_circular_transfer_labels`), everything else passes through 1:1.
        labels2 = _circular_transfer_labels(g, g2, fd.region_labels)
        return CircularCombo{SoergelPoly}(CircularDecorated(g2, labels2, fd.outer_label))
    end
    term2 = dotted_term(endA, endB)
    term3 = dotted_term(endB, endA)

    return term1 + term2 + term3
end



"""
    find_circular_d4_match(m::CircularMorphismGraph) -> Union{Nothing, Tuple{Int,Int}}

Searches for a D4 match by this criterion: go through all circular dots
(`arm_count == 1`). Check whether there is a neighbouring **region** with an
edge of the same colour that has a smaller distance to the marking.

Distances are computed over **regions** (not faces): an edge separates two
regions exactly when it borders no dot. This prevents two image-separated
boundary areas that coincide under the face notion from necessarily carrying
the same distance.

DEVIATION FROM `find_d4_match`: the dots are not
searched in node/edge order, but first sorted by the distance of their
**region** to the marking — DESCENDING.
Reason: the termination measure `Σ_dots dist(region)` should strictly
decrease per D4 step (see `reduce_to_circular_leave`). Within one dot, the
neighbour regions/edges are sorted by edge index.
"""
function find_circular_d4_match(m::CircularMorphismGraph)
    g = m.graph
    dist = circular_region_distances(m)
    all(==(-1), dist) && return nothing              # no marking
    adj, _ = circular_region_adjacency(g)

    # region of each dot
    dot_region = Dict{Int, Int}()
    for (ni, nd) in enumerate(g.nodes)
        arm_count(nd) == 1 || continue
        dot_region[ni] = circular_region_of_dot(g, ni)
    end

    # largest region distance first; node index as tiebreak
    candidates = sort(collect(dot_region); by = t -> (-dist[t[2]], t[1]))
    for (ni, R) in candidates
        nd = g.nodes[ni]
        col = arm_colour(nd, 1)
        for (R2, ei) in adj[R]
            dist[R2] < dist[R] || continue
            g.edges[ei].colour == col || continue
            return (ni, ei)
        end
    end
    CIRCULAR_D4_WORD_TIEBREAK[] || return nothing
    return _find_circular_d4_word_tiebreak(m, dist, adj, candidates)
end

"""
    CIRCULAR_D4_WORD_TIEBREAK

**D4 without a distance gain**, decided by the REGION WORD. Default `true`;
`false` restricts D4 to a region of strictly smaller distance.

**Why.** A dot can sit in a region whose neighbours across its own colour all
have the SAME distance — then the first pass finds no gain and the term stays a
circular leaf, even though its region word says otherwise. The case that shows
it: a 2-dot in a region of distance 1 whose region word is `22` (not reduced),
between two arms of a 5-armed `mono[2,2,2,2,2]`. Dot fusion hands exactly this
shape (`case == :direct`) back to D4, so without the tiebreak nobody acts.

**The rule.** Only when the first pass found nothing: among the same candidates
(largest distance first), take the dots whose region word
(`circular_region_words`, the level-BFS reading) is NOT reduced, and fire towards
the neighbour across an edge of the dot colour with the SAME distance whose
region word is strictly SHORTER — smallest length first, then smallest edge
index.

**Termination.** Lexicographic `(Σ_dots dist, Σ_dots |word|)`: this pass keeps
`Σ dist` and lowers the word length of the moved dot (the new dot terms sit at
the same or a closer region, the trivalent term has one dot fewer), while the
first pass lowers `Σ dist`.
"""
const CIRCULAR_D4_WORD_TIEBREAK = Ref(true)

"""
    _find_circular_d4_word_tiebreak(m, dist, adj, candidates) -> Union{Nothing, Tuple{Int,Int}}

The second pass of [`find_circular_d4_match`](@ref) — see
[`CIRCULAR_D4_WORD_TIEBREAK`](@ref) for what it decides and why.
"""
function _find_circular_d4_word_tiebreak(m::CircularMorphismGraph, dist, adj, candidates)
    g = m.graph
    words, _ = circular_region_words(m)
    for (ni, R) in candidates
        w = words[R]
        (w === nothing || is_reduced(w)) && continue
        col = arm_colour(g.nodes[ni], 1)
        best = nothing
        for (R2, ei) in adj[R]
            dist[R2] == dist[R] || continue
            g.edges[ei].colour == col || continue
            w2 = words[R2]
            (w2 === nothing || length(w2) >= length(w)) && continue
            key = (length(w2), ei)
            (best === nothing || key < best[1]) && (best = (key, ei))
        end
        best === nothing || return (ni, best[2])
    end
    return nothing
end

"""
    find_circular_d4_any(m::CircularMorphismGraph) -> Union{Nothing, NamedTuple}

The COMBINED D4 candidate search over both variants: the edge variant
[`find_circular_d4_match`](@ref) and the node variant
[`find_circular_d4_node_match`](@ref). Picks the candidate with the **largest
distance of the dot region**; ties go to the **edge variant** (the proven
one).

Returns `(kind = :edge, dot_node, ei)` resp. `(kind = :node, dot_node, v, p,
q, insert_after)`, or `nothing` if neither matches.

Why this is already the global order: each of the two searches already sorts
its own dots by DESCENDING region distance and returns the best one — so the
maximum of the two winners is the maximum over all candidates.
"""
function find_circular_d4_any(m::CircularMorphismGraph)
    g = m.graph
    dist = circular_region_distances(m)
    e_hit = find_circular_d4_match(m)
    n_hit = find_circular_d4_node_match(m)
    e_dist = e_hit === nothing ? -1 : dist[circular_region_of_dot(g, e_hit[1])]
    n_dist = n_hit === nothing ? -1 : dist[n_hit.R]
    if e_hit !== nothing && (n_hit === nothing || e_dist >= n_dist)
        return (kind = :edge, dot_node = e_hit[1], ei = e_hit[2])
    elseif n_hit !== nothing
        return (kind = :node, dot_node = n_hit.dot_node, v = n_hit.v,
                p = n_hit.p, q = n_hit.q, insert_after = n_hit.insert_after)
    end
    return nothing
end

"Applies the match returned by [`find_circular_d4_any`](@ref)."
_apply_circular_d4_hit(fd::CircularDecorated, h) =
    h.kind === :edge ? apply_circular_d4(fd, h.dot_node, h.ei) :
                       apply_circular_d4_node(fd, h.dot_node, h.v, h.p, h.q, h.insert_after)
