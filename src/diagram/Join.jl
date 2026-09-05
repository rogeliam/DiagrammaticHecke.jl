# Join.jl  —  gluing two diagrams along a shared reversed subword
#
# `join_graphs(g1, A, g2, B)` glues arc A of g1's boundary to arc B of g2's boundary
# (matched reversed, colour by colour): leaf A[t] fuses with leaf B[end-t+1]. The
# glued leaves vanish and their two edges become one. This one operation gives
# morphism composition, tensor (empty arc) and partial gluing.
#
# Named `join_graphs` (NOT `join` — clashes with Base.join). GOTCHA-worthy.

# internal marker (NOT exported, never appears in a finished WordGraph's edges): an
# end of an edge that sits ON a glued leaf. `k` is the index of the g1 leaf of the
# glued PAIR (the g2 partner is translated to the same `k`), so the two edge-ends
# meeting at one glue point carry the SAME marker. See join_graphs.
struct _JoinGlue <: Port
    k::Int
end

# ---- the core gluing -------------------------------------------------------

"""
    join_graphs(g1, A, g2, B) -> WordGraph

Glue arc `A` (leaf indices of `g1`) to arc `B` (leaf indices of `g2`), matched in
REVERSE: `A[t]` ↔ `B[end-t+1]`, requiring equal colours. Glued leaves disappear and
the edge chain running through them fuses into ONE edge. The result boundary is g1's
surviving leaves (from just after A, cyclically) followed by g2's surviving leaves
(from just after B). A chain that closes up without meeting a genuine port becomes a
free monochrome circle (`Edge(colour, Circle, Circle)`), on which `GRAPH_RULES`'
`:free_circle` fires — of any length, not just the two-edge loop.
"""
function join_graphs(g1::WordGraph, A::Vector{Int}, g2::WordGraph, B::Vector{Int})
    length(A) == length(B) || error("join_graphs: arcs differ in length")
    colsA = [leaf_colour(g1, a) for a in A]
    colsB = reverse([leaf_colour(g2, b) for b in B])
    colsA == colsB || error("join_graphs: arc colours mismatch $colsA vs $colsB")

    n1 = length(g1.word); n2 = length(g2.word)
    gluedA = Set(A); gluedB = Set(B)
    # partner map: A[t] (g1 leaf) ↔ B[end-t+1] (g2 leaf)
    partner = Dict{Int,Int}()                      # g1 leaf → its g2 partner leaf
    for (t, a) in enumerate(A); partner[a] = B[length(B) - t + 1]; end

    # surviving-leaf renumbering. Surviving g1 leaves in cyclic order starting just
    # after the LAST arc leaf, then surviving g2 leaves likewise. (Any deterministic
    # order works; check_valid is the invariant.)
    surv1 = [k for k in _cyclic_from(n1, isempty(A) ? 1 : last(A) + 1) if !(k in gluedA)]
    surv2 = [k for k in _cyclic_from(n2, isempty(B) ? 1 : last(B) + 1) if !(k in gluedB)]
    new1 = Dict(k => i for (i, k) in enumerate(surv1))
    new2 = Dict(k => (length(surv1) + i) for (i, k) in enumerate(surv2))

    off = length(g1.nodes)                          # g2 node index shift
    nodes = vcat(g1.nodes, g2.nodes)

    # port translation. A surviving leaf → renumbered Leaf; a glued leaf → a marker so
    # we can fuse the two edges. Node ports just shift (g2 by off).
    GLUE = -1                                       # tag: this port is a glued leaf
    # translate a port; return (translated_port_or_original_leaf, tag, glued_leaf_index)
    tr1(p::NodePort) = (NodePort(p.node, p.slot), 0, 0)
    tr1(p::Leaf)     = haskey(new1, p.k) ? (Leaf(new1[p.k]), 0, 0) : (p, GLUE, p.k)
    tr2(p::NodePort) = (NodePort(p.node + off, p.slot), 0, 0)
    tr2(p::Leaf)     = haskey(new2, p.k) ? (Leaf(new2[p.k]), 0, 0) : (p, GLUE, p.k)

    # ---- fusing edges: CHAINS, not a single splice -------------------------
    #
    # Every glued leaf pair is a JOINT where exactly two edge ends meet (the g1 edge
    # at A[t] and the g2 edge at B[end-t+1]); each leaf carries exactly one edge in
    # its own graph. The fusion is therefore an edge chain: walk from a genuine port
    # across the joints until another genuine port comes up.
    #
    # ⚠ Handling only the SINGLE splice (one g1 edge between two joints whose other
    # neighbours are already genuine ports) plus the closed two-edge circle is not
    # enough. A LONGER chain — a g1 edge between two joints, next to it again a g2
    # edge between two joints — would write the internal marker as a port into an
    # edge and break later in the face tracer. That is exactly what happens when
    # composing two double leaves with z = ε: whole strand chains run through the
    # glue seam. The chain walk below covers both simpler cases as special cases.
    g1partner_of_g2 = Dict(v => k for (k, v) in partner)

    # segments: (colour, endA, endB); an end is a genuine port or a `_JoinGlue(k)`
    # with `k` = index of the joint's g1 leaf.
    segs = Tuple{Int,Port,Port}[]
    for e in g1.edges
        (pa, ta, ka) = tr1(e.a); (pb, tb, kb) = tr1(e.b)
        push!(segs, (e.colour,
                     ta == GLUE ? _JoinGlue(ka) : pa,
                     tb == GLUE ? _JoinGlue(kb) : pb))
    end
    for e in g2.edges
        (pa, ta, ka) = tr2(e.a); (pb, tb, kb) = tr2(e.b)
        push!(segs, (e.colour,
                     ta == GLUE ? _JoinGlue(g1partner_of_g2[ka]) : pa,
                     tb == GLUE ? _JoinGlue(g1partner_of_g2[kb]) : pb))
    end

    # joint → the (at most two) segment ends meeting there.
    at_joint = Dict{Int,Vector{Tuple{Int,Int}}}()   # k → [(segment, end 1|2), …]
    for (i, (_, ea, eb)) in enumerate(segs)
        ea isa _JoinGlue && push!(get!(at_joint, ea.k, Tuple{Int,Int}[]), (i, 1))
        eb isa _JoinGlue && push!(get!(at_joint, eb.k, Tuple{Int,Int}[]), (i, 2))
    end
    for (k, ends) in at_joint
        length(ends) == 2 || error("join_graphs: joint at leaf $k has " *
                                   "$(length(ends)) edge ends instead of 2")
    end

    # The ORDER of the emitted edges is part of the interface: edges are addressed
    # elsewhere by their INDEX (rules, tests, stored examples). So we collect
    # (position, edge) here and sort at the end by the LARGEST segment index of the
    # chain — for the simple cases (unglued g1 edge; g1 edge + g2 edge; pure g2 edge)
    # that is exactly the natural position. Without it the indices shift and e.g.
    # `circular_2parallel_apply(gD, 5, 15)` misses its target.
    out = Tuple{Int,Edge}[]
    used = falses(length(segs))
    # the other end of a segment
    other_end(i, side) = side == 1 ? segs[i][3] : segs[i][2]
    # The neighbour at joint `k`, seen from segment `cur`. The two ends of a joint
    # ALWAYS belong to different segments (one from g1, one from g2), so the
    # comparison is unambiguous.
    neighbour(k, cur) = (e = at_joint[k]; e[1][1] == cur ? e[2] : e[1])

    # OPEN CHAINS: starting from every segment with at least one genuine end.
    for i in 1:length(segs)
        used[i] && continue
        col, ea, eb = segs[i]
        (ea isa _JoinGlue && eb isa _JoinGlue) && continue     # later (circles)
        start_port = ea isa _JoinGlue ? eb : ea
        far        = ea isa _JoinGlue ? ea : eb
        used[i] = true
        cur_seg = i
        pos = i
        while far isa _JoinGlue
            j, side = neighbour(far.k, cur_seg)
            used[j] && error("join_graphs: chain runs into an already used segment")
            used[j] = true
            far = other_end(j, side)
            cur_seg = j
            pos = max(pos, j)
        end
        push!(out, (pos, Edge(col, start_port, far)))
    end

    # CLOSED CHAINS: whatever is left now hangs on joints only and closes up into a
    # free monochrome circle. It is represented like the K1 marker
    # (`Edge(colour, Circle, Circle)`, diagram/Graph.jl) so that `GRAPH_RULES`'
    # `:free_circle` can fire on it.
    for i in 1:length(segs)
        used[i] && continue
        col, ea, _ = segs[i]
        used[i] = true
        far = ea
        cur_seg = i
        pos = i
        while true
            j, side = neighbour(far.k, cur_seg)
            used[j] && break                       # circle closed
            used[j] = true
            far = other_end(j, side)
            far isa _JoinGlue || error("join_graphs: closed chain with a genuine port")
            cur_seg = j
            pos = max(pos, j)
        end
        push!(out, (pos, Edge(col, Circle(col), Circle(col))))
    end

    sort!(out, by = first)
    edges = Edge[e for (_, e) in out]

    # new boundary word from the surviving leaves, in the new numbering order
    total = length(surv1) + length(surv2)
    wordcols = zeros(Int, total)
    for k in surv1; wordcols[new1[k]] = leaf_colour(g1, k); end
    for k in surv2; wordcols[new2[k]] = leaf_colour(g2, k); end
    return WordGraph(CircularWord(wordcols), nodes, edges)
end

# cyclic index list 1..n starting at `start` (1-based, wraps).
_cyclic_from(n::Int, start::Int) = [mod1(start + j, n) for j in 0:(n - 1)]

# ---- tensor: side-by-side (empty glue) -------------------------------------

"""
    tensor(g1, g2) -> WordGraph

Horizontal composition: disjoint union with boundaries concatenated (g1's leaves 1..n
then g2's). No leaves are glued. (The general insertion-gap variant can
be added later; this is the default gap-after-n.)
"""
tensor(g1::WordGraph, g2::WordGraph) = join_graphs(g1, Int[], g2, Int[])

# ---- MorphismGraph composition ---------------------------------------------

"""
    compose(f::MorphismGraph, g::MorphismGraph) -> MorphismGraph

Compose morphisms `f : a → b` and `g : b → c` into `a → c` by gluing f's TOP arc to
g's BOTTOM arc (top(f) == bottom(g) required). The result's bottom is f's bottom, its
top is g's top.
"""
function compose(f::MorphismGraph, g::MorphismGraph)
    top(f) == bottom(g) || error("compose: codomain ≠ domain: $(top(f)) vs $(bottom(g))")
    Atop = _top_leaves(f)          # f's top arc (forward); glued reversed to g's bottom
    Bbot = _bottom_leaves(g)       # g's bottom arc (forward)
    # join_graphs matches A[t] ↔ B[end-t+1]. top(f) = reverse(colours of Atop); we need
    # top(f)==bottom(g), and bottom(g)=colours of Bbot. So glue Atop reversed to Bbot:
    joined = join_graphs(f.graph, reverse(Atop), g.graph, reverse(Bbot))
    # new cuts: the result boundary = f's surviving (bottom) leaves then g's surviving
    # (top) leaves. Bottom of result = f.bottom (all of f's survivors); top = g.top.
    nb = length(bottom(f))
    ntot = length(joined.word)
    # reduce mod ntot: if g's top is EMPTY, nb == ntot, and cut2 must equal cut1 (=0)
    # RAW for the "whole boundary, top = ε" special case in `_bottom_leaves`/
    # `_top_leaves` (MorphismGraph.jl) to trigger — those check `cut1 == cut2` without
    # reducing mod n, so a literal `cut2 = ntot` (rather than the equivalent `0`) is
    # silently WRONG. The same degeneracy arises in `tensor` above and in
    # `dot_morphism`/`cap_morphism` (LightLeaves.jl), and is handled there too.
    cut2 = ntot == 0 ? 0 : mod(nb, ntot)
    return MorphismGraph(joined, 0, cut2)     # cut1 after last leaf (=0), cut2 after bottom
end

"""
    _mirror_graph(m::MorphismGraph) -> WordGraph

Shared geometric core of `flip` and `hflip`: reflect the whole boundary across the
axis running THROUGH `m`'s two cut points (`axis = cut1+cut2+1`), i.e.
`Leaf(k) ↦ Leaf(mod1(axis - k, n))`, together with the matching braid-slot remap
(`opposite_slot`, `mod1(slot + nd.m, 2*nd.m)`) and, for `:braid` nodes with ODD `m`,
the colour swap `[s,t] → [t,s]` that keeps the parity-flipped slot's colour correct
(see `flip`'s docstring for the full derivation). `flip` and `hflip` differ ONLY in
which cut ends up `cut1` vs `cut2` afterwards, and in how they handle the
`cut1==cut2` degenerate case — both handled by the respective callers, not here.

Returns `g` UNCHANGED (n == 0, nothing to reflect) when the boundary is empty.

The `_EMPTY_BOTTOM_CUT` sentinel (`cut1 == cut2 == -1`, "bottom=ε/top=whole") is
normalised to `0` for the axis computation: both `cut1==cut2==0` ("bottom=whole/
top=ε") and `cut1==cut2==_EMPTY_BOTTOM_CUT` describe the SAME physical boundary
arc (all `n` leaves, read as one whole arc) with only the bottom/top ROLE
differing, so they must reflect through the SAME axis to be inverse to each other
(needed for `flip`/`hflip`'s degenerate branches).
"""
function _mirror_graph(m::MorphismGraph)
    g = m.graph
    n = length(g.word)
    n == 0 && return g
    c1 = m.cut1 == _EMPTY_BOTTOM_CUT ? 0 : m.cut1
    c2 = m.cut2 == _EMPTY_BOTTOM_CUT ? 0 : m.cut2
    axis = c1 + c2 + 1
    # MIRRORING, not rotation. The boundary is mirrored (`Leaf(k) ↦ Leaf(axis − k)`,
    # below), so a genuine mirror must also REVERSE the arm order at every braid
    # node, not merely rotate it: a rotation (`opposite_slot`, `slot + m`) preserves
    # the orientation instead of flipping it, and `flip` then produced a node of the
    # opposite handedness to the rest of the diagram.
    #
    # In analogy with the boundary rule `i ↦ n+1−i`, the slot rule is
    # `slot ↦ 2m − slot` (cyclically): slot `2m` stays fixed, the rest reverse.
    # THE CONSTANT: a braid slot colour depends only on the PARITY of the
    # slot (`_slot_colour`, diagram/Graph.jl), and `slot ↦ c − slot` preserves parity
    # exactly for EVEN `c`. With `c = 2m` all slot colours therefore stay put, so the
    # colour swap `[s,t] → [t,s]` for odd `m`, which would compensate the parity
    # shift of `+m`, is unnecessary.
    newnodes = g.nodes
    function remap(p::Port)
        if p isa Leaf
            return Leaf(mod1(axis - p.k, n))
        elseif p isa NodePort
            nd = g.nodes[p.node]
            # TRIVALENTS MUST BE MIRRORED TOO. A trivalent does have an
            # orientation: the CYCLIC ORDER of its three arms. Leaving it alone made
            # exactly this node run against the rest of the diagram after `flip`
            # (arms 1→2→3 counter-clockwise while the unmirrored half has them
            # clockwise). Irrelevant for `display_tutte` (which recomputes the
            # embedding), but every renderer reading the stored slot order as a
            # rotation system then turns the node.
            #
            # Rule `slot ↦ mod1(3 − slot, 3)`: 1↔2, 3 fixed. That reverses the
            # orientation while keeping the ROLES intact — in `merge_morphism` slots
            # 1,2 are the two lower legs and slot 3 the upper output, and the
            # merge/cap semantics hang on that role. (`slot ↦ 4 − slot` would swap a
            # leg with the output and is therefore wrong.) A `:dot` has one arm only
            # and stays unchanged.
            nd.kind === :trivalent && return NodePort(p.node, mod1(3 - p.slot, 3))
            nd.kind === :braid || return p          # dot: nothing to mirror
            return NodePort(p.node, mod1(2 * nd.m - p.slot, 2 * nd.m))
        else
            return p                                # Circle
        end
    end
    newedges = [Edge(e.colour, remap(e.a), remap(e.b)) for e in g.edges]
    oldcols  = [leaf_colour(g, k) for k in 1:n]
    newword  = CircularWord([oldcols[mod1(axis - i, n)] for i in 1:n])
    return WordGraph(newword, newnodes, newedges)
end

"""
    flip(m::MorphismGraph) -> MorphismGraph

The VERTICAL mirror of `m`: reflect the diagram top↔bottom. If `m : a → b`
(drawn bottom = `a`, top = `b`), then `flip(m) : b → a` with the LETTERS of each
boundary word kept in order — `bottom(flip(m)) == top(m)` and `top(flip(m)) ==
bottom(m)`. `flip` is an involution (on endpoints).

Why it is NOT just a cut swap: `bottom` reads
its arc forward but `top` reads its arc reversed (morphism/PathToGraph.jl), so merely swapping
`cut1`/`cut2` reads the same physical arc the other way and REVERSES the word.
A true vertical mirror needs THREE things:
  1. reflect the cyclic BOUNDARY across the axis through the two cut points —
     leaf `k ↦ (cut1+cut2+1) − k (mod n)` — this fixes the boundary WORDS;
  2. reflect each BRAID node about its own axis — slot `j ↦ j+m (mod 2m)` =
     `opposite_slot`. A braid's slot colour follows the slot PARITY
     (odd → colours[1], even → colours[2]). When `m` is ODD, `+m` flips the
     parity, so slot `j+m` carries the WRONG colour unless the node's colours
     are also SWAPPED, `[s,t] → [t,s]`. When `m` is even the parity is preserved and
     no swap is needed. This colour swap is what keeps every glued edge
     colour-consistent;
  3. reflect each TRIVALENT node about its own axis — slot `j ↦ 3 − j (mod 3)`,
     i.e. swap the two bottom legs (slots 1,2) and fix the top (slot 3). A
     trivalent's orientation IS the cyclic order of its three arms, so leaving it
     unmirrored would run it the wrong way round relative to the mirrored rest
     of the diagram. Only `:dot` nodes, with their single arm, have no
     orientation to mirror.
Then swap the cuts. `flip` keeps `is_planar_braid` holding, EVERY edge colour
matching its endpoint slot/leaf colours, each slot/leaf used once, and satisfies
`flip∘flip = id`.

**The degenerate branch** (`cut1==cut2`, bottom=ε/top=whole or
bottom=whole/top=ε, `n>1`): `bottom` reads its arc forward but `top` reads its
arc REVERSED (`MorphismGraph.jl`), so even though the physical leaf set is the
same, returning the SAME graph unmirrored would make `top(flip(m))` read out
`reverse(bottom(m))` instead of `bottom(m)` itself — the required contract
`top(flip(m)) == bottom(m)` (letters in the SAME order) only holds that way
when the word happens to be a palindrome. So the axis
reflection (`_mirror_graph`) DOES need to run here too — it performs exactly the
letter-order reversal `top`'s reversed read requires to cancel out — using the
axis `cut1+cut2+1` normalised so the `_EMPTY_BOTTOM_CUT` sentinel (`-1`) and the
plain `0` degenerate cut agree on the SAME axis (see `_mirror_graph`'s own
docstring). Only the CUT handling still differs from the generic branch: swapping
`cut1`/`cut2` raw is a no-op when they're equal, so the `_EMPTY_BOTTOM_CUT`
sentinel is toggled instead, on top of the (always-run) reflection.
`bottom(flip(m)) == top(m)` and `top(flip(m)) == bottom(m)` hold letter-for-letter
on the degenerate cases too (including round-tripping both sentinel directions),
and `flip(hflip(m)) == hflip(flip(m))` holds there as well.
"""
function flip(m::MorphismGraph)
    n = length(m.graph.word)
    n == 0 && return MorphismGraph(m.graph, m.cut2, m.cut1)
    newgraph = _mirror_graph(m)
    if m.cut1 == m.cut2
        # degenerate input (bottom=whole/top=ε, or bottom=ε/top=whole): swapping
        # cut1/cut2 raw is a NO-OP here (they're equal), so it would silently return
        # the SAME degenerate reading instead of the other one. Toggle the
        # `_EMPTY_BOTTOM_CUT` sentinel instead. The reflection itself (`newgraph`,
        # above) DOES need to have run — see the note above; it is what makes
        # `top`'s reversed read recover the original letter order.
        newcut = m.cut1 == _EMPTY_BOTTOM_CUT ? 0 : _EMPTY_BOTTOM_CUT
        return MorphismGraph(newgraph, newcut, newcut)
    end
    return MorphismGraph(newgraph, m.cut2, m.cut1)
end

"""
    hflip(m::MorphismGraph) -> MorphismGraph

The HORIZONTAL mirror of `m`: reflect the diagram left↔right (as opposed to `flip`'s
top↔bottom). If `m : x → y` (bottom = `x`, top = `y`), then `hflip(m) : reverse(x) →
reverse(y)` — bottom stays bottom and top stays top, only each word's LETTER ORDER
reverses:

    bottom(hflip(m)) == reverse(bottom(m))
    top(hflip(m))    == reverse(top(m))

A genuine left↔right mirror
must fix bottom-is-bottom and top-is-top and only reverse each word's letters —
merely swapping `cut1`/`cut2` does not achieve this: it swaps which word is
bottom and which is top, i.e. it is `flip` composed with a true horizontal
mirror, not a horizontal mirror on its own.

Geometrically: `flip` reflects across the axis BETWEEN the two cuts (through the
gaps, `axis = cut1+cut2+1`), swapping which cut becomes `cut1` vs `cut2`; `hflip`
reflects across the axis THROUGH the two cuts, keeping each cut's role. Both use the
SAME boundary reflection (`_mirror_graph`, defined above): `Leaf(k) ↦
Leaf(mod1(axis - k, n))`, the braid slot remap `opposite_slot` (`mod1(slot + nd.m,
2*nd.m)`), and the colour swap `[s,t] → [t,s]` for `:braid` nodes with ODD `m`
(needed because `+m` flips the slot's parity, and a braid's slot colour follows slot
parity — see `flip`'s own derivation above, unchanged). `hflip` needs this FULL
remap: because it actually MOVES the boundary (reflecting through the cuts moves
every leaf that isn't one of the two cut points themselves), omitting the remap
would leave edge colours inconsistent with their new endpoints. The two mirrors
differ ONLY in the final cut assignment:
`flip` returns `MorphismGraph(newgraph, cut2, cut1)` (cuts swapped), `hflip` returns
`MorphismGraph(newgraph, cut1, cut2)` (cuts kept).

Degenerate case (`cut1 == cut2`, "bottom=whole/top=ε" or "bottom=ε/top=whole",
possibly via the `_EMPTY_BOTTOM_CUT` sentinel): the remap DOES need to run here
too — bottom stays bottom and top stays top, but `reverse(bottom(m))` for a
nonempty, non-palindromic `bottom(m)` is a genuinely different letter sequence,
and only the axis reflection produces it (see `_mirror_graph`'s
sentinel-normalisation note). The cuts themselves are NOT toggled (unlike
`flip`) — bottom stays bottom, top stays top, and reversing `ε` is still `ε`,
so the sentinel (if present) passes through as-is; only the underlying graph
changes.

Verified on all A₃ `double_leaves` of word pairs of length ≤ 3, plus
the degenerate `flip(dot_morphism(1))`: `bottom(hflip(m)) == reverse(bottom(m))`,
`top(hflip(m)) == reverse(top(m))`, `hflip` is an involution,
`is_planar_braid(hflip(m).graph)` holds, every edge colour still matches its
endpoint slot/leaf colour, and `degree(hflip(m)) == degree(m)` (dots/trivalents
untouched).

`flip(hflip(m)) == hflip(flip(m))` holds in ALL cases checked, including the
`cut1==cut2` degenerate boundary.
"""
function hflip(m::MorphismGraph)
    n = length(m.graph.word)
    n == 0 && return MorphismGraph(m.graph, m.cut1, m.cut2)
    newgraph = _mirror_graph(m)
    return MorphismGraph(newgraph, m.cut1, m.cut2)
end

