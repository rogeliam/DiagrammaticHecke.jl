# test/circularbraidedges.jl — the edge patterns of the braid rules.
#
# This file builds NO rule. It checks two facts about the braid rules:
#
#   §3  C6 and C14 are ONE relation. A dot at the 8-armed `:gbraid` can be computed
#       via a PREIMAGE (`expand_gbraid`, C13⁻¹, then the existing rules WITHOUT C14),
#       and that gives the same result as C14.
#
#   §1  The unfolding route is a valid CROSS-CHECK ORACLE for C15: the same pattern
#       computed once via `_fr_braid_on_gen12` and once via the preimages gives the
#       same result.
#
# One site per fact, not the full series.

using Test
using DiagrammaticHecke
using DiagrammaticHecke: Edge, NodePort, Leaf, CIRCULAR_RULES, expand_gbraid,
                     _fr_dot_on_gen12, _fr_gen12_merge, _fr_braid_on_gen12,
                     _fr_dot_into_braid, CircularComboR

# ---- shared helpers -----------------------------------------------------------

_fbe_graph(t) = t isa CircularGraph ? t : t.graph
_fbe_key(c) = sort([(string(co), circular_canonical_key(_fbe_graph(t))) for (t, co) in pairs_of(c)],
                   by = string)
_fbe_without(names...) = [r for r in CIRCULAR_RULES if !(r.name in names)]

"Full reduction of every term of a combination (the basis for comparison)."
function _fbe_voll(c)
    out = CircularComboR()
    for (t, co) in pairs_of(c)
        out = out + co * reduce_circular_full(t)
    end
    return out
end

"The planar unfoldings of the 8-armed `:gbraid` `v` (4 preimages × 8 rotations)."
function _fbe_preimages(g::CircularGraph, v::Int)
    out = CircularGraph[]
    for variant in 1:4, rot in 0:7
        e = expand_gbraid(g, v, variant; rot = rot)
        e === nothing && continue
        (is_wired(e) && euler(e) == 2 && isempty(check_wiring(e))) || continue
        push!(out, e)
    end
    return out
end

"""
The preimage route at site `v`: reduce every planar unfolding with the existing rules
WITHOUT `raus`, keep only the weight-lowering ones, and fully re-reduce each
result. Returns `(count planar, count lowering, set of keys)` — the set must be
a SINGLETON, otherwise the right-hand side is not independent
of the preimage.
"""
function _fbe_preimage_route(g::CircularGraph, v::Int, raus::Symbol...)
    w0 = circular_weight(g)
    rules = _fbe_without(raus...)
    exps = _fbe_preimages(g, v)
    keys = Set(); lowering = 0
    for e in exps
        c, _ = reduce_circular(circular_decorated(e); rules = rules)
        isempty(pairs_of(c)) && continue
        all(circular_weight(t.graph) < w0 for (t, _) in pairs_of(c)) || continue
        lowering += 1
        push!(keys, _fbe_key(_fbe_voll(c)))
    end
    return (length(exps), lowering, keys)
end

# ---- the fixtures -------------------------------------------------------------

"The preimage of the 8-armed `:gbraid`: two `:braid` + two trivalents (as in test/circulardotongen12.jl)."
_fbe_preimage8() = CircularGraph(CircularWord([2, 3, 2, 3, 2, 3, 2, 3]),
    CircularNode[circular_node([2, 3, 2, 3, 2, 3]), circular_node([2, 3, 2, 3, 2, 3]),
            circular_node([3, 3, 3]), circular_node([3, 3, 3])],
    Edge[Edge(3, NodePort(1, 4), NodePort(3, 1)), Edge(3, NodePort(3, 3), NodePort(2, 2)),
         Edge(2, NodePort(1, 5), NodePort(2, 1)), Edge(3, NodePort(1, 6), NodePort(4, 1)),
         Edge(3, NodePort(4, 2), NodePort(2, 6)), Edge(2, Leaf(1), NodePort(1, 1)),
         Edge(3, Leaf(2), NodePort(4, 3)), Edge(2, Leaf(3), NodePort(2, 5)),
         Edge(3, Leaf(4), NodePort(2, 4)), Edge(2, Leaf(5), NodePort(2, 3)),
         Edge(3, Leaf(6), NodePort(3, 2)), Edge(2, Leaf(7), NodePort(1, 3)),
         Edge(3, Leaf(8), NodePort(1, 2))])

"Dot on boundary leaf `k`."
function _fbe_dot(g::CircularGraph, k::Int)
    lets = letters(g.word); c = lets[k]
    nodes = vcat(copy(g.nodes), CircularNode[circular_node([c])]); di = length(nodes)
    shift(p) = p isa Leaf ? (p.k == k ? NodePort(di, 1) : Leaf(p.k > k ? p.k - 1 : p.k)) : p
    CircularGraph(CircularWord(vcat(lets[1:k-1], lets[k+1:end])), nodes,
             Edge[Edge(e.colour, shift(e.a), shift(e.b)) for e in g.edges])
end

"The bare 8-armed `:gbraid` (as in test/circularbraidongen12.jl)."
function _fbe_bare(s::Int, t::Int)
    a = [isodd(i) ? s : t for i in 1:8]; lf(i) = i == 1 ? 1 : 10 - i
    word = [a[i] for i in 1:8][invperm([lf(i) for i in 1:8])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:8])
end

"Attaches a 6-armed `:braid` to the leaves `k, k+1, k+2` (as in test/circularbraidongen12.jl)."
function _fbe_with_braid(g::CircularGraph, k::Int)
    lets = letters(g.word); n = length(lets)
    ks = [mod1(k + i, n) for i in 0:2]
    c1, c2 = lets[ks[1]], lets[ks[2]]
    for start in 1:6, dir in (1, -1), newdir in (1, -1)
        barms = [isodd(i) ? c1 : c2 for i in 1:6]
        sl = [mod1(start + dir * (i - 1), 6) for i in 1:3]
        all(barms[sl[i]] == lets[ks[i]] for i in 1:3) || continue
        rest = [s for s in 1:6 if !(s in sl)]
        ns = newdir == 1 ? rest : reverse(rest)
        rw = [lets[mod1(k + 2 + i, n)] for i in 1:(n - 3)]
        neu = vcat(rw, [barms[s] for s in ns])
        leafmap = Dict(mod1(k + 2 + i, n) => i for i in 1:(n - 3))
        nodes = vcat(copy(g.nodes), CircularNode[circular_node(barms)]); bi = length(nodes)
        slot_of_old = Dict(ks[i] => sl[i] for i in 1:3)
        rp(p) = p isa Leaf ? (haskey(slot_of_old, p.k) ? NodePort(bi, slot_of_old[p.k]) :
                                                         Leaf(leafmap[p.k])) : p
        edges = Edge[Edge(e.colour, rp(e.a), rp(e.b)) for e in g.edges]
        for (j, s) in enumerate(ns)
            push!(edges, Edge(barms[s], Leaf(length(neu) - 3 + j), NodePort(bi, s)))
        end
        h = CircularGraph(CircularWord(neu), nodes, edges)
        (is_wired(h) && euler(h) == 2 && isempty(check_wiring(h))) && return h
    end
    return nothing
end

# ⚠ the DOT POLICY blocks C14 and the preimage stage on an
# 8-armed node with ONE dot. This block tests the edge patterns of the
# surgery, not the gate — hence with the switch turned off.
# The §3 block measures that all 32 unfoldings lower the weight and equal C14.
# That is the world `CIRCULAR_TRIVALENT_MERGE = false`: C14 live, C19 dead, arm
# weight. Under `true` none of those unfoldings lowers the region measure any more.
_fbe_old_policy = CIRCULAR_DOT_POLICY[]
_fbe_old_direction = CIRCULAR_TRIVALENT_MERGE[]
CIRCULAR_DOT_POLICY[] = false
CIRCULAR_TRIVALENT_MERGE[] = false
try

@testset "edge-pattern measurements on the braid rules" begin

    G8 = first(pairs_of(_fr_gen12_merge(_fbe_preimage8())))[1]

    @testset "§3 — C6 and C14 are ONE relation (preimage route == C14)" begin
        # Over all eight leaves there are 32 planar unfoldings each, all 32
        # weight-lowering, with ONE `circular_canonical_key` each, and that key is
        # C14's. Leaf 1 stands here as the regression anchor.
        gd = _fbe_dot(G8, 1)
        @test gd.nodes[1].kind === :braid && arm_count(gd.nodes[1]) == 8

        c14 = _fr_dot_on_gen12(gd)
        @test c14 !== nothing
        @test length(collect(pairs_of(c14))) == 1          # C14: ONE term

        (planar, lowering, keys) = _fbe_preimage_route(gd, 1, :dot_on_gen12)
        @test planar == 32
        @test lowering == 32
        @test length(keys) == 1                            # PREIMAGE-INDEPENDENT
        @test only(keys) == _fbe_key(_fbe_voll(c14))       # and == C14

        # The apparent difference "C6 two terms, C14 one term" is an artefact of
        # where the count is taken: C6 acts BEFORE, C14 AFTER the further reduction.
        # Fully reduced, the two routes agree.
        u = _fbe_dot(_fbe_preimage8(), 1)                    # dot on the PREIMAGE ⇒ C6
        c6 = _fr_dot_into_braid(u)
        @test c6 !== nothing
        @test length(collect(pairs_of(c6))) == 2           # C6: TWO terms
        @test _fbe_key(reduce_circular_full(u)) == _fbe_key(_fbe_voll(c14))
    end

    @testset "§1 — the unfolding route as a cross-check oracle for C15" begin
        # Over both colour pairs × all eight sites there are 32 planar unfoldings
        # each, of which 24 lower the weight, with ONE key each, equal to C15's
        # right-hand side. One site stands here.
        g = _fbe_with_braid(_fbe_bare(1, 2), 1)
        @test g !== nothing
        c15 = _fr_braid_on_gen12(g)
        @test c15 !== nothing

        (planar, lowering, keys) = _fbe_preimage_route(g, 1, :braid_on_gen12, :gen12_on_gen12)
        @test planar == 32
        @test lowering == 24
        @test length(keys) == 1
        @test only(keys) == _fbe_key(_fbe_voll(c15))
        # C15 is the workhorse; the unfolding route stands ALONGSIDE it as a
        # cross-check, which is the role of this test.
    end
end

finally
    CIRCULAR_DOT_POLICY[] = _fbe_old_policy
    CIRCULAR_TRIVALENT_MERGE[] = _fbe_old_direction
end

