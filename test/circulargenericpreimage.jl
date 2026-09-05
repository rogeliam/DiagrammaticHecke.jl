# test/circulargenericpreimage.jl — the GENERIC PREIMAGE STAGE in the circular driver
# (src/circular/rules/CircularDriver.jl, `_fr_generic_preimage`/`circular_generic_preimage`).
#
# The stage: if NO rule from `CIRCULAR_RULES` matches, an 8-armed `:gbraid` is unfolded via
# `expand_gbraid` (C13⁻¹) and the full registry is computed on it; only what is strictly
# lighter and label-free is taken (the three acceptance conditions of
# `circular_via_preimage`).
#
# It does NOT stand in `CIRCULAR_RULES` but is the fallback in
# `_circular_first_rule_match` — which is why the tests below check the two separately: the
# stage directly (it should reproduce the right-hand sides of the special rules) and the
# driver (where NOTHING may change, every existing rule keeping precedence).

using Test
using DiagrammaticHecke
using DiagrammaticHecke: _fr_generic_preimage, _fr_dot_on_gen12, _circular_gen12_expansions,
                     _CIRCULAR_PREIMAGE_DEPTH, Edge, NodePort, Leaf, CIRCULAR_RULES

"The bare 8-armed `:gbraid` (as in test/circulargen12twoedges.jl)."
function _gp_bare(s::Int, t::Int)
    a = [isodd(i) ? s : t for i in 1:8]; lf(i) = i == 1 ? 1 : 10 - i
    word = [a[i] for i in 1:8][invperm([lf(i) for i in 1:8])]
    CircularGraph(CircularWord(word), CircularNode[circular_node(a)],
             Edge[Edge(a[i], Leaf(lf(i)), NodePort(1, i)) for i in 1:8])
end

"""
Attach a node with arms `barms` to the `j` consecutive leaves `k, …, k+j-1` (as
`_c18_with_partner` in test/circulargen12twoedges.jl; with `j = 1` and `barms = [c]` this is
exactly a DOT at leaf `k`). The first planar
wiring found, else `nothing`.
"""
function _gp_with_partner(g::CircularGraph, k::Int, j::Int, barms::Vector{Int})
    lets = letters(g.word); n = length(lets); m = length(barms)
    ks = [mod1(k + i - 1, n) for i in 1:j]
    for start in 1:m, dir in (1, -1), newdir in (1, -1)
        sl = [mod1(start + dir * (i - 1), m) for i in 1:j]
        all(barms[sl[i]] == lets[ks[i]] for i in 1:j) || continue
        rest = [s for s in 1:m if !(s in sl)]
        ns = newdir == 1 ? rest : reverse(rest)
        rw = [lets[mod1(ks[end] + i, n)] for i in 1:(n - j)]
        neu = vcat(rw, [barms[s] for s in ns])
        leafmap = Dict(mod1(ks[end] + i, n) => i for i in 1:(n - j))
        nodes = vcat(copy(g.nodes), CircularNode[circular_node(barms)]); bi = length(nodes)
        slot_of_old = Dict(ks[i] => sl[i] for i in 1:j)
        rp(p) = p isa Leaf ? (haskey(slot_of_old, p.k) ? NodePort(bi, slot_of_old[p.k]) :
                                                         Leaf(leafmap[p.k])) : p
        edges = Edge[Edge(e.colour, rp(e.a), rp(e.b)) for e in g.edges]
        for (jj, s) in enumerate(ns)
            push!(edges, Edge(barms[s], Leaf(n - j + jj), NodePort(bi, s)))
        end
        h = CircularGraph(CircularWord(neu), nodes, edges)
        (is_wired(h) && euler(h) == 2 && isempty(check_wiring(h))) && return h
    end
    return nothing
end

"Dot at leaf `k` of the bare `:gbraid` — the C6/C14 site."
_gp_dot(k::Int) = _gp_with_partner(_gp_bare(1, 2), k, 1, [letters(_gp_bare(1, 2).word)[k]])

# ⚠ the DOT POLICY (`circular_dots_reducible`) blocks C14 and the
# generic preimage stage on a `2k >= 8`-node with ONE dot — there the dot
# stays put. This block tests the SURGERY, not the gate, and hence
# runs with the switch turned off. The gate itself has its own
# tests in `test/circulardotmerge.jl`.
# The block checks the preimage stage against C14 as ground truth: arm weights
# (8,2,9) and (6,2,9), and the expectation that the driver takes `:dot_on_gen12`.
# That is the world `CIRCULAR_TRIVALENT_MERGE = false`: C14 live, C19 dead, arm
# weight. Under `true` the C14 step raises the region measure and the driver
# leaves the single dot where it is.
_old_policy = CIRCULAR_DOT_POLICY[]
_old_direction = CIRCULAR_TRIVALENT_MERGE[]
CIRCULAR_DOT_POLICY[] = false
CIRCULAR_TRIVALENT_MERGE[] = false
try

@testset "generic preimage stage (CircularDriver.jl)" begin

    @testset "no site ⇒ nothing" begin
        # no `:gbraid` in the diagram: the stage touches nothing at all.
        @test _fr_generic_preimage(circular(braid(1, 2; m = 3))) === nothing
        @test circular_generic_preimage(circular(braid(1, 2; m = 3))) === nothing
    end

    @testset "a site exists but nothing lighter ⇒ a fixed point stays a fixed point" begin
        # The bare `:gbraid` is a normal form; all 32 planar unfoldings lead back to it
        # (C13 rebuilds it) and none lowers the weight. The guard refuses — the stage is
        # a no-op here.
        g = _gp_bare(1, 2)
        @test length(_circular_gen12_expansions(g, 1)) == 32
        @test _fr_generic_preimage(g) === nothing
        c, hist = reduce_circular(circular_decorated(g))
        @test isempty(hist)
        @test length(collect(pairs_of(c))) == 1
    end

    @testset "C6/C14 is CONTAINED in the preimage route" begin
        # The claim "one relation" stands or falls with the preimage route giving
        # the same right-hand side as the special rule C14. One leaf is written out,
        # the other seven only counted.
        h = _gp_dot(1)
        @test h !== nothing
        @test circular_weight(h) == (8, 2, 9)

        res = _fr_generic_preimage(h)
        @test res !== nothing
        @test length(collect(pairs_of(res))) == 1
        @test circular_weight(first(pairs_of(res))[1]) == (6, 2, 9)     # = the C14 computation
        c14 = _fr_dot_on_gen12(h)
        @test c14 !== nothing
        @test circular_canonical_key(first(pairs_of(res))[1]) ==
              circular_canonical_key(first(pairs_of(c14))[1])

        # all eight leaves: ONE term each, and each equal to C14
        matches = 0
        for k in 1:8
            hk = _gp_dot(k); hk === nothing && continue
            r = _fr_generic_preimage(hk); a = _fr_dot_on_gen12(hk)
            (r === nothing || a === nothing) && continue
            sort([circular_canonical_key(t) for (t, _) in pairs_of(r)], by = string) ==
                sort([circular_canonical_key(t) for (t, _) in pairs_of(a)], by = string) &&
                (matches += 1)
        end
        @test matches == 8
    end

    @testset "C14 fires directly through the driver" begin
        # With the dot policy off (top of the file), C14 is reachable directly through
        # `CIRCULAR_RULES`; the driver fires it without falling back to the generic
        # preimage stage. The stage remains available as a fallback — that it delivers
        # the same key is shown in the test block above.
        _, hist = reduce_circular(circular_decorated(_gp_dot(1)))
        @test :dot_on_gen12 in hist
    end

    @testset "a restricted registry does not see the stage" begin
        # The fallback hangs on `rules === CIRCULAR_RULES`: a restricted registry gets
        # exactly the restricted computation, with no silent change of meaning.
        h = _gp_dot(1)
        without14 = [r for r in CIRCULAR_RULES if r.name !== :dot_on_gen12]
        _, hist = reduce_circular(circular_decorated(h); rules = without14)
        @test !(:generic_preimage in hist)
    end

    @testset "recursion guard: the depth stays 1" begin
        # The counter is 0 before and after each run, and the stage's intermediate
        # reduction cannot call it again (a fresh rule array).
        @test _CIRCULAR_PREIMAGE_DEPTH[] == 0
        reduce_circular(circular_decorated(_gp_bare(1, 2)))
        @test _CIRCULAR_PREIMAGE_DEPTH[] == 0
        _fr_generic_preimage(_gp_dot(1))
        @test _CIRCULAR_PREIMAGE_DEPTH[] == 0
    end

    @testset "maxsites caps it (and `0` is the off switch)" begin
        h = _gp_dot(1)
        @test _fr_generic_preimage(h; maxsites = 0) === nothing
        @test _fr_generic_preimage(h; maxsites = 1) !== nothing
        alt = CIRCULAR_GENERIC_PREIMAGE_MAXSITES[]
        try
            CIRCULAR_GENERIC_PREIMAGE_MAXSITES[] = 0
            @test circular_generic_preimage(h) === nothing
        finally
            CIRCULAR_GENERIC_PREIMAGE_MAXSITES[] = alt
        end
        @test CIRCULAR_GENERIC_PREIMAGE_MAXSITES[] == 4     # the default value
    end

    @testset "circular_via_preimage accepts `except_rule = nothing`" begin
        # Without `except_rule` the combinator computes with the FULL registry,
        # which the generic stage needs, since a rule name would also switch off the
        # non-circular path of the same rule (see the file header of
        # CircularViaPreimage.jl).
        h = _gp_dot(1)
        r1 = circular_via_preimage(h, [1], v -> _circular_gen12_expansions(h, v))
        r2 = _fr_generic_preimage(h)
        @test r1 !== nothing && r2 !== nothing
        @test sort([circular_canonical_key(t) for (t, _) in pairs_of(r1)], by = string) ==
              sort([circular_canonical_key(t) for (t, _) in pairs_of(r2)], by = string)
    end
end

finally
    CIRCULAR_DOT_POLICY[] = _old_policy
    CIRCULAR_TRIVALENT_MERGE[] = _old_direction
end

