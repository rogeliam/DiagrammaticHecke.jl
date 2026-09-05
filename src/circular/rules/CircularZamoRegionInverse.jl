# circular/rules/CircularZamoRegionInverse.jl — Zamo triangles/regions (the
# corner-continuation count), the backward direction choice (Z1⁻¹/Z2⁻¹, its acceptance
# guard and log), and region-signature diagnostics.
#
# Included directly after CircularZamoRegionStep.jl.

# ---- Zamo triangles and Zamo regions ---------------------------------------
#
# THE DIFFERENCE, IN WORDS:
#
#   A "Zamo triangle" is a triangle in which every edge colour occurs once —
#   that definition stays. But a Zamo REGION is not just any connected
#   collection of Zamo triangles — the nodes further out must also continue
#   the SAME pattern. A Zamo triangle whose 1-edge and 3-edge meet at a
#   `[1,1,3,3]` node does NOT count towards the region. This fixes the
#   overcounting in `Zamo(1,8)`.
#
# So: `zamo_triangles` = the purely local condition (inner region, boundary
# length 3, three distinct colours). `zamo_regions` = the subset of those whose
# corners continue the pattern.
#
# THE EXCLUSION CRITERION. A `[1,1,3,3]` node is a node of the 1/3 world whose
# arm sequence splits into TWO colour blocks — unlike the alternating
# crossing shape `[1,3,1,3]` (four blocks). The difference is exactly that
# between "two strands cross" (the pattern continues outward) and "two
# strands are merged here and turn back" (P1-merged, the pattern ends). If a
# triangle's 1-edge meets its 3-edge at such a node, the triangle doesn't
# count.
#
# The count that matters is that of the END TERM, not of the fixture:
#
#   | Graph                            | inner reg. | triangles (naive) | Zamo regions (corrected) |
#   |-----------------------------------|------------|--------------------|----------------------------|
#   | fixtures (1,8)/(2,9)/(8,1)/(9,2)   |      6     |         6          |             6              |
#   | Zamo(1,8), end term               |      8     |         8          |           **6**            |
#   | Zamo(1,7), end term               |      8     |         8          |             5              |
#   | Zamo(1,10), end term (without direction choice) | 13 |  13         |            10              |
#
# (With the direction choice, `Zamo(1,10)` has no such end term: 10 ≥ 7 ⇒ Z1⁻¹
# fires, then `CIRCULAR_RULES`, giving 6 terms with ≤ 2 Zamo regions.)
#
# The fixtures stay at 6 — the criterion cuts nothing genuine — while
# `Zamo(1,8)` falls from 8 to 6.

# Colour blocks of the cyclic arm sequence: [1,1,3,3] -> 2, [1,3,1,3] -> 4,
# [1,1,1] -> 1.
function _fzr_colour_blocks(a::Vector{Int})
    n = length(a)
    n <= 1 && return 1
    return count(k -> a[k] != a[mod1(k + 1, n)], 1:n)
end

# A "`[1,1,3,3]`-like" node: two-coloured from the 1/3 world, arms in EXACTLY
# TWO blocks. `[1,3,1,3]` (a crossing) is NOT one, nor is `[1,1,1]`
# (single-coloured).
_fzr_grouped_mixed(nd::CircularNode) =
    nd.kind === :mixed && length(unique(arms(nd))) == 2 && _fzr_colour_blocks(arms(nd)) == 2

"""
    zamo_triangles(g::CircularGraph) -> Vector{Tuple{Int,Vector{CircularBoundaryLetter}}}

The **Zamo triangles**: inner regions of boundary length 3 whose three
boundary letters carry three DIFFERENT colours (each edge colour exactly
once). This is the definition kept unchanged from the original —
purely local, without looking at the corners.

For the count the direction choice depends on, see [`zamo_regions`](@ref):
it discards the triangles whose pattern ends at a corner.
"""
function zamo_triangles(g::CircularGraph)
    out = Tuple{Int,Vector{CircularBoundaryLetter}}[]
    for (R, ls) in _fzr_inner(g)
        length(ls) == 3 || continue
        length(unique(l.colour for l in ls)) == 3 || continue
        push!(out, (R, ls))
    end
    return out
end

# Does the triangle continue its pattern outward at ALL three corners?
function _fzr_continues_outward(g::CircularGraph, ls::Vector{CircularBoundaryLetter})
    n = length(ls)
    for j in 1:n
        c = ls[j].corner
        c.vertex[1] === :node || return false      # corner at a leaf: no pattern outward
        _fzr_grouped_mixed(g.nodes[c.vertex[2]]) || continue
        # The grouped node is only a problem if the two DIFFERENTLY coloured
        # triangle edges meet here — exactly the "1-edge and 3-edge meet at a
        # [1,1,3,3] node" case.
        ls[j].colour == ls[mod1(j + 1, n)].colour && continue
        return false
    end
    return true
end

"""
    zamo_regions(g) -> Vector{Int}

The **Zamo regions**: the region numbers of [`zamo_triangles`](@ref) whose
corners continue the pattern outward. A triangle drops out as soon as, at one
of its corners, its two differently-coloured edges meet at a `[1,1,3,3]`-like
`:mixed` node — there the pattern is P1-merged and does not continue.

**Counted per triangle, not per connected component.** "At least 7 Zamo
regions" refers to this count: a single Zamo side has 6, at 7 or more there
is more than one cluster. The components are given by
[`zamo_region_components`](@ref).

The measurement table (fixtures 6, `Zamo(1,8)` end term 8 → 6) is in the
block comment above this function.
"""
zamo_regions(g::CircularGraph) =
    Int[R for (R, ls) in zamo_triangles(g) if _fzr_continues_outward(g, ls)]

zamo_regions(g::WordGraph) = zamo_regions(circular(g))
zamo_regions(m::MorphismGraph) = zamo_regions(circular(m.graph))
zamo_regions(m::CircularMorphismGraph) = zamo_regions(m.graph)

"""
    zamo_region_components(g::CircularGraph) -> Vector{Vector{Int}}

The [`zamo_regions`](@ref), grouped into connected components: two Zamo
regions lie in the same component if they share an edge. A single Zamo side
is exactly ONE component of 6 triangles.

For diagnostics/notebooks only — the direction choice depends on
`length(zamo_regions(g))`.
"""
function zamo_region_components(g::CircularGraph)
    Rs = zamo_regions(g)
    isempty(Rs) && return Vector{Int}[]
    pos = Dict(R => k for (k, R) in enumerate(Rs))
    adj, _ = circular_region_adjacency(g)
    parent = collect(1:length(Rs))
    find(x) = parent[x] == x ? x : (parent[x] = find(parent[x]))
    for R in Rs, (R2, _) in adj[R]
        haskey(pos, R2) || continue
        a, b = find(pos[R]), find(pos[R2])
        a == b || (parent[a] = b)
    end
    groups = Dict{Int,Vector{Int}}()
    for (k, R) in enumerate(Rs)
        push!(get!(groups, find(k), Int[]), R)
    end
    return sort([sort(v) for v in values(groups)], by = first)
end

# ---- the direction choice ---------------------------------------------------
#
# Policy: if Z1 or Z2 can be applied, do so. If there are at least 7 Zamo
# regions and neither applies, try Z1⁻¹, then Z2⁻¹. If a different rule then
# fires afterwards, things get simpler; if only Z1/Z2 fire afterwards, that is
# accepted as-is.
#
# So: (a) forward (`circular_zamo_region_step`, both pairs); (b) if
# nothing fires AND `length(zamo_regions(g)) >= 7`, then Z1⁻¹, then Z2⁻¹;
# (c) a backward result is ACCEPTED if a non-Zamo rule fires on EVERY one of
# its terms (`CIRCULAR_RULES` incl. the generic preimage stage, P1, 2parallel, D4)
# — otherwise it is discarded, the term stays a fixed point, and the case is
# recorded in `zamo_inverse_log()` ("only Z1/Z2 fire afterwards").
#
# WHY "every term": the backward image carries, by construction, the FORWARD
# pattern (Z1⁻¹ gives `rev(Zamo(8,1))`, whose region pattern is exactly the
# LHS of Z1). If the term ran unguarded into the recursion, Z1 would fire
# straight back forward, the original would again have ≥ 7 regions, and
# backward would fire again — an infinite loop. Hence two guards together:
#
#   1. here: accepted only if a non-Zamo rule fires (termination is then
#      carried by that rule's measures: circular_weight, edge count, Σ dist);
#   2. in `reduce_to_circular_leave` (`after_inverse = true` for the recursive call
#      on the backward result): in the FIRST round of that call, NO Zamo step
#      fires (neither forward nor backward) — the non-Zamo rule from guard 1
#      is guaranteed to go first. After that the term has changed, and from
#      round 2 the normal direction choice applies again.

"""
    ZamoInverseLogEntry

A log entry of the direction choice (see [`zamo_inverse_log`](@ref)):
`pair` the backward pair tried, `status ∈ (:accepted, :rejected)`
(accepted/rejected), `word` the term's boundary word, `nodes`/`zamo_regions`
its counts, `fires` per term of the backward result the first non-Zamo rule
that fires (`:circular_rules`, `:p1`, `:two_parallel`, `:d4`) or `nothing`.
"""
struct ZamoInverseLogEntry
    pair::Tuple{Int,Int}
    status::Symbol
    word::Vector{Int}
    nodes::Int
    zamo_regions::Int
    fires::Vector{Union{Nothing,Symbol}}
end

Base.show(io::IO, e::ZamoInverseLogEntry) = print(io,
    "Z", e.pair == (1, 8) ? "1" : "2", "⁻¹ ", e.status, " · boundary word ",
    join(e.word), " · ", e.nodes, " nodes · ", e.zamo_regions,
    " Zamo regions · fires afterwards: ", e.fires)

const _FZR_INVERSE_LOG = ZamoInverseLogEntry[]

"""
    zamo_inverse_log() -> Vector{ZamoInverseLogEntry}
    zamo_inverse_log_clear!()

The LOG of the direction choice: every backward attempt from
`reduce_to_circular_leave`, accepted or rejected. The rejected ones
(`status == :rejected`) are exactly the "only Z1/Z2 fire afterwards" cases.
A caller running many pairs reads and clears the log per pair.
"""
zamo_inverse_log() = copy(_FZR_INVERSE_LOG)
zamo_inverse_log_clear!() = (empty!(_FZR_INVERSE_LOG); nothing)

"""
    zamo_nonzamo_rule_fires(dd::CircularDecorated, cut1, cut2) -> Union{Nothing,Symbol}

Which NON-Zamo rule fires on the term `dd` (morphism with cuts `cut1`/`cut2`)
— the check order is the driver's: `:circular_rules` (`_circular_first_rule_match`, i.e.
`CIRCULAR_RULES` including the generic preimage stage), `:p1`
(`circular_parallel_merge_step`), `:two_parallel` (`circular_2parallel_step`), `:d4`
(`find_circular_d4_any`); `nothing` if none does. This is the acceptance criterion
of the direction choice (section above).
"""
function zamo_nonzamo_rule_fires(dd::CircularDecorated, cut1::Int, cut2::Int)
    _circular_first_rule_match(dd, CIRCULAR_RULES) === nothing || return :circular_rules
    m = CircularMorphismGraph(dd.graph, cut1, cut2)
    fdm = CircularDecoratedMorphism(m, dd.region_labels, dd.outer_label)
    circular_parallel_merge_step(fdm) === nothing || return :p1
    circular_2parallel_step(fdm) === nothing || return :two_parallel
    find_circular_d4_any(m) === nothing || return :d4
    return nothing
end

"""
    circular_zamo_direction_step(g::CircularGraph, cut1, cut2; log = true)
        -> Union{Nothing, Tuple{CircularComboR, Bool}}

THE DIRECTION CHOICE (section above): first forward Z1/Z2
(`circular_zamo_region_step` for `(1,8)`, then `(2,9)`) — result `(combo, false)`;
if nothing fires and `g` has at least 7 [`zamo_regions`](@ref), then Z1⁻¹,
then Z2⁻¹, each accepted only if a non-Zamo rule fires on EVERY term of the
result ([`zamo_nonzamo_rule_fires`](@ref)) — result `(combo, true)`. Otherwise
`nothing`. Every backward attempt is recorded in
[`zamo_inverse_log`](@ref) (`log = false` disables that).

The backward attempt runs ONLY at ≥ 7 Zamo regions — a single Zamo side has
6, only from 7 does more than one cluster exist.
"""
function circular_zamo_direction_step(g::CircularGraph, cut1::Int, cut2::Int; log::Bool = true)
    for pair in ((1, 8), (2, 9))
        out = circular_zamo_region_step(g; pair = pair)
        out === nothing || return (out, false)
    end
    nreg = length(zamo_regions(g))
    nreg >= 7 || return nothing
    for pair in ((1, 8), (2, 9))
        out = circular_zamo_region_step(g; pair = pair, inverse = true)
        out === nothing && continue
        fires = Union{Nothing,Symbol}[zamo_nonzamo_rule_fires(dd, cut1, cut2)
                                      for (dd, _) in pairs_of(out)]
        ok = !isempty(fires) && all(!isnothing, fires)
        log && push!(_FZR_INVERSE_LOG, ZamoInverseLogEntry(
            pair, ok ? :accepted : :rejected, collect(letters(g.word)),
            length(g.nodes), nreg, fires))
        ok && return (out, true)
    end
    return nothing
end

# ---- diagnostics -------------------------------------------------------------

"""
    zamo_region_signature(g; pair-less) -> Vector

The signature of `g`'s region pattern, for by-hand and test comparisons: per
inner region the **lexicographically smallest rotation** of its colour
sequence, plus the multiset of edge colours between inner regions. This is
the coarse sieve — it does not fully distinguish the
four Zamo sides (LHS(1,8) and RHS(1,8) have the same shape); the isomorphism
check in `circular_zamo_region_matches` handles that.
"""
function zamo_region_signature(g::CircularGraph)
    inner = _fzr_inner(g)
    iset = Set(R for (R, _) in inner)
    words = String[]
    for (_, ls) in inner
        n = length(ls)
        cs = [l.colour for l in ls]
        push!(words, minimum(join(cs[mod1(j + r, n)] for j in 1:n) for r in 0:(n - 1)))
    end
    adj, _ = circular_region_adjacency(g)
    cols = Int[]
    for (R, _) in inner, (R2, ei) in adj[R]
        (R2 in iset && R2 > R) && push!(cols, g.edges[ei].colour)
    end
    return (boundary_words = sort(words), edge_colours = sort(cols))
end

zamo_region_signature(g::WordGraph) = zamo_region_signature(circular(g))
zamo_region_signature(m::MorphismGraph) = zamo_region_signature(circular(m.graph))
