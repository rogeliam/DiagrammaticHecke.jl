# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/Combo.jl  —  formal linear combinations of diagrams + a rewrite scaffold
#
# GOAL. On a FIXED boundary word we want to form formal linear
# combinations of the possible diagrams on it, and later rewrite/simplify them by
# LOCAL rules. This file builds the CONSTRUCT ONLY — the concrete rewrite relations
# are supplied separately.
#
# DESIGN
#   * Coefficients: **ℤ** (`Int`). The type is parametrised `DiagramCombo{T}`
#     so swapping to `R = K[α_s]` later costs little (need `zero(T)`, `+`, `*`,
#     `==` on T). Default constructor uses `Int`.
#   * Diagram = `WordGraph` (the planar diagram). THE KEY REQUIREMENT is a
#     canonical form so ISOMORPHIC diagrams collapse to one term in a combination.
#     `canonical_key(g)` is that hook; see its caveats.
#   * A combination is a SPARSE map canonical_key ⇒ (representative, coeff), with
#     zero coefficients dropped. Two diagrams with the same key are "the same term".
#   * Rules: `DiagramRule` = a local match→combination replacement. `simplify`
#     applies rules to a fixpoint. NO concrete rules are defined here — they are
#     supplied separately. `simplify` with an empty rule set just normalises terms.

# ---- canonical key: collapse isomorphic diagrams ---------------------------
#
# The key is a STRUCTURAL fingerprint: the boundary word's normal form, plus the
# sorted node fingerprints, plus the sorted edge fingerprints. Edge fingerprints
# use a CANONICAL node numbering — a function of the graph structure only, not of
# array order — computed by `_dart_canonical_key` (diagram/CanonicalKey.jl, shared
# with `circular_canonical_key`). That numbering is what makes two structurally
# identical nodes at swapped array positions one and the same term.
#
# ROTATION. The key takes the raw slot NUMBER, so two diagrams that differ by
# rotating a node are different terms. That follows from the data model: a
# `WordGraph` stores no rotation system (diagram/Graph.jl §REP), and `_slot_colour`
# reads a node's colour off the slot PARITY (§PLANAR SLOT CONVENTION), so the slot
# number carries meaning that a rotation-invariant key would have to throw away.
# Making the key rotation-invariant means changing that convention, not the key:
# `_port_fp`, the numbering and `_node_fp` would all have to normalise cyclically
# together, and `colours[1]` — which `_slot_colour`, `opposite_slot`,
# `is_planar_braid` and the rules read — would lose its meaning.
#
# `_slot_bases` and `_min_edge_fp` below are the groundwork for such a convention
# (tested, currently unused); `slot_rotation_keys` answers rotation-equality
# directly for callers that need it.
#
# An edge is a colour + an unordered pair of port fingerprints. `renumber` maps the
# ORIGINAL node array index to its CANONICAL rank; Leaf/Circle ports already carry
# an intrinsic, index-independent identity and ignore it.
_port_fp(p::Leaf, renumber, base) = (0, p.k, 0)
_port_fp(p::NodePort, renumber, base) =
    (1, renumber[p.node], base === nothing ? p.slot : base[p.node][p.slot])
_port_fp(p::Circle, renumber, base) = (2, p.colour, 0)      # free-circle marker

function _edge_fp(e::Edge, renumber, base = nothing)
    a, b = _port_fp(e.a, renumber, base), _port_fp(e.b, renumber, base)
    lo, hi = a <= b ? (a, b) : (b, a)
    return (e.colour, lo, hi)
end

# GROUNDWORK, UNUSED (see the block above). The candidate slot origins of
# a node: all rotations `r` making the arm sequence minimal. Usually exactly one; for
# a symmetric arm sequence (e.g. `[s,t,s,t]` of an m=2 crossing) several, and then
# minimising the whole key decides. Returns a vector of slot↦rank tables (1-based).
function _slot_bases(nd::Node)
    d = _node_degree(nd)
    arms = [_slot_colour(nd, k) for k in 1:d]
    rots = [[arms[mod1(k + r, d)] for k in 1:d] for r in 0:d-1]
    best = minimum(rots)
    return [[mod1(k - r, d) for k in 1:d] for r in 0:d-1 if rots[r+1] == best]
end

# GROUNDWORK, UNUSED (see the block above). The edge fingerprints minimised
# over all combinations of slot origins (`_slot_bases`) — the part of the key that
# divides out the slot numbering; `renumber` (the node numbering) stays fixed.
function _min_edge_fp(g::WordGraph, renumber)
    best = nothing
    for combo in Iterators.product((_slot_bases(nd) for nd in g.nodes)...)
        base = collect(combo)
        t = Tuple(sort([_edge_fp(e, renumber, base) for e in g.edges]))
        (best === nothing || t < best) && (best = t)
    end
    return best
end

_node_fp(nd::Node) = (nd.kind, Tuple(nd.colours), nd.m)

# ---- canonical node numbering (dart walk) -----------------------------------
#
# The canonical-numbering machinery lives in the generic core
# diagram/CanonicalKey.jl, shared with `circular_canonical_key`
# (circular/CircularGraph.jl). `canonical_key` below supplies `WordGraph`'s specifics: the
# node fingerprint `_node_fp`, the degree `_node_degree`, and `free_origin = false`
# (no slot-origin minimisation — the open rotation-invariance gap above).

"""
    canonical_key(g::WordGraph)

A hashable key that is EQUAL for diagrams we consider the same, so that isomorphic
diagrams collapse into one term of a linear combination. A structural fingerprint
(boundary normal form + sorted node fingerprints + sorted edge fingerprints, the
latter using a CANONICAL node numbering — see the note at the top of this file —
rather than raw array indices, so the key is independent of how the nodes happen
to be ordered in the array).

The key MISSES an isomorphism that needs a planar rotation system: a `WordGraph`
stores none (diagram/Graph.jl §REP), so diagrams that agree only up to rotating a
node are different keys. See the ROTATION note at the top of this file.
"""
function canonical_key(g::WordGraph)
    bnd = letters(g.word)                              # already in cyclic normal form
    # `free_origin = false`: `Node` has a fixed slot layout, slot numbers enter the
    # key RAW — the dart walk only replaces the node numbering, not the slot
    # semantics, so the rotation-invariance gap above stays open.
    # Head of the flat Int certificate: kind id + colour count + colours + m.
    head! = (out, nd, _o) -> begin
        push!(out, _ck_kind_id(nd.kind), length(nd.colours))
        append!(out, nd.colours)
        push!(out, nd.m)
    end
    key = _dart_canonical_key(g.nodes, g.edges, bnd, _node_fp, _node_degree,
                              head!, false)
    return key
end

"""
    slot_rotation_keys(g::WordGraph) -> Set

The `canonical_key`s of ALL renumberings of `g` obtained by cyclically rotating the
slots of each node. Two diagrams are equal *modulo slot rotation* exactly when these
sets intersect:

    !isempty(intersect(slot_rotation_keys(g), slot_rotation_keys(h)))

Purpose: `canonical_key` carries the absolute slot NUMBER, although only the cyclic
ARM ORDER at a node has meaning (see the rotation-invariance note above). Diagrams
that are the same therefore count as different — the defect that showed up as the
mirrored embedding in `_r8_trivalent_past_crossing`.

The set has `prod(deg(n))` elements — usable for small diagrams (fixtures, tests),
NOT a replacement for `canonical_key` as the `DiagramCombo` key.
"""
function slot_rotation_keys(g::WordGraph)
    degs = [_node_degree(nd) for nd in g.nodes]
    keys = Set{Any}()
    for shifts in Iterators.product((0:d-1 for d in degs)...)
        rot(p) = p isa NodePort ?
                 NodePort(p.node, mod1(p.slot + shifts[p.node], degs[p.node])) : p
        h = WordGraph(g.word, g.nodes,
                      Edge[Edge(e.colour, rot(e.a), rot(e.b)) for e in g.edges])
        push!(keys, canonical_key(h))
    end
    return keys
end

# ---- the linear combination ------------------------------------------------

"""
    DiagramCombo{T}

A formal `T`-linear combination of diagrams (`WordGraph`s), stored sparsely as
`canonical_key ⇒ (representative_diagram, coefficient)` with zero coefficients
dropped. `T` is the coefficient type (default `Int`; later `R = K[α_s]`). Terms with
the same `canonical_key` are merged — this is where isomorphic diagrams combine.
"""
struct DiagramCombo{T}
    terms::Dict{Any,Tuple{WordGraph,T}}
end

DiagramCombo{T}() where {T} = DiagramCombo{T}(Dict{Any,Tuple{WordGraph,T}}())
DiagramCombo() = DiagramCombo{Int}()

"A single diagram as a combination with coefficient `one(T)`."
function DiagramCombo{T}(g::WordGraph) where {T}
    c = DiagramCombo{T}()
    c.terms[canonical_key(g)] = (g, one(T))
    return c
end
DiagramCombo(g::WordGraph) = DiagramCombo{Int}(g)

Base.length(c::DiagramCombo) = length(c.terms)
Base.isempty(c::DiagramCombo) = isempty(c.terms)
coefficient(c::DiagramCombo{T}, g::WordGraph) where {T} =
    haskey(c.terms, canonical_key(g)) ? c.terms[canonical_key(g)][2] : zero(T)

"The (representative, coefficient) pairs, zero coefficients already excluded."
pairs_of(c::DiagramCombo) = [(v[1], v[2]) for v in values(c.terms)]

# add `coeff * g` into `c` in place, merging isomorphic terms and dropping zeros.
function _add!(c::DiagramCombo{T}, g::WordGraph, coeff::T) where {T}
    coeff == zero(T) && return c
    k = canonical_key(g)
    if haskey(c.terms, k)
        rep, old = c.terms[k]
        new = old + coeff
        new == zero(T) ? delete!(c.terms, k) : (c.terms[k] = (rep, new))
    else
        c.terms[k] = (g, coeff)
    end
    return c
end

function Base.:+(a::DiagramCombo{T}, b::DiagramCombo{T}) where {T}
    out = DiagramCombo{T}(copy(a.terms))
    for (g, coeff) in pairs_of(b)
        _add!(out, g, coeff)
    end
    return out
end

Base.:-(a::DiagramCombo{T}) where {T} = (zero(T) - one(T)) * a
Base.:-(a::DiagramCombo{T}, b::DiagramCombo{T}) where {T} = a + (-b)

function Base.:*(s::T, a::DiagramCombo{T}) where {T}
    out = DiagramCombo{T}()
    s == zero(T) && return out
    for (g, coeff) in pairs_of(a)
        _add!(out, g, s * coeff)
    end
    return out
end
Base.:*(a::DiagramCombo{T}, s::T) where {T} = s * a

Base.:(==)(a::DiagramCombo, b::DiagramCombo) =
    Set(keys(a.terms)) == Set(keys(b.terms)) &&
    all(a.terms[k][2] == b.terms[k][2] for k in keys(a.terms))

function Base.show(io::IO, c::DiagramCombo{T}) where {T}
    if isempty(c)
        print(io, "0 (empty DiagramCombo{$T})"); return
    end
    ps = sort(pairs_of(c); by = p -> letters(p[1].word))
    print(io, "DiagramCombo{$T} with ", length(ps), " term(s):")
    for (g, coeff) in ps
        print(io, "\n  ", coeff, " · ", g)
    end
end

# ---- rewrite rules (scaffold; concrete relations added separately) ---

"""
    DiagramRule

A local rewrite: `apply(g) -> Union{Nothing, DiagramCombo{T}}` returns the linear
combination a single diagram rewrites to (e.g. "needle = 0", "two dots = α_s · id"),
or `nothing` if the rule does not fire on `g`. These are supplied separately; this file
only defines the interface and the fixpoint driver.
"""
struct DiagramRule{T}
    name::Symbol
    apply::Function          # WordGraph -> Union{Nothing, DiagramCombo{T}}
end

"""
    simplify(c::DiagramCombo{T}, rules; maxpasses = 1000) -> DiagramCombo{T}

Apply `rules` to every term repeatedly until nothing changes (or `maxpasses`).
A rule that fires on a term replaces that term (× its coefficient) by the rule's
output combination; results are re-merged, so isomorphic diagrams collapse. With an
empty rule set this just returns `c` (terms already merged on construction).

CONJECTURE: local rules suffice to reach a normal form. Not proven — if this
loops or stalls on some input, that input is exactly what to inspect by hand.
"""
function simplify(c::DiagramCombo{T}, rules::AbstractVector{<:DiagramRule};
                  maxpasses::Int = 1000) where {T}
    isempty(rules) && return c
    cur = c
    for _ in 1:maxpasses
        nxt = DiagramCombo{T}()
        changed = false
        for (g, coeff) in pairs_of(cur)
            fired = false
            for r in rules
                res = r.apply(g)
                if res !== nothing
                    nxt = nxt + (coeff * res)
                    fired = true; changed = true
                    break
                end
            end
            fired || _add!(nxt, g, coeff)
        end
        cur = nxt
        changed || return cur
    end
    return cur
end
