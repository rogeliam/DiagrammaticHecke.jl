# circular/CircularGraphKey.jl — rotation canonicalisation of a CircularNode and the
# CircularGraph canonical key (circular_canonical_key, ==, hash). Included after
# CircularGraph.jl.
#
# ---- rotation canon — the core --------------------------------------------------
#
# Equality holds ONLY up to rotation, NOT up to reflection: the arm sequence is
# never reversed, so the diagram's orientation is preserved (consistent with
# the clockwise slot convention and with `flip`). `_min_rotation` is the same
# notion as `_canonical_rotation` in words/CircularWord.jl, but for arbitrary
# Int vectors — that code is left untouched, this is a local reimplementation
# (an O(n²) scan is fine, node degrees stay small).

"The lexicographically smallest cyclic rotation of `v`. `v` is not mutated."
function _min_rotation(v::Vector{Int})
    n = length(v)
    n <= 1 && return copy(v)
    best = copy(v)
    for k in 1:(n - 1)
        rot = vcat(v[k+1:end], v[1:k])
        rot < best && (best = rot)
    end
    return best
end

"""
    _min_rotation_shift(v) -> Int

The offset `k` (`0 ≤ k < n`) for which `vcat(v[k+1:end], v[1:k])` is the
smallest rotation. Old slot `s` thereby lands on new slot `mod1(s - k, n)`.

⚠️ Returns only the FIRST such offset. If the arm sequence itself is
rotation-symmetric (`[1,3,1,3]` under rotation by 2, `[1,1,1]` under any),
there are SEVERAL equally good ones — then the choice is arbitrary and the
key depends on it. Using this instead of `_min_rotation_shifts` lets
rotations of the same diagram change the key. Anyone building the key must
use `_min_rotation_shifts` (plural) and minimise over ALL of them.
"""
function _min_rotation_shift(v::Vector{Int})
    n = length(v)
    n <= 1 && return 0
    best, bestk = copy(v), 0
    for k in 1:(n - 1)
        rot = vcat(v[k+1:end], v[1:k])
        if rot < best
            best, bestk = rot, k
        end
    end
    return bestk
end

"""
    _min_rotation_shifts(v) -> Vector{Int}

ALL offsets `k` for which `vcat(v[k+1:end], v[1:k])` is the lexicographically
smallest rotation — the stabiliser orbit of the arm sequence. Usually a single
one; several for a symmetric arm sequence (`[1,3,1,3]` ⇒ `[0,2]`, `[1,1,1]` ⇒
`[0,1,2]`).

The key must minimise over all of them, or an arbitrary choice ends up
deciding it — see `_min_rotation_shift` and `circular_canonical_key`.
"""
function _min_rotation_shifts(v::Vector{Int})
    n = length(v)
    n <= 1 && return [0]
    rots = [vcat(v[k+1:end], v[1:k]) for k in 0:(n - 1)]
    best = minimum(rots)
    return [k for k in 0:(n - 1) if rots[k + 1] == best]
end

"Canonical key of a node: kind + lexicographically smallest rotation."
circular_key(nd::CircularNode) = (nd.kind, Tuple(_min_rotation(nd.arms)))

# Fast path: identical arm sequence (the common case, e.g. right after the
# converter) skips the O(n²) scan.
function Base.:(==)(a::CircularNode, b::CircularNode)
    a.arms == b.arms && return true
    return circular_key(a) == circular_key(b)
end
Base.hash(nd::CircularNode, h::UInt) = hash(circular_key(nd), h)

# ---- graph key --------------------------------------------------------------
#
# THREE degrees of freedom must be normalised away, or equal diagrams get
# different keys:
#
# 1. SLOT ROTATION. `CircularNode` does not distinguish a slot 1, but a
#    `NodePort(i, s)` names a concrete slot. So every slot is converted to its
#    node's canonical numbering before fingerprinting: if `k` is the offset
#    from `_min_rotation_shift`, then `rot[j] = arms[mod1(j+k, deg)]`, so old
#    slot `s` lands on `mod1(s - k, deg)`.
#
# 2. NODE INDEX. Exactly as with `canonical_key` (diagram/Combo.jl), the raw
#    array indices must not enter the key: two structurally identical nodes at
#    swapped positions would otherwise give different keys for isomorphic
#    diagrams. Weisfeiler–Leman avoids this there; we reuse the same
#    technique via `_circular_canonical_numbering`.
#
# 3. LEAF ROTATION. `g.word` is a `CircularWord` — by definition only fixed up
#    to cyclic rotation, see words/CircularWord.jl. `Leaf(k)`, however, names
#    an ABSOLUTE position in the word. Two `CircularGraph`s differing only in
#    where the boundary leaf numbering starts are therefore the SAME diagram,
#    but would otherwise get different keys because `p.k` enters the
#    fingerprint raw. The key must therefore be minimised over all rotations
#    `r ∈ 0:n-1` that rename `Leaf(k)` to `Leaf(mod1(k-r, n))` — BUT ONLY over
#    those that leave `letters(g.word)` unchanged (the actual symmetries of the
#    boundary word): a rotation that changes the letter sequence would
#    otherwise wrongly identify GENUINELY different diagrams. `n == 0` has no
#    leaves and needs no rotation.

"""
    _circular_leaf_rotations(w::CircularWord) -> Vector{Int}

The leaf rotations `r` (`0 ≤ r < n`) for which `Leaf(k) -> Leaf(mod1(k-r, n))`
is a genuine symmetry of `w`, i.e. `letters(w)` shifted cyclically by `r`
gives the same sequence. `r = 0` (identity) is always included. Empty word ⇒
`[0]` (nothing to rotate). Unused in `circular_canonical_key` — the
word there is fixed, see below — kept as a building block for a possible
future redesign.
"""
function _circular_leaf_rotations(w::CircularWord)
    lets = letters(w)
    n = length(lets)
    n == 0 && return [0]
    return [r for r in 0:(n - 1) if all(j -> lets[mod1(j + r, n)] == lets[j], 1:n)]
end

# ---- canonicalisation lives in the generic core -------------------------------
#
# The actual dart-walk canonicalisation (colour-free, see the header of
# diagram/CanonicalKey.jl) lives in the GENERIC core diagram/CanonicalKey.jl —
# shared with `canonical_key` (diagram/Combo.jl). `circular_canonical_key`
# below calls it with `CircularGraph`'s specifics: `circular_key` as the node
# fingerprint (already rotation-normalised, see above), `_circular_degree`
# (arm count) as the degree, `_min_rotation_shifts` as the candidates for the
# slot origin (the GENUINE combinatorial search that `canonical_key`'s
# analogue does NOT make — see CanonicalKey.jl point 2), and
# `mod1(slot - shift, deg)` as the canonical slot.
"""
    circular_canonical_key(g::CircularGraph)

Hashable key, invariant under ALL THREE degrees of freedom of the type:
slot rotation at each node, the order of nodes in the array (via a dart
walk, like `canonical_key` in diagram/Combo.jl), and leaf rotation (`g.word`
is only defined up to rotation). Consists of the boundary word in normal
form, the sorted node keys, and the sorted edge fingerprints — minimised
over all equally good slot origins.

THE WORD IS FIXED: no minimisation over `_circular_leaf_rotations(g.word)` —
the left marking is part of the identity.

ALL equally good slot origins per node are considered, not just one
arbitrarily chosen offset (`_min_rotation_shift`). For a rotation-symmetric
arm sequence (`[1,3,1,3]`, `[1,1,1]`) there are several, and leaving the
choice arbitrary would let it leak into the key — a rotation of the diagram
could then produce a different key even though the diagram is the same.
"""
function circular_canonical_key(g::CircularGraph)
    # Memoised in the lazy field `_ckey`: `==`/`hash` run through here, and
    # without the cache the key was recomputed on every `_add!` addition and
    # every memo lookup — on a dict hit, even for both sides.
    k = g._ckey[]
    k === nothing || return k
    # `free_origin = true`: each node's arm origin is FREE and fixed by the
    # dart traversal, not guessed via `_min_rotation_shifts` and minimised
    # over all combinations (see the header of diagram/CanonicalKey.jl for
    # the cost of that approach). A node's head in the certificate is
    # kind id + arm count + arm sequence STARTING AT THE ORIGIN (flat Int
    # certificate, see `_dart_certificate`).
    head! = (out, nd, o) -> begin
        d = length(nd.arms)
        push!(out, _ck_kind_id(nd.kind), d)
        for t in 1:d
            push!(out, nd.arms[mod1(o + t - 1, d)])
        end
    end
    k = _dart_canonical_key(g.nodes, g.edges, letters(g.word), circular_key, _circular_degree,
                            head!, true)
    g._ckey[] = k
    return k
end

Base.:(==)(a::CircularGraph, b::CircularGraph) = circular_canonical_key(a) == circular_canonical_key(b)
Base.hash(g::CircularGraph, h::UInt) = hash(circular_canonical_key(g), h)

"""
    _circular_leaf_colour(g::CircularGraph, k::Int) -> Int

Colour of the edge at leaf `k` (0 if none). Analogue of `leaf_colour`
(morphism/MorphismGraph.jl), but `CircularGraph`-typed. Lives here, not
render-local in CircularTutte.jl, because `circular/rules/CircularRules.jl`
(§C8, `_fr_braid_relation`) needs it too.
"""
function _circular_leaf_colour(g::CircularGraph, k::Int)
    # Byte-identical to `leaf_colour` above, under a separate name that code
    # outside this file group depends on.
    for e in g.edges
        e.a isa Leaf && e.a.k == k && return e.colour
        e.b isa Leaf && e.b.k == k && return e.colour
    end
    return 0
end

