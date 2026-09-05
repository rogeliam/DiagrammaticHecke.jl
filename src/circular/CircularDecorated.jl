# circular/CircularDecorated.jl — decorated CircularGraphs.
#
# Modeled on `DecoratedDiagram` (diagram/Decorated.jl): a CircularGraph together
# with a polynomial per inner cell. The cell infrastructure already exists
# (CircularFaces.jl: face_count, inner_faces, face_of_port, gap_cell, sector_cell) —
# here only the label vector and its linear combination `CircularComboR` get added
# (CircularCombo.jl already defines `CircularCombo{T}`, but `CircularComboR = CircularCombo{
# SoergelPoly}` needs `CircularDecorated` as its element type — analogous to how
# `DiagramComboR` is NOT `DiagramCombo{SoergelPoly}`, but runs over
# `DecoratedDiagram`).
#
# THE OUTER LABEL — the circular side; the plain side lives in
# diagram/Decorated.jl, with the derivation there.
# `regions(g)` only counts INNER regions; the outer one has no id
# (`cell_of_face[outer_face] = 0`, `_circular_dart_region` skips it). Treating a
# polynomial that lands there as the term's COEFFICIENT is only correct in the
# distance-0 case (there the polynomial really is a global scalar, see
# `circular_extract_scalars`); in general it must stay in place outside. So
# `CircularDecorated` carries the outer label in its OWN field, so
# `region_labels` stays indexed by `R.id` with no shift.
# `==`/`hash`/`circular_canonical_key` include it: `key(α·D) != key(D)` for an
# α sitting in place outside.
# The distance-0 case applies to EVERY marked morphism
# (`circular_outer_region_distance` is constantly 0, the boundary circle is a
# wall). The field remains necessary anyway: an UNMARKED `CircularDecorated`
# has no distance, and that is exactly where the label keeps sitting in place.

"""
    CircularDecorated

A `CircularGraph` together with a vector `region_labels` of polynomials
(`SoergelPoly`), **one per REGION** (order of `regions`, index = `R.id`),
**plus the label of the OUTER region** (`outer_label`, default `1`). Label `1`
means unlabeled.

The third argument is optional:
`CircularDecorated(g, labels)` is `CircularDecorated(g, labels, one(R))`. See the file
header for why.

**Labels attach to REGIONS, not faces**: regions
are the FINER decomposition (`region_count ≥ face_count`) and the same notion
`find_circular_d4_match` decides on. A single
diagram can have one face but several regions, and a label naming "the"
region of that face would be ambiguous — see the `Region` docstring
(diagram/Faces.jl) for the general distinction.
"""
struct CircularDecorated
    graph::CircularGraph
    region_labels::Vector{SoergelPoly}
    outer_label::SoergelPoly
    # Lazy cache for `circular_canonical_key(::CircularDecorated)` — same rationale as
    # the `_ckey` field of `CircularGraph`: a pure function of the other fields, not
    # part of `==`/`hash`/`show`.
    _ckey::Base.RefValue{Any}

    function CircularDecorated(g::CircularGraph, labels::Vector{SoergelPoly},
                          outer::SoergelPoly = one(SoergelPoly))
        length(labels) == region_count(g) || throw(ArgumentError(
            "expected $(region_count(g)) labels (one per region), got $(length(labels))"))
        return new(g, labels, outer, Ref{Any}(nothing))
    end
end

"""
    circular_decorated(g::CircularGraph) -> CircularDecorated
    circular_decorated(g::CircularGraph, labels::Vector{SoergelPoly}[, outer]) -> CircularDecorated

`g` with region labels: `1` everywhere if not given. `outer` is the label of
the OUTER region (default `1`).
"""
circular_decorated(g::CircularGraph) = CircularDecorated(g, fill(one(SoergelPoly), region_count(g)))
circular_decorated(g::CircularGraph, labels::Vector{SoergelPoly},
              outer::SoergelPoly = one(SoergelPoly)) = CircularDecorated(g, labels, outer)

circular_boundary(d::CircularDecorated) = circular_boundary(d.graph)
Base.length(d::CircularDecorated) = length(d.graph)
face_count(d::CircularDecorated) = face_count(d.graph)
inner_faces(d::CircularDecorated) = inner_faces(d.graph)
region_count(d::CircularDecorated) = region_count(d.graph)
regions(d::CircularDecorated) = regions(d.graph)

"Label of region `i` (index = `R.id` from `regions`)."
region_label(d::CircularDecorated, i::Int) = d.region_labels[i]

"""
    outer_label(d::CircularDecorated) -> SoergelPoly

The label of the OUTER region (see file header). `1` means
unlabeled. The outer region has no region id — it is therefore NOT reachable
via `region_label(d, i)`.
"""
outer_label(d::CircularDecorated) = d.outer_label

# If any region carries the label 0, the whole diagram is 0.
# The method for `CircularGraph` (always `false`) lives in circular/CircularCombo.jl, where
# `_add!` also calls it; the docstring there explains origin and effect.
# The OUTER region counts too — a `0` outside
# annuls the term just like one inside.
_is_zero_diagram(d::CircularDecorated) = any(iszero, d.region_labels) || iszero(d.outer_label)

function Base.:(==)(a::CircularDecorated, b::CircularDecorated)
    a.region_labels == b.region_labels || return false
    a.outer_label == b.outer_label || return false
    return a.graph == b.graph
end

Base.hash(d::CircularDecorated, h::UInt) =
    hash((circular_canonical_key(d.graph), d.region_labels, d.outer_label),
         hash(:CircularDecorated, h))

"""
    circular_canonical_key(d::CircularDecorated)

Key for `CircularComboR`: the structural `circular_canonical_key` of the graph, the
label sequence in region order, **and the outer label**.

The outer label is part of this: if
an `α` sits in place outside, `α·D` is a DIFFERENT term from `D`.
"""
function circular_canonical_key(d::CircularDecorated)
    k = d._ckey[]
    k === nothing || return k
    k = (circular_canonical_key(d.graph), Tuple(d.region_labels), d.outer_label)
    d._ckey[] = k
    return k
end

function Base.show(io::IO, d::CircularDecorated)
    print(io, "CircularDecorated(", d.graph)
    for i in 1:region_count(d)
        f = d.region_labels[i]
        isone(f) && continue
        print(io, ", R", i, " ↦ ", f)
    end
    isone(d.outer_label) || print(io, ", [outer] ↦ ", d.outer_label)
    print(io, ")")
end

# ---- CircularComboR over CircularDecorated -----------------------------------------------
#
# `CircularCombo{T}` (circular/CircularCombo.jl) is already parametric — here just the
# `CircularDecorated` version of `CircularComboR` (instead of a raw `CircularCombo{SoergelPoly}`
# over `CircularGraph`). Analogous to how `DiagramComboR` runs over
# `DecoratedDiagram`, not over `WordGraph`.

"""
    CircularCombo{T}(d::CircularDecorated) -> CircularCombo{T}

A single decorated CircularGraph as a combination with coefficient `one(T)`, keyed
by `circular_canonical_key(d::CircularDecorated)` (this includes the labels — two
structurally equal but differently decorated graphs are different terms).
"""
function CircularCombo{T}(d::CircularDecorated) where {T}
    c = CircularCombo{T}()
    # A diagram with a 0-labeled region IS 0 — the combination stays empty.
    # See `_is_zero_diagram`.
    _is_zero_diagram(d) && return c
    c.terms[circular_canonical_key(d)] = (d, one(T))
    return c
end

coefficient(c::CircularCombo{T}, d::CircularDecorated) where {T} =
    haskey(c.terms, circular_canonical_key(d)) ? c.terms[circular_canonical_key(d)][2] : zero(T)

function _add!(c::CircularCombo{T}, d::CircularDecorated, coeff::T) where {T}
    coeff == zero(T) && return c
    k = circular_canonical_key(d)
    if haskey(c.terms, k)
        rep, old = c.terms[k]
        new = old + coeff
        new == zero(T) ? delete!(c.terms, k) : (c.terms[k] = (rep, new))
    else
        c.terms[k] = (d, coeff)
    end
    return c
end
_add!(c::CircularComboR, d::CircularDecorated, coeff::Integer) = _add!(c, d, coeff * one(SoergelPoly))

"""
    _circular_nodes_preserved(g_old, g_new) -> Bool

`true` if `g_new` contains the nodes of `g_old` **unchanged and in the same
numbering** (new nodes appended at the end are allowed). Only then may a port
`NodePort(n, s)` be read as the SAME place in both graphs.

Used as a gatekeeper for the anchor-based region mapping
([`_circular_gapless_region_map`](@ref)): `_circular_delete_nodes` renumbers on node
deletion, and after a renumbering a port could accidentally match a different
node. In that case the anchor switches off, returning an empty map.
"""
function _circular_nodes_preserved(g_old::CircularGraph, g_new::CircularGraph)
    length(g_new.nodes) >= length(g_old.nodes) || return false
    return all(i -> arms(g_new.nodes[i]) == arms(g_old.nodes[i]),
               1:length(g_old.nodes))
end

"""
    _circular_gapless_region_map(g_old, g_new) -> Dict{Int,Int}

Maps OLD regions to NEW ones, **without using the boundary** — via an anchor
in the interior: an edge that occurs identically in both graphs (same colour,
same two ports, same orientation). For such an edge, the dart `e.a → e.b` is
the same place in both graphs, so the region left of it in `g_old` belongs to
the one left of it in `g_new`; likewise for the reverse dart.

**Why this anchor**: region addressing must not go through
boundary gaps. A region without a boundary gap (lens, barbell, an enclosed
area) has no boundary gap to address it through — but it does have edges, and
edges outside the rule stay put piece by piece under a LOCAL rewiring. The
anchor is exactly the local description the rule needs.

**Preconditions — and where mapping is withheld:**
- `_circular_nodes_preserved` must hold (no renumbering), otherwise an empty dict.
- The anchor edge must occur in `g_new` **exactly once** (otherwise ambiguous).
- If two anchors for the same old region give different new regions, the old
  region has fallen apart; it is then **not** mapped (the caller throws or
  drops it, depending on the label) — better a clear error than a guess.
"""
function _circular_gapless_region_map(g_old::CircularGraph, g_new::CircularGraph)
    _circular_nodes_preserved(g_old, g_new) || return Dict{Int,Int}()
    dr_old, t_old = _circular_dart_region(g_old)
    dr_new, t_new = _circular_dart_region(g_new)

    cnt = Dict{Tuple{Int,Any,Any},Int}()
    for e in g_new.edges
        k = (e.colour, e.a, e.b)
        cnt[k] = get(cnt, k, 0) + 1
    end

    out  = Dict{Int,Int}()
    bad  = Set{Int}()
    for e in g_old.edges
        get(cnt, (e.colour, e.a, e.b), 0) == 1 || continue
        (haskey(t_old.port_dart, e.a) && haskey(t_new.port_dart, e.a)) || continue
        do1, dn1 = t_old.port_dart[e.a], t_new.port_dart[e.a]
        for (dO, dN) in ((do1, dn1), (t_old.darts[do1].rev, t_new.darts[dn1].rev))
            R = get(dr_old, dO, 0)
            S = get(dr_new, dN, 0)
            (R == 0 || S == 0) && continue
            if haskey(out, R) && out[R] != S
                push!(bad, R)                      # region fell apart ⇒ don't guess
            else
                out[R] = S
            end
        end
    end
    for R in bad
        delete!(out, R)
    end
    return out
end

# ---- label transfer -----------------------------------------------------------

"""
    _circular_transfer_labels(g_old, g_new, labels_old; overrides = Dict{Int,SoergelPoly}(),
                         outer = nothing)
        -> Vector{SoergelPoly}

Transfers the **region** labels of `g_old` onto the graph `g_new` (produced by
node deletion/rewiring): each old region is mapped via its first boundary gap
(`gaps` are boundary-indexed, the boundary word never changes). `overrides`
maps NEW region id → new label and wins over the transfer.

**Regions WITHOUT a boundary gap.** For these there is no boundary gap to map
through; they are mapped via an INNER anchor ([`_circular_gapless_region_map`](@ref))
— an edge that occurs identically in both graphs. Only if that fails too does
the default apply (trivial ⇒ drop, otherwise throw) — **unless `outer`
is set**.

**`outer` — the target slot for anchor-less regions.** A
`Ref{SoergelPoly}`: if a borderless region finds no anchor, its label is
**multiplied into** it instead of throwing. The call site assigns the result
as the `outer_label` of the new `CircularDecorated` (`_fr_barbell`,
`circular_extract_floating_components`) and usually seeds the ref with the
previous outer label (`Ref(fd.outer_label)`), so it is not lost while the
graph is rebuilt. Without `outer` (the default, and hence all other call sites) the
behaviour is unchanged: a non-trivial label with no target is an error, not a
silent loss.

**Regions, not faces**: the transfer runs over regions because `find_circular_d4_match`
decides over regions, not faces — the two notions diverge exactly where a dot
disappears.

**MERGING — the non-trivial label wins.** Regions
are FINER than faces, so two old regions more often land in one new region;
the rule for what happens then:

| old labels | new label |
|---|---|
| both `1` | `1` |
| exactly one ≠ `1` | that one |
| both ≠ `1` | **error** — genuine information would be lost here |

**A NEW region with no predecessor** gets `1`, unless an `overrides` entry
claims it (this is how `apply_circular_d4` sets the `α_i` into the split cell).
"""
function _circular_transfer_labels(g_old::CircularGraph, g_new::CircularGraph,
                              labels_old::Vector{SoergelPoly};
                              overrides::Dict{Int,SoergelPoly} = Dict{Int,SoergelPoly}(),
                              outer::Union{Nothing,Base.RefValue{SoergelPoly}} = nothing)
    rn = regions(g_new)
    labels_new = fill(one(SoergelPoly), length(rn))
    # Which old region has already claimed this new region? (only for the
    # error message — the decision itself depends only on the labels.)
    claimed_by = Dict{Int,Int}()
    ro = regions(g_old)
    # Only compute the anchor if there is a region without a boundary gap at
    # all (`_circular_dart_region` would otherwise run twice for nothing).
    anchor = any(R -> isempty(R.gaps), ro) ?
             _circular_gapless_region_map(g_old, g_new) : Dict{Int,Int}()
    for R in ro
        f = labels_old[R.id]
        local k
        if isempty(R.gaps)
            # Cell without a boundary gap (free circle, lens, barbell): not
            # mappable via the boundary — but mappable via an INNER anchor
            # (`_circular_gapless_region_map`: identical edge ⇒ identical dart ⇒
            # same region). If that finds no anchor either, the default
            # applies: trivial ⇒ nothing to save, otherwise throw.
            kk = get(anchor, R.id, 0)
            if kk == 0
                if !isone(f)
                    outer === nothing && error(
                        "_circular_transfer_labels: old region $(R.id) without boundary gaps carries " *
                        "the non-trivial label $f — transfer not possible")
                    outer[] = outer[] * f
                end
                continue
            end
            k = kk
        else
            kg = findfirst(S -> R.gaps[1] in S.gaps, rn)
            kg === nothing && error(
                "_circular_transfer_labels: gap $(R.gaps[1]) of old region $(R.id) " *
                "is in no region of the new graph — unexpected")
            k = kg
        end
        haskey(overrides, k) && continue        # override wins, transfer skipped
        if haskey(claimed_by, k)
            alt = labels_new[k]
            if isone(alt)
                labels_new[k] = f               # (1, f) ⇒ f  (also f == 1)
            elseif !isone(f) && f != alt
                error("_circular_transfer_labels: old regions $(claimed_by[k]) and $(R.id) " *
                      "land in the same new region $k, both with " *
                      "non-trivial labels ($alt resp. $f) — label bookkeeping unclear")
            end
            # (f, alt) with f == 1 or f == alt: `alt` stays.
        else
            labels_new[k] = f
        end
        claimed_by[k] = R.id
    end
    for (k, lab) in overrides
        1 <= k <= length(rn) ||
            error("_circular_transfer_labels: override for invalid new region $k")
        labels_new[k] = lab
    end
    return labels_new
end

# ---- fusion ---------------------------------------------------------------------

"""
    _circular_fuse_edge(fd::CircularDecorated, ei::Int, fside::Symbol) -> Union{Nothing, CircularComboR}

Fusion rule at edge `ei` (colour `i`) separating two different **regions**
`l`/`r` (a "genuine strand"): `f · (strand) · g = (strand, label 1 left,
g·act(i,f) right) + (broken strand, two fresh degree-1 dots(i), merged region
with g·demazure(i,f))`. `fside` determines which side is the f-side. `nothing`
if both sides see the same region.

**Darts, not faces.** The two sides come from `_circular_dart_region`
(circular/CircularRegion.jl): the region left of the dart `e.a -> e.b` and left of its
reverse dart. Same left/right convention as the face-based `g.cells.port_cell`, just
without the detour through faces: when several regions lie in the same face, a
face-based lookup cannot disambiguate them, even though the dart assignment
already knows the pair.
"""
function _circular_fuse_edge(fd::CircularDecorated, ei::Int, fside::Symbol)
    g = fd.graph
    e = g.edges[ei]
    dart_region, t = _circular_dart_region(g)
    haskey(t.port_dart, e.a) || return nothing
    d = t.port_dart[e.a]         # dart e.a -> e.b; the region to its LEFT is the f side
    l = get(dart_region, d, 0)
    r = get(dart_region, t.darts[d].rev, 0)
    fside === :right && ((l, r) = (r, l))
    (l >= 1 && r >= 1 && l != r) || return nothing
    return _circular_fuse_edge_regions(fd, ei, l, r)
end

"""
    _circular_fuse_edge_regions(fd::CircularDecorated, ei::Int, l::Int, r::Int) -> CircularComboR

Same fusion step, but with the regions given **explicitly**: `l` is the
f-side, `r` the g-side. Used by [`circular_fusion_step`](@ref), which already knows
the region pair from `circular_region_adjacency` anyway.

**Why this is needed.** The `fside` variant derives the two sides from
`g.cells.port_cell`, i.e. from FACES, and cannot determine them once a face
contains several regions — a diagram can have all of its regions lie in a
single face, making the face-based lookup inherently ambiguous even though
the adjacency already knows the right pair.
"""
function _circular_fuse_edge_regions(fd::CircularDecorated, ei::Int, l::Int, r::Int)
    g = fd.graph
    e = g.edges[ei]
    i = e.colour
    f = fd.region_labels[l]
    glab = fd.region_labels[r]

    # ---- Term 1: polynomial pushed through — strand stays -----------------------
    labels1 = copy(fd.region_labels)
    labels1[l] = one(SoergelPoly)
    labels1[r] = glab * act(i, f)
    # Fusion does not touch the outer region (`l`/`r` are both inner regions,
    # callers require `l >= 1 && r >= 1`) — the outer label is carried through
    # unchanged.
    term1 = CircularDecorated(g, labels1, fd.outer_label)

    # ---- Term 2: strand broken — two fresh dots, regions merge ------------------
    gplus = CircularGraph(g.word, vcat(g.nodes, [circular_node([i]), circular_node([i])]), g.edges)
    dA, dB = length(gplus.nodes) - 1, length(gplus.nodes)
    keep = [ee for (j, ee) in enumerate(gplus.edges) if j != ei]
    push!(keep, Edge(i, e.a, NodePort(dA, 1)))
    push!(keep, Edge(i, e.b, NodePort(dB, 1)))
    g2 = _circular_delete_nodes(gplus, Set{Int}(), keep)
    # `dead` is empty, so node numbering in `g2` matches `gplus` — `dA`/`dB`
    # remain valid node indices.
    #
    # ---- ADDRESSING the merged region ------------------------------------------
    # Local address: the merged region is exactly the one the two freshly
    # placed dots `dA`/`dB` lie in — the strand that separated `l` from `r` is
    # gone, and dots are transparent for the region decomposition
    # (`circular_region_adjacency` skips dot edges). `circular_region_of_dot`
    # (circular/CircularRegion.jl) gives this region directly from planarity
    # (`_circular_dart_region`), with no reference to the boundary. Both dots MUST
    # see the same region; if they don't, `l`/`r` are still separated by a
    # second edge and the fusion term would not be well-defined at all — then
    # better to throw.
    nidA = circular_region_of_dot(g2, dA)
    nidB = circular_region_of_dot(g2, dB)
    nidA == nidB || error(
        "_circular_fuse_edge: the two new dots lie in different regions " *
        "($nidA resp. $nidB) — edge $ei does not separate $l/$r alone")
    nid = nidA
    labels2 = _circular_transfer_labels(g, g2, fd.region_labels;
                                   overrides = Dict(nid => glab * demazure(i, f)))
    term2 = CircularDecorated(g2, labels2, fd.outer_label)

    return CircularCombo{SoergelPoly}(term1) + CircularCombo{SoergelPoly}(term2)
end

"""
    _fr_fusion(fd::CircularDecorated, ei::Int) -> Union{Nothing, CircularComboR}

Fusion rule at edge `ei` (f-side = left, as in `fuse_edge`).
"""
_fr_fusion(fd::CircularDecorated, ei::Int) = _circular_fuse_edge(fd, ei, :left)
