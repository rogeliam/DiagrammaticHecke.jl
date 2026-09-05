# CORE — the shared diagram/word layer under the circular pipeline.
# diagram/CanonicalKey.jl — the ONE canonical-key implementation, shared by
# `canonical_key(::WordGraph)` (diagram/Combo.jl) and `circular_canonical_key(::CircularGraph)`
# (circular/CircularGraph.jl).
#
# ══ THE METHOD: A DART WALK, NOT WEISFEILER–LEMAN ══
#
# Fixing the node numbering by WL colour refinement instead would need to resolve
# two ambiguities by EXHAUSTIVE SEARCH: the arm origin of each circular node (several
# when the arm sequence is rotation-symmetric), as an `Iterators.product` OVER ALL
# NODES, and every permutation within each WL class. Both multiply, giving a
# combinatorial explosion in the node count.
#
# A `CircularGraph` is a ribbon graph: every node carries a CYCLIC arm sequence. For such
# graphs the canonical form is a DART WALK, not a minimum over origin combinations.
# A dart is an edge end, i.e. a pair (node, arm). Walking the graph from a start
# dart fixes BOTH quantities that an exhaustive search would otherwise guess or
# minimise:
#
#   * the node number = the order of discovery,
#   * the arm origin  = the arm through which the node was first reached.
#
# Nothing is chosen, everything follows. Over all start darts we take the
# lexicographically smallest certificate: 2·|edges| walks instead of a product over
# the nodes — LINEAR in the edges instead of multiplicative in the node symmetries.
#
# COMPLETENESS. The walk is a genuine canonical form for connected graphs: two
# graphs are equal iff some start dart yields the same certificate — stronger than
# colour refinement plus class permutations, which could in principle merge
# non-isomorphic graphs.
#
# WHY THE TWO CALLERS ARE PARAMETRISED DIFFERENTLY.
#
#   1. NODE FINGERPRINT (`node_fp`). A `Node` (diagram/Graph.jl) has a fixed slot
#      layout per `kind`: slot 1 of a trivalent node is a named place, so
#      `_node_fp(nd) = (kind, colours, m)` describes the node completely and needs
#      no rotation normalisation. A `CircularNode` has an open arm count and no
#      layout to name a first arm by, so two arm sequences that are rotations of
#      each other describe the same node; `circular_key(nd) = (kind,
#      min_rotation(arms))` divides that symmetry out before the walk sees it.
#
#   2. FREE OR FIXED ORIGIN (`free_origin`). The same difference one level up, for
#      the slots the walk records. On a `CircularGraph` an arm has no absolute
#      number, only a position relative to some other arm, so the origin is free:
#      the walk fixes it at the arm through which the node was first reached and
#      counts slots from there. On a `WordGraph` the slot number IS part of the
#      node's meaning, so the walk records it raw and rotating a node gives a
#      different key — which is what the fixed layout of `Node` asks for.
#
#   3. BOUNDARY PORTS. Both take the raw boundary position `p.k`. The boundary is
#      absolute (the word is fixed in normal form) and hence untouched by the walk —
#      it is also what usually pins down connected diagrams immediately.
#
# FORMAT. Rely only on the equality the key induces, never on the numeric values.

# ---- ports: the part INDEPENDENT of the walk --------------------------------
#
# Leaf and circle already carry their own index-independent identity.

_dart_port_desc(p::Leaf)   = (0, p.k, 0)
_dart_port_desc(p::Circle) = (2, p.colour, 0)

"""
    _dart_arm_table(nodes, edges, degree) -> (arm, lose)

`arm[i][s]` = all `(edge colour, opposite port)` at slot `s` of node `i` (a list,
since a loop may have both ends at the same slot). `lose` = sorted fingerprints of
the edges attached to NO node (leaf–leaf, circles) — they belong to no component and
enter the key separately.
"""
function _dart_arm_table(nodes::AbstractVector, edges::Vector{Edge}, degree)
    n = length(nodes)
    arm = [[Tuple{Int,Any}[] for _ in 1:degree(nodes[i])] for i in 1:n]
    lose = Tuple[]
    for e in edges
        at = e.a isa NodePort
        bt = e.b isa NodePort
        at && push!(arm[e.a.node][e.a.slot], (e.colour, e.b))
        bt && push!(arm[e.b.node][e.b.slot], (e.colour, e.a))
        if !at && !bt
            x, y = _dart_port_desc(e.a), _dart_port_desc(e.b)
            push!(lose, (e.colour, min(x, y), max(x, y)))
        end
    end
    return arm, sort(lose)
end

# FLAT INT CERTIFICATE. Every record starts with a NEGATIVE sentinel (−1 node, −2
# arm) and all payload is ≥ 0, which makes the encoding self-delimiting and
# injective, so lexicographic `isless(::Vector{Int}, ::Vector{Int})` replaces dynamic
# dispatch and sorting of stringified tuples. The node kind enters as a small
# number; the arm count sits immediately after the kind id in the certificate
# (`head!` in CircularGraph.jl), so 6- and 8-armed braids stay distinguishable.
const _CK_KIND_ID = Dict{Symbol,Int}(
    :dot => 1, :trivalent => 2, :braid => 3, :mono => 4, :mixed => 5)
_ck_kind_id(k::Symbol) =
    get(_CK_KIND_ID, k) do
        error("CanonicalKey: unknown node kind :$k — add it to _CK_KIND_ID")
    end

"""
    _dart_certificate(nodes, arm, degree, head!, free_origin, start_node, start_slot)
        -> (certificate::Vector{Int}, component)

The walk starting from the dart `(start_node, start_slot)`. Nodes are numbered in
discovery order; the arm through which a node is first reached is its origin. Under
`free_origin` the arms are traversed cyclically from that origin and slots counted
relative to it; otherwise from slot 1 with raw slot numbers. `component` is the list
of visited nodes.

`head!(cert, nd, o)` writes the head of a node into the certificate as a sequence of
Ints (values ≥ 0 only!) — see the sentinel block above.
"""
function _dart_certificate(nodes::AbstractVector, arm, degree, head!, free_origin::Bool,
                           start_node::Int, start_slot::Int)
    n = length(nodes)
    num    = zeros(Int, n)          # 0 = not yet discovered
    origin = zeros(Int, n)
    order  = Int[]

    function discover!(i::Int, o::Int)
        num[i] = length(order) + 1
        origin[i] = free_origin ? o : 1
        push!(order, i)
    end

    discover!(start_node, start_slot)
    cert = Int[]
    pos = 1
    while pos <= length(order)
        i = order[pos]; pos += 1
        d = degree(nodes[i])
        o = origin[i]
        push!(cert, -1)
        head!(cert, nodes[i], o)
        for t in 0:(d - 1)
            s = free_origin ? mod1(o + t, d) : t + 1
            for (ecol, q) in arm[i][s]
                if q isa NodePort
                    num[q.node] == 0 && discover!(q.node, q.slot)
                    dq = degree(nodes[q.node])
                    cslot = free_origin ? mod1(q.slot - origin[q.node] + 1, dq) : q.slot
                    push!(cert, -2, t, ecol, 1, num[q.node], cslot)
                else
                    x = _dart_port_desc(q)
                    push!(cert, -2, t, ecol, x[1], x[2], x[3])
                end
            end
        end
    end
    return cert, order
end

"""
    _dart_canonical_key(nodes, edges, boundary_letters, node_fp, degree, head,
                        free_origin) -> Tuple

The full canonical key: boundary letters (already in normal form from the caller) +
sorted node fingerprints + the lexicographically smallest dart certificate per
component (sorted) + the node-less edges.

`head!(cert, nd, origin)` writes the head of a node into the certificate as a
sequence of Ints (values ≥ 0, see the sentinel block at `_dart_certificate`),
relative to the chosen origin: `CircularGraph` writes kind id + arm count + the arm
sequence rotated FROM THE ORIGIN (the origin is structurally determined, hence NOT
the minimal rotation), `WordGraph` writes kind id + colour count + colours + `m`.
"""
function _dart_canonical_key(nodes::AbstractVector, edges::Vector{Edge}, boundary_letters,
                             node_fp, degree, head!, free_origin::Bool;
                             anchored::Bool = true)
    arm, lose = _dart_arm_table(nodes, edges, degree)
    n = length(nodes)
    nds = sort([node_fp(nd) for nd in nodes])

    # BOUNDARY ANCHOR. The boundary is absolute (the word is fixed in normal form,
    # `Leaf(k)` carries its `k` raw in the certificate), so for every component
    # touching the boundary the canonical start dart is FORCED: the smallest leaf of
    # the component, entering along its edge. ONE walk instead of 2·|E|. Correctness:
    # equal certificates ⟺ there is an iso preserving the (absolute) leaf numbers ⟹
    # it maps the smallest leaf to the smallest leaf, so the anchor darts correspond
    # and the anchored certificates agree — the INDUCED equality is the same as under
    # full minimisation. Only boundary-less components still minimise over all
    # darts. `anchored = false` keeps the full minimisation as a reference.
    leafanch = NTuple{3,Int}[]           # (k, node, slot), ascending in k
    if anchored
        for e in edges
            if e.a isa Leaf && e.b isa NodePort
                push!(leafanch, (e.a.k, e.b.node, e.b.slot))
            elseif e.b isa Leaf && e.a isa NodePort
                push!(leafanch, (e.b.k, e.a.node, e.a.slot))
            end
        end
        sort!(leafanch)
    end

    seen = falses(n)
    certs = Vector{Int}[]
    for i in 1:n
        seen[i] && continue
        # walk once, any way, just to learn the component. A node without arms has
        # no dart — then the single walk from slot 1 is already the answer.
        best, comp = _dart_certificate(nodes, arm, degree, head!, free_origin, i, 1)
        anchor = nothing
        if anchored && !isempty(leafanch)
            incomp = falses(n)
            for j in comp
                incomp[j] = true
            end
            ix = findfirst(a -> incomp[a[2]], leafanch)
            ix === nothing || (anchor = leafanch[ix])
        end
        if anchor !== nothing
            # boundary component: the anchor fixes the start, nothing is minimised.
            best, _ = _dart_certificate(nodes, arm, degree, head!, free_origin,
                                        anchor[2], anchor[3])
        else
            # boundary-less component: try every dart as a start and keep the
            # smallest certificate (lexicographic `isless` on `Vector{Int}`).
            for j in comp, s in 1:degree(nodes[j])
                c, _ = _dart_certificate(nodes, arm, degree, head!, free_origin, j, s)
                isless(c, best) && (best = c)
            end
        end
        for j in comp
            seen[j] = true
        end
        push!(certs, best)
    end
    return (Tuple(boundary_letters), Tuple(nds), sort!(certs), Tuple(lose))
end
