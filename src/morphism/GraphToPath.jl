# morphism/GraphToPath.jl — peel a WordGraph outside-in into a reduction path.

# ---- graph → path : peel vertices outside-in -------------------------------

# an "open port" carries its colour and its Port (leaf or node connector).
struct _OpenPort
    colour::Int
    port::Port
end

"""
    graph_to_path(g::WordGraph) -> (PathDiagram, status::Symbol)

Peel `g` into a reduction path (layer/vertex sequence) from the boundary inward.
`status ∈ (:ok, :stalled)`. Paths are not unique; this returns ONE. Stalls (rather
than guessing) when no leaf-cap and no contiguous-window vertex is available — the
planarity approximation.
"""
graph_to_path(g::WordGraph) = (r = _graph_to_path_full(g); (r[1], r[2]))

# Full peel including per-step render data: (PathDiagram, status, rings, steps), where
# rings/steps are in the `_radial_core` format (see `_shell_layers`).
function _graph_to_path_full(g::WordGraph)
    n = length(g.word)
    # the open ring: one entry per current boundary strand, cyclically ordered.
    open = _OpenPort[_OpenPort(leaf_colour(g, k), Leaf(k)) for k in 1:n]
    layers = Vector{Int}[[p.colour for p in open]]
    verts  = Vertex[]
    # per-step render data for the radial/shell renderer (see `_radial_core`):
    # steps[i] = (wins, surv, repl_idx). Built alongside the peel; ignored by callers
    # that only want (PathDiagram, status).
    steps = Tuple[]
    remaining = Set(1:length(g.nodes))
    # a fast port→(edge, other port) lookup that ignores already-consumed edges.
    consumed = Set{Int}()
    other_end(p::Port) = begin
        for (ei, e) in enumerate(g.edges)
            ei in consumed && continue
            if _same_port(e.a, p); return (ei, e.b, e.colour); end
            if _same_port(e.b, p); return (ei, e.a, e.colour); end
        end
        return nothing
    end

    guard = 0
    while !isempty(open) && guard < 10_000
        guard += 1
        progressed = false

        # (1) a cap: two ADJACENT open ports joined directly by one edge (ii → ε).
        for i in 1:length(open)
            j = mod1(i + 1, length(open))
            i == j && continue
            pi = open[i].port
            oe = other_end(pi)
            oe === nothing && continue
            (ei, far, col) = oe
            _same_port(far, open[j].port) || continue
            # emit ii → ε at window position i
            push!(consumed, ei)
            push!(verts, Vertex(VDot, [col], i, 2, Int[], Int[]))
            # render tracking: the two capped positions vanish, nothing is produced;
            # every OTHER old position k maps to its index after deleting {i,j}.
            gone = Set(sort([i, j]))
            surv = Dict{Int,Int}()
            newpos = 0
            for k in 1:length(open)
                k in gone && continue
                newpos += 1
                surv[k] = newpos
            end
            push!(steps, ([(i, 2, -2, Int[])], surv, [Int[]]))
            deleteat!(open, sort([i, j]))
            push!(layers, [p.colour for p in open])
            progressed = true
            break
        end
        progressed && continue

        # (2) a node covering a CONTIGUOUS cyclic window of open ports. Prefer
        # dot < trivalent < braid, smallest window start.
        best = nothing   # (score, nodeindex, window_positions, far_ports_in_order)
        for v in remaining
            nd = g.nodes[v]
            deg = _node_degree(nd)
            # find the open-ports that connect to node v, and their window positions
            hits = Int[]     # positions in `open` that connect to v
            portslot = Dict{Int,Int}()   # open position → v's slot
            for (pos, op) in enumerate(open)
                oe = other_end(op.port)
                oe === nothing && continue
                (_, far, _) = oe
                if far isa NodePort && far.node == v
                    push!(hits, pos); portslot[pos] = far.slot
                end
            end
            isempty(hits) && continue
            # We do NOT require ALL of v's open legs to be adjacent — only that SOME
            # contiguous block of adjacent ring positions, all belonging to v, is a valid
            # peel. E.g. legs at ring positions [1,3,4,5] still let a
            # braid peel using the adjacent sub-block {3,4,5}. `_max_contiguous_run` finds
            # the longest such run of v's legs; `wpos` is that run.
            wpos = _max_contiguous_run(sort(hits), length(open))
            wpos === nothing && continue
            # the window must be a VALID peel for this node kind:
            #   dot        window 1                    (i → ε)
            #   trivalent  window 1 or 2               (i → ii  /  ii → i)
            #   braid      ANY contiguous window 1..2m of its legs:
            #     m legs = a genuine braid move; 2m legs = the whole vertex is exposed
            #     (e.g. 121212 → ε); a partial window unfolds the rest inward. Rings may
            #     GROW — fine, finitely many vertices still drive toward ε.
            valid_win =
                nd.kind === :dot        ? length(wpos) == 1 :
                nd.kind === :trivalent  ? length(wpos) in (1, 2) :
                nd.kind === :braid      ? 1 <= length(wpos) <= 2 * nd.m :
                false
            valid_win || continue
            # prefer the peel that grows the ring LEAST: a window of
            # size w on a degree-`deg` vertex changes the open-ring size by deg − 2w, so
            # a full vertex (w = deg) shrinks it most, a half-braid (w = deg/2) is neutral,
            # a small window grows it. Order by that delta (most-shrinking first), then
            # dot < trivalent < braid, then smallest start. Short rings ⇒ nicer concentric
            # picture; brief growth is fine.
            deg = _node_degree(nd)
            delta = deg - 2 * length(wpos)          # ring-size change (want most negative)
            score = nd.kind === :dot ? 1 : nd.kind === :trivalent ? 2 : 3
            key = (delta, score, minimum(wpos))
            if best === nothing || key < best[1]
                best = (key, v, wpos, portslot)
            end
        end
        if best !== nothing
            (_, v, wpos, portslot) = best
            nd = g.nodes[v]
            deg = _node_degree(nd)
            s = wpos[1]; pl = length(wpos)
            # consume the window edges; build replacement = v's OTHER slots (the ones
            # not in the window), wired as new open ports where the window was.
            window_slots = Set(portslot[p] for p in wpos)
            for p in wpos
                oe = other_end(open[p].port)
                oe === nothing && continue
                push!(consumed, oe[1])
            end
            rest_slots = [sl for sl in 1:deg if !(sl in window_slots)]
            repl = Int[_slot_colour(nd, sl) for sl in rest_slots]
            newopen = [_OpenPort(_slot_colour(nd, sl), NodePort(v, sl)) for sl in rest_slots]
            # window letters (for classification)
            window = [open[p].colour for p in wpos]
            kind, cols = _classify(window, repl)
            push!(verts, Vertex(kind, cols, s, pl, repl, Int[]))
            # render tracking: rebuild the next ring as [replacement block] ++
            # [survivors in cyclic order after the window]. This fixed order lets us map
            # old→new positions exactly. `wpos` is a contiguous (maybe wrapping) window;
            # `wend` = its last position in cyclic order (the one whose +1 is a survivor).
            rl = length(newopen)
            oldlen = length(open)
            winset = Set(wpos)
            wend = wpos[end]
            for p in wpos
                (mod1(p + 1, oldlen) in winset) || (wend = p)
            end
            # replacement ports take new positions 1..rl; survivors follow at rl+1, …
            surv = Dict{Int,Int}()
            survivors_ordered = _OpenPort[]
            pos = mod1(wend + 1, oldlen); np = rl; cnt = 0
            while cnt < oldlen - length(wpos)
                if !(pos in winset)
                    np += 1; surv[pos] = np
                    push!(survivors_ordered, open[pos])
                end
                pos = mod1(pos + 1, oldlen); cnt += 1
            end
            replidx = collect(1:rl)
            open = vcat(newopen, survivors_ordered)
            push!(steps, ([(1, pl, rl - pl, repl)], surv, [replidx]))
            delete!(remaining, v)
            push!(layers, [p.colour for p in open])
            progressed = true
        end

        progressed || return (PathDiagram(layers, verts), :stalled, layers, steps)
    end
    return (PathDiagram(layers, verts), :ok, layers, steps)
end

# positions (sorted) form a contiguous cyclic window in a ring of length L?
# return the window as an ordered position list (start..) or nothing.
# The longest run of CONSECUTIVE ring positions (mod L) all contained in `pos`.
# Unlike `_contiguous_window` (which needs the WHOLE set contiguous), this finds the
# best adjacent sub-block — so a node whose legs are scattered can still peel the
# adjacent ones. Returns the run as an ordered position vector, or nothing if empty.
function _max_contiguous_run(pos::Vector{Int}, L::Int)
    isempty(pos) && return nothing
    S = Set(pos)
    length(S) == L && return collect(1:L)          # everything open — whole ring
    best = Int[]
    for start in pos
        # a run can only START where the previous position is NOT in S
        (mod1(start - 1, L) in S) && continue
        run = Int[]
        p = start
        while p in S && length(run) < L
            push!(run, p); p = mod1(p + 1, L)
        end
        length(run) > length(best) && (best = run)
    end
    return isempty(best) ? nothing : best
end

function _contiguous_window(sorted_pos::Vector{Int}, L::Int)
    k = length(sorted_pos)
    k == 0 && return nothing
    k == L && return sorted_pos
    # try each start; a window is contiguous if positions are consecutive mod L
    for start in sorted_pos
        w = [mod1(start + j, L) for j in 0:(k-1)]
        Set(w) == Set(sorted_pos) && return w
    end
    return nothing
end

# replace the contiguous positions `wpos` (ordered) in `open` by `newopen`.
function _splice_window!(open::Vector{_OpenPort}, wpos::Vector{Int}, newopen::Vector{_OpenPort})
    # wpos may wrap; delete then insert at the (adjusted) start.
    start = wpos[1]
    sorted = sort(wpos)
    # if the window is contiguous without wrap, simple splice
    if sorted == collect(start:start+length(wpos)-1)
        splice!(open, start:start+length(wpos)-1, newopen)
    else
        # wrapping window: delete highest-first, then insert newopen at position 1
        for p in sort(wpos; rev = true); deleteat!(open, p); end
        for (i, np) in enumerate(newopen); insert!(open, i, np); end
    end
    return open
end

