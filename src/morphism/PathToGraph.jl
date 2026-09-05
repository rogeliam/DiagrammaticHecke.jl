# morphism/PathToGraph.jl — wire a reduction path back into a WordGraph.
# Also contains `path_from_words`.

# ---- path → graph : wire the layer/vertex sequence -------------------------

"""
    path_to_graph(d::PathDiagram) -> WordGraph

Materialise a reduction path into a WordGraph on boundary `d.layers[1]`. `ii → ε`
becomes a plain cap edge (no node); `i → ε` a dot; `ii → i` a trivalent; a braid /
commutation a braid node. Result satisfies `check_valid` when the path ends at ε.
"""
function path_to_graph(d::PathDiagram)
    n = length(d.layers[1])
    word = CircularWord(copy(d.layers[1]))
    # open ports waiting to be wired, in ring order; start = the boundary leaves.
    open = Port[Leaf(k) for k in 1:n]
    opencol = copy(d.layers[1])
    nodes = Node[]
    edges = Edge[]

    for v in d.vertices
        L = length(open)
        wpos = [mod1(v.pos + j, L) for j in 0:(v.pl-1)]
        winports = [open[p] for p in wpos]
        wincols  = [opencol[p] for p in wpos]
        rl = length(v.repl)
        if v.pl == 1 && rl == 0                                   # i → ε : dot
            push!(nodes, Node(:dot, [wincols[1]], 0))
            nn = length(nodes)
            push!(edges, Edge(wincols[1], winports[1], NodePort(nn, 1)))
            _splice_ports!(open, opencol, wpos, Port[], Int[])
        elseif v.pl == 2 && rl == 0                               # ii → ε : cap edge
            push!(edges, Edge(wincols[1], winports[1], winports[2]))
            _splice_ports!(open, opencol, wpos, Port[], Int[])
        elseif v.pl == 2 && rl == 1                               # ii → i : trivalent
            c = v.repl[1]
            push!(nodes, Node(:trivalent, [c], 0)); nn = length(nodes)
            push!(edges, Edge(c, winports[1], NodePort(nn, 1)))
            push!(edges, Edge(c, winports[2], NodePort(nn, 2)))
            _splice_ports!(open, opencol, wpos, Port[NodePort(nn, 3)], Int[c])
        elseif v.pl == 1 && rl == 2                               # i → ii : split (trivalent)
            c = wincols[1]
            push!(nodes, Node(:trivalent, [c], 0)); nn = length(nodes)
            push!(edges, Edge(c, winports[1], NodePort(nn, 3)))
            _splice_ports!(open, opencol, wpos,
                                  Port[NodePort(nn, 1), NodePort(nn, 2)], Int[c, c])
        elseif v.pl == rl && v.pl in (2, 3) &&
               length(unique(vcat(wincols, v.repl))) == 2         # genuine braid / comm move
            cols = sort(unique(vcat(wincols, v.repl)))
            s, t = cols[1], cols[2]
            m = v.pl                                              # 13↔31 (m=2) vs 121↔212 (m=3)
            push!(nodes, Node(:braid, [s, t], m)); nn = length(nodes)
            # wire the window ports to the first pl slots (colour-matching)
            for (j, p) in enumerate(wpos)
                push!(edges, Edge(wincols[j], winports[j], NodePort(nn, j)))
            end
            # replacement = the remaining slots pl+1 … 2m
            restslots = collect((v.pl + 1):(2m))
            newp = Port[NodePort(nn, sl) for sl in restslots]
            newc = Int[_slot_colour(nodes[nn], sl) for sl in restslots]
            _splice_ports!(open, opencol, wpos, newp, newc)
        elseif v.pl != rl && (v.pl + rl) in (4, 6) && iseven(v.pl + rl) &&
               length(unique(vcat(wincols, v.repl))) == 2 &&
               all(vcat(wincols, reverse(v.repl))[j] !=
                   vcat(wincols, reverse(v.repl))[mod1(j + 1, v.pl + rl)]
                   for j in 1:(v.pl + rl))
            # LENGTH-CHANGING braid move: e.g. 1 ↔ 313 (m=2),
            # 13 ↔ 31 is length-neutral but 1↔313 / 3↔131 are the "seen from the side"
            # braid. One braid(s,t;m) vertex, 2m = pl+rl legs: the pl window strands on
            # slots 1..pl, the rl replacement strands on slots pl+1..2m (going around the
            # node so colours alternate). No leg is left over — all 2m are used.
            m = (v.pl + rl) ÷ 2
            allcols = vcat(wincols, v.repl)
            cols = sort(unique(allcols)); s, t = cols[1], cols[2]
            push!(nodes, Node(:braid, [s, t], m)); nn = length(nodes)
            # window strands → slots 1..pl (colour-matching)
            for (j, p) in enumerate(wpos)
                push!(edges, Edge(wincols[j], winports[j], NodePort(nn, j)))
            end
            # replacement strands → slots pl+1..2m
            restslots = collect((v.pl + 1):(2m))
            newp = Port[NodePort(nn, sl) for sl in restslots]
            newc = Int[_slot_colour(nodes[nn], sl) for sl in restslots]
            _splice_ports!(open, opencol, wpos, newp, newc)
        elseif rl == 0 && v.pl == 3 && length(unique(wincols)) == 1
            # trivalent CLOSURE: three same-colour strands iii → ε meet at ONE trivalent
            # vertex whose 3 legs are exactly these strands.
            c = wincols[1]
            push!(nodes, Node(:trivalent, [c], 0)); nn = length(nodes)
            for (j, p) in enumerate(wpos)
                push!(edges, Edge(c, winports[j], NodePort(nn, j)))
            end
            _splice_ports!(open, opencol, wpos, Port[], Int[])
        elseif rl == 0 && iseven(v.pl) && v.pl >= 4 &&
               length(unique(wincols)) == 2 &&
               all(wincols[j] == wincols[mod1(j + 2, v.pl)] for j in 1:v.pl)
            # braid CLOSURE: an alternating window s t s t … (length 2m) closes into ONE
            # braid(m) vertex whose 2m legs are exactly these strands (window → ε). This
            # is the "big braid" that ends R6 (212121 → ε), not a length-neutral move.
            m = v.pl ÷ 2
            s, t = wincols[1], wincols[2]
            push!(nodes, Node(:braid, [s, t], m)); nn = length(nodes)
            for (j, p) in enumerate(wpos)
                push!(edges, Edge(wincols[j], winports[j], NodePort(nn, j)))
            end
            _splice_ports!(open, opencol, wpos, Port[], Int[])
        else
            # a COMPOSITE reducer move (e.g. 131→3, 1212→ε from EXTRA_REDUCERS) is not a
            # single diagram generator; we cannot express it as one node. This is the
            # documented limitation: path_to_graph handles generator moves
            # (dot / cap / trivalent / braid / comm) only.
            error("path_to_graph: move (pl=$(v.pl), rl=$rl, repl=$(v.repl)) is not a " *
                  "single generator — composite reducer move, not supported.")
        end
    end
    return WordGraph(word, nodes, edges)
end

"""
    path_from_words(words) -> PathDiagram

Build a `PathDiagram` from an EXPLICIT sequence of words (each a `Vector{Int}`), one
per line of a reduction. Between consecutive words the single changed CONTIGUOUS
window is found and classified (dot / trivalent / braid / commutation) — so the path
can describe a diagram simply by writing its reduction path, and `path_to_graph` then
gives the unique graph. ("if I describe a path there is a unique
graph".) The words are non-cyclic here (the window must be contiguous, not wrapping);
that covers the light-leaf / stacked-morphism reading. Consecutive words must differ
by exactly one elementary move.
"""
function path_from_words(words::Vector{Vector{Int}})
    isempty(words) && error("path_from_words: need at least one word")
    layers = [copy(w) for w in words]
    verts = Vertex[]
    for i in 1:(length(words) - 1)
        a = words[i]; b = words[i + 1]
        la, lb = length(a), length(b)
        p = 0
        while p < min(la, lb) && a[p + 1] == b[p + 1]; p += 1; end
        s = 0
        while s < min(la, lb) - p && a[la - s] == b[lb - s]; s += 1; end
        window = a[(p + 1):(la - s)]
        repl   = b[(p + 1):(lb - s)]
        kind, cols = _classify(window, repl)
        ridx = collect((p + 1):(lb - s))
        push!(verts, Vertex(kind, cols, p + 1, length(window), copy(repl), ridx))
    end
    return PathDiagram(layers, verts)
end

# splice ports+colours together (parallel arrays).
function _splice_ports!(open::Vector{Port}, opencol::Vector{Int},
                               wpos::Vector{Int}, newp::Vector{Port}, newc::Vector{Int})
    L = length(open)
    start = wpos[1]
    if sort(wpos) == collect(start:start+length(wpos)-1)
        splice!(open, start:start+length(wpos)-1, newp)
        splice!(opencol, start:start+length(wpos)-1, newc)
    else
        for p in sort(wpos; rev = true); deleteat!(open, p); deleteat!(opencol, p); end
        for (i, (np, nc)) in enumerate(zip(newp, newc))
            insert!(open, i, np); insert!(opencol, i, nc)
        end
    end
    return open
end

# ---------------------------------------------------------------------------
# Rich display for MorphismGraph: draw the underlying WordGraph, same as for a
# plain graph (WordGraph's show methods in render/Graph.jl). The cuts are
# structural bookkeeping; the picture is the diagram itself.
# ---------------------------------------------------------------------------

# the rich-display `show(::MIME, ::MorphismGraph)` methods live in
# render/tutte/Display.jl (Tutte is the default layout).
Base.showable(::MIME"image/svg+xml", ::MorphismGraph) = true
Base.showable(::MIME"image/png", ::MorphismGraph) = true
