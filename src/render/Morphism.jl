# render/Morphism.jl  —  automatic radial rendering of a word reduction (a morphism)
#
# Renders a reduction of a circular word (a chosen down-path to ε) as a RADIAL flow
# SVG: centre = ε, each reduction step = one concentric ring further out, outermost
# ring = the start word on its boundary circle. The word stays a circle on every
# ring, so there is NO seam and no cut.
#
# Input is generated internally by `reduction_with_moves`, so every ring carries
# reliable move info (which cyclic window, which rule) — no log parsing, no
# normal-form/rotation mismatch. Angular re-spreading, neutral moves and crowding
# are not handled.

# ---- reduction path with move info -----------------------------------------

# One reduction step: the word, the cyclic window (pos, patlen) rewritten, the rule
# that fired, and the resulting word.
struct RedStep
    word::CircularWord
    pos::Int
    patlen::Int
    rule::Symbol
    result::CircularWord
end

const _EXT_BY_NAME = Dict(r.name => r for r in EXTENDED_RULES)

"""
    reduction_with_moves(w; maxsteps = 200) -> Vector{RedStep}

A greedy strictly-reducing path from `w` down to ε (biggest length drop first),
each step carrying the exact cyclic window and rule that fired. Self-consistent:
`step.result` is the next step's `word`, ending at ε (empty result) if `w` reduces.
Returns the steps bottom→top (start first). Empty if `w` has no reducing move.
"""
function reduction_with_moves(w::CircularWord; maxsteps::Int = 200)
    steps = RedStep[]
    cur = w
    for _ in 1:maxsteps
        isempty(cur) && break
        best = nothing
        for m in moves(cur, EXTENDED_RULES)
            m.sign < 0 || continue                     # strictly reducing only
            if best === nothing || length(m.result) < length(best.result)
                best = m
            end
        end
        best === nothing && break
        r = _EXT_BY_NAME[best.rule]
        pl = best.forward ? length(r.from) : length(r.to)
        push!(steps, RedStep(cur, best.pos, pl, best.rule, best.result))
        cur = best.result
    end
    return steps
end

# ---- colours (real RGB for SVG) --------------------------------------------

const _RGB = Dict(1 => "#007bff", 2 => "#ff4136", 3 => "#2e8b2e")
const _RGBMIX = Dict((1, 2) => "#8a2be2", (2, 3) => "#ff851b", (1, 3) => "#7ea628")
_rgbmix(s, t) = _RGBMIX[minmax(s, t)]

# colour of the node a rule produces, from the window's colours
function _node_rgb(cols)
    length(cols) == 2 && return _rgbmix(cols[1], cols[2])
    return _RGB[cols[1]]
end

# ---- weighted colour mixer --------------------------------------------------
#
# A 1/3-coloured node is mixed by ARM SHARE in RGB space, not by a single fixed mix
# colour: with `n` arms of colour 1 and `m` of colour 3 the fill is
# `(n·colour1 + m·colour3)/(n+m)`, so a `[1,1,1,3]` node does not look like a
# `[1,3,3,3]` one.
# The hue then shows which colour dominates the node.
#
# The two-colour special case n == m does NOT return exactly `_RGBMIX` (those are
# hand-picked, more saturated tones) — the mixer is its own continuous scale, not
# an interpolation between those constants.

_hex_rgb(h::AbstractString) = (parse(Int, h[2:3], base = 16),
                               parse(Int, h[4:5], base = 16),
                               parse(Int, h[6:7], base = 16))

"""
    _rgb_weighted(cols) -> String

Mix colour for the MULTISET `cols` of arm colours. A single colour returns
`_RGB[c]`.

TWO colours: piecewise linear between the two pure colours WITH THE
HAND-PICKED MIX COLOUR `_RGBMIX` IN THE MIDDLE. The purely computed average of
`#007bff` and `#2e8b2e` is `#178396`, which reads as blue (blue contributes
B = 255); the midpoint is instead `#7ea628`, the same tone the thin renderer
has always used via `_node_rgb` — so circular and thin pictures colour the same
node the same way. The arm share stays legible: `[1,1,1,3]` sits between blue
and the midpoint, `[1,3,3,3]` between the midpoint and green.

Three or more colours: the weighted average in RGB space, as before.
"""
function _rgb_weighted(cols)
    isempty(cols) && return "#666666"
    us = sort(unique(cols))
    length(us) == 1 && return get(_RGB, us[1], "#666666")
    if length(us) == 2 && haskey(_RGBMIX, (us[1], us[2]))
        s, t = us[1], us[2]
        p = count(==(s), cols) / length(cols)          # share of the SMALLER colour
        lo  = _hex_rgb(_RGB[t])                        # p = 0   -> pure t
        mid = _hex_rgb(_RGBMIX[(s, t)])                # p = 0.5 -> the mix colour
        hi  = _hex_rgb(_RGB[s])                        # p = 1   -> pure s
        comp(i) = p <= 0.5 ? lo[i] + (mid[i] - lo[i]) * (p / 0.5) :
                             mid[i] + (hi[i] - mid[i]) * ((p - 0.5) / 0.5)
        return _rgb_hex(comp(1), comp(2), comp(3))
    end
    r = g = b = 0.0
    for c in cols
        (cr, cg, cb) = _hex_rgb(get(_RGB, c, "#666666"))
        r += cr; g += cg; b += cb
    end
    n = length(cols)
    return _rgb_hex(r / n, g / n, b / n)
end

_rgb_hex(r, g, b) = string("#", _hex2(r), _hex2(g), _hex2(b))
_hex2(x) = string(round(Int, clamp(x, 0, 255)), base = 16, pad = 2)

# ---- guided reduction on RAW letter vectors --------------------------------
#
# For drawing we work on the raw (un-normalized) letter vector, so positions stay
# stable from ring to ring — CircularWord's minimal-rotation normal form would
# silently rotate the word between steps and the strand positions would not line up.
#
# A purely greedy strictly-reducing pass gets STUCK (it can pick reductions that
# strand the rest in a non-reducible remnant like "32"). So we take ONE move per
# layer, allowing REDUCING and NEUTRAL (braid/commutation) moves, and we only step
# into a remnant that STILL `reduces_to_empty_down` — the oracle steers us onto a
# path that actually reaches ε. Neutral moves are needed to bring letters into a
# reducible position; a `seen` set stops neutral cycles.

# One raw move on `v`: window at cyclic start `s`, length `pl`, replaced by `repl`.
# Position-faithful: emits the result letters, `surv[k]` = result index of untouched
# source position k, `wpos` = the covered source positions, `ridx` = result indices
# of the replacement letters. dsign = |repl|-|window| (≤0: reducing or neutral).
struct _RawMove
    result::Vector{Letter}
    surv::Dict{Int,Int}
    wpos::Vector{Int}
    repl::Vector{Letter}
    ridx::Vector{Int}
    dsign::Int
end

# All reducing-or-neutral raw moves applicable to `v` (every rule, both directions
# with |to| ≤ |from|, every cyclic start).
function _raw_moves(v::Vector{Letter})
    n = length(v); out = _RawMove[]
    n == 0 && return out
    for r in EXTENDED_RULES
        for (from, to) in ((r.from, r.to), (r.to, r.from))
            length(to) <= length(from) || continue     # reducing or neutral only
            pl = length(from)
            (pl == 0 || pl > n) && continue
            for s in 1:n
                _cyclic_window(v, s, pl) == from || continue
                wpos = [mod1(s + j, n) for j in 0:(pl-1)]
                covered = Set(wpos)
                result = Letter[]; surv = Dict{Int,Int}(); ridx = Int[]
                for k in 1:n
                    if k in covered
                        if k == s                       # window start: emit replacement
                            for c in to; push!(result, c); push!(ridx, length(result)); end
                        end
                    else
                        push!(result, v[k]); surv[k] = length(result)
                    end
                end
                push!(out, _RawMove(result, surv, wpos, collect(to), ridx,
                                    length(to) - length(from)))
            end
        end
    end
    return out
end

# Pick a next move: prefer strictly reducing (most negative dsign), require the
# remnant to still reduce to ε (unless it IS ε), and skip already-seen remnants.
function _guided_next(v::Vector{Letter}, seen::Set{CircularWord})
    mv = _raw_moves(v)
    sort!(mv, by = m -> m.dsign)                        # reducing first
    for m in mv
        nw = CircularWord(m.result)
        nw in seen && continue
        if isempty(m.result) || reduces_to_empty_down(nw)
            return m
        end
    end
    return nothing
end

# Ring radii, OUTER first. Up to 5 rings share the picture evenly; every further
# ring shrinks the whole inner picture by 20% and adds a new outermost circle.
function _ring_radii(n::Int, R::Real)
    n <= 5 && return [R * (n - i + 1) / n for i in 1:n]
    return vcat(R, 0.8 .* _ring_radii(n - 1, R))
end

# ---- shared: build the batched reduction layers ----------------------------

# Run the batched reduction and return (rings, steps):
#   rings[i]  = raw letter vector on layer i (layer 1 = start word)
#   steps[i]  = (windows, surv, repl_idx) taking layer i to layer i+1 (or to ε)
# `steps` may be shorter than `rings` if the word stalls before ε.
function _reduction_layers(w::CircularWord)
    rings = Vector{Letter}[]
    steps = Tuple[]
    cur = collect(letters(w))
    seen = Set{CircularWord}([CircularWord(cur)])
    while !isempty(cur) && length(rings) < 200
        push!(rings, copy(cur))
        m = _guided_next(cur, seen)
        m === nothing && break
        # one move per layer → a single-window `wins` entry (start, patlen, sign, repl)
        s = m.wpos[1]                                  # window start position (wpos built s..s+pl-1)
        wins = [(s, length(m.wpos), m.dsign, m.repl)]
        surv = m.surv
        repl_idx = [m.ridx]
        push!(steps, (wins, surv, repl_idx))
        push!(seen, CircularWord(m.result))
        cur = m.result
    end
    return rings, steps
end

# ---- radial SVG ------------------------------------------------------------

"""
    radial_svg(w; size = 600, r_outer = 250) -> String

Render the reduction of `w` as a radial-flow SVG string: the start word on the
outermost ring, each ring inward = one PARALLEL batch of reducing moves (as many
disjoint windows as possible fire at once). On every ring the letters sit at
regular angles (n letters → every 2π/n); strands connect each letter to its actual
position on the next ring. ε is not drawn — the picture simply ends where the word
empties. At most 5 rings share the radius evenly; each further ring shrinks the
inner picture by 20%, so the image never grows. Colours 1/2/3 = blue/tomato/green;
a step's node uses the mix colour of its window.
"""
function radial_svg(w::CircularWord; size::Int = 600, r_outer::Real = 250)
    isempty(w) && return _empty_svg(size, "ε")
    rings, steps = _reduction_layers(w)
    isempty(steps) && return _empty_svg(size, compact(w) * " (no reduction)")
    return _radial_core(rings, steps; size = size, r_outer = r_outer)
end

"""
    _radial_core(rings, steps; size, r_outer) -> String

The radial-flow drawing kernel, factored out of `radial_svg` so ANY reduction path
— not just one produced by the guided reducer — can be rendered the same way. A
`WordGraph` / rule side is turned into `(rings, steps)` by `_path_layers` and drawn
here. `rings[i]` = the letter vector on ring i (outer = start); `steps[i] =
(wins, surv, repl_idx)` where `wins = [(start, patlen, dsign, repl)]` is the vertex
between ring i and i+1, `surv :: k=>k2` maps surviving positions, `repl_idx[wi]` are
the produced positions on ring i+1.
"""
function _radial_core(rings::Vector, steps::Vector; size::Int = 600, r_outer::Real = 250)
    cx = cy = size / 2
    nr = length(rings)
    radii = _ring_radii(nr, r_outer)
    ang(k, n) = 2π * (k - 1) / n - π/2
    px(r, θ) = cx + r * cos(θ);  py(r, θ) = cy + r * sin(θ)

    io = IOBuffer()
    println(io, """<svg viewBox="0 0 $size $size" xmlns="http://www.w3.org/2000/svg" font-family="monospace">""")
    for r in radii                                     # guide circles
        println(io, """<circle cx="$cx" cy="$cy" r="$r" fill="none" stroke="#e2e7ee" stroke-dasharray="3" stroke-width="1"/>""")
    end
    for (i, (wins, surv, repl_idx)) in enumerate(steps)
        lc = rings[i]; n = length(lc)
        R1 = radii[i]
        R2 = i < nr ? radii[i + 1] : radii[i] * 0.35   # last step ends in nodes, no ε ring
        m  = i < nr ? length(rings[i + 1]) : 0
        # survivors: strand from ring-i position k to its true position on ring i+1
        for (k, k2) in surv
            θ1 = ang(k, n); θ2 = ang(k2, m); c = lc[k]
            θm = _mean_angle([θ1, θ2]); rm = (R1 + R2) / 2
            println(io, """<path d="M$(px(R1,θ1)),$(py(R1,θ1)) Q$(px(rm,θm)),$(py(rm,θm)) $(px(R2,θ2)),$(py(R2,θ2))" stroke="$(_RGB[c])" stroke-width="2.4" fill="none" opacity="0.9"/>""")
        end
        # windows: strands converge to a node; node feeds its replacement letters.
        # A DOT is a degree-1 vertex: exactly one strand touches it (it ends a strand,
        # pl==1 & no repl, OR it starts one, no window & one repl). We draw a dot as a
        # small filled disc at the STRAND END, not as a merge fan. A trivalent/braid is
        # a merge point where several strands meet and (some) leave.
        for (wi, (s, pl, _, repl)) in enumerate(wins)
            wpos = [mod1(s + j, n) for j in 0:(pl-1)]
            rl = length(repl)
            is_dot = (pl == 1 && rl == 0) || (pl == 0 && rl == 1) ||
                     (pl == 1 && rl == 1)          # a dot sitting on a single strand
            wac = isempty(wpos) ? (isempty(repl_idx[wi]) ? -π/2 : ang(repl_idx[wi][1], m)) :
                                  _mean_angle([ang(k, n) for k in wpos])
            Rn = (R1 + 2R2) / 3
            nodex = px(Rn, wac); nodey = py(Rn, wac)
            if is_dot
                # incoming strand (if any) runs to the dot; outgoing strand (if any)
                # leaves it. The dot itself is a small solid disc = a genuine cap/cup.
                col = isempty(wpos) ? _RGB[repl[1]] : _RGB[lc[wpos[1]]]
                for k in wpos                                   # strand into the dot
                    θ = ang(k, n)
                    println(io, """<path d="M$(px(R1,θ)),$(py(R1,θ)) Q$(px((R1+Rn)/2,θ)),$(py((R1+Rn)/2,θ)) $nodex,$nodey" stroke="$col" stroke-width="2.4" fill="none"/>""")
                end
                for k2 in repl_idx[wi]                          # strand out of the dot
                    θ2 = ang(k2, m)
                    println(io, """<path d="M$nodex,$nodey Q$(px((Rn+R2)/2,θ2)),$(py((Rn+R2)/2,θ2)) $(px(R2,θ2)),$(py(R2,θ2))" stroke="$col" stroke-width="2.4" fill="none"/>""")
                end
                println(io, """<circle cx="$nodex" cy="$nodey" r="5" fill="$col" stroke="white" stroke-width="1"/>""")
                continue
            end
            for k in wpos
                θ = ang(k, n)
                println(io, """<path d="M$(px(R1,θ)),$(py(R1,θ)) Q$(px((R1+Rn)/2,θ)),$(py((R1+Rn)/2,θ)) $nodex,$nodey" stroke="$(_RGB[lc[k]])" stroke-width="2.4" fill="none"/>""")
            end
            for (j, k2) in enumerate(repl_idx[wi])
                θ2 = ang(k2, m)
                println(io, """<path d="M$nodex,$nodey Q$(px((Rn+R2)/2,θ2)),$(py((Rn+R2)/2,θ2)) $(px(R2,θ2)),$(py(R2,θ2))" stroke="$(_RGB[repl[j]])" stroke-width="2.4" fill="none"/>""")
            end
            ncol = _node_rgb(sort(unique(lc[k] for k in wpos)))
            println(io, """<circle cx="$nodex" cy="$nodey" r="6" fill="$ncol"/>""")
        end
        if i == 1                                      # labels on the outer ring only
            for k in 1:n
                θ = ang(k, n); c = lc[k]
                println(io, """<circle cx="$(px(R1,θ))" cy="$(py(R1,θ))" r="9" fill="white" stroke="$(_RGB[c])" stroke-width="2"/>""")
                println(io, """<text x="$(px(R1,θ))" y="$(py(R1,θ)+4)" text-anchor="middle" font-size="11" fill="$(_RGB[c])">$c</text>""")
            end
        end
    end
    # a word that stalls before ε keeps its last ring as plain letters
    if length(steps) < nr
        lc = rings[nr]; n = length(lc); R = radii[nr]
        for k in 1:n
            θ = ang(k, n); c = lc[k]
            println(io, """<circle cx="$(px(R,θ))" cy="$(py(R,θ))" r="5" fill="$(_RGB[c])"/>""")
        end
    end
    println(io, "</svg>")
    return String(take!(io))
end

# ---- shell (concentric) SVG of a WordGraph / MorphismGraph -----------------
#
# `shell_svg` IS the Tutte layout (`_shell_svg_planar` below): the planar
# Tutte solver (`render/tutte/Tutte.jl`) draws every diagram cleanly and is
# already the shared layout for `WordGraph`/`CircularGraph`. There is only one
# planar renderer, reached under two names for the two call sites
# (`display_shell` vs. `display_tutte`, `render/tutte/Display.jl`).
# `_graph_to_path_full`/`_radial_core` stay in use for `radial_svg`/`linear_svg`
# (an actual reduction PATH, a different input shape) and are otherwise unused.

"""
    shell_svg(g; size = 600, r_outer = 250) -> String

Draw a `WordGraph` (or a `MorphismGraph`'s graph) via the planar Tutte layout
(`_shell_svg_planar`, an alias for `tutte_svg`). `r_outer` maps to `tutte_svg`'s
`r`.
"""
function shell_svg(g::WordGraph; size::Int = 600, r_outer::Real = 250)
    isempty(g.word) && isempty(g.nodes) && isempty(g.edges) && return _empty_svg(size, "ε")
    return _shell_svg_planar(g; size = size, r_outer = r_outer)
end

"""
    _shell_svg_planar(g; size = 600, r_outer = 250) -> String

The drawing kernel behind `shell_svg`: delegates straight to `tutte_svg` (the
unified planar solver, `render/tutte/Tutte.jl`).
"""
function _shell_svg_planar(g::WordGraph; size::Int = 600, r_outer::Real = 250)
    return tutte_svg(g; size = size, r = r_outer)
end

# ---- linear (stacked) SVG --------------------------------------------------

"""
    linear_svg(w; width = 640, row_h = 70, pad = 40) -> String

Render the reduction of `w` as a STACKED flow SVG (the alternative to `radial_svg`):
the start word on the BOTTOM row, each row upward = one parallel batch of reducing
moves, ε at the top (not drawn — the picture ends where the word empties). On every
row the letters sit at regular x-positions; strands connect each letter to its true
position on the row above. Same batching, colours, and node/mix logic as the radial
form; here the layout is a left-to-right cut-open circle instead of a ring. Height
grows with the number of rows; width is fixed.
"""
function linear_svg(w::CircularWord; width::Int = 640, row_h::Real = 70, pad::Real = 40)
    isempty(w) && return _empty_svg_wh(width, 120, "ε")
    rings, steps = _reduction_layers(w)
    isempty(steps) && return _empty_svg_wh(width, 120, compact(w) * " (no reduction)")

    nr = length(rings)
    height = round(Int, 2pad + (nr - 1) * row_h + 40)
    xpos(k, n) = pad + (n == 1 ? (width - 2pad) / 2 : (k - 1) * (width - 2pad) / (n - 1))
    yrow(i) = height - pad - (i - 1) * row_h          # row 1 = bottom (start word)

    io = IOBuffer()
    println(io, """<svg viewBox="0 0 $width $height" xmlns="http://www.w3.org/2000/svg" font-family="monospace">""")
    for i in 1:nr                                      # faint guide baselines
        y = yrow(i)
        println(io, """<line x1="$pad" y1="$y" x2="$(width-pad)" y2="$y" stroke="#e2e7ee" stroke-dasharray="3" stroke-width="1"/>""")
    end
    for (i, (wins, surv, repl_idx)) in enumerate(steps)
        lc = rings[i]; n = length(lc)
        y1 = yrow(i)
        y2 = i < nr ? yrow(i + 1) : yrow(i) - row_h * 0.5
        m  = i < nr ? length(rings[i + 1]) : 0
        ym = (y1 + y2) / 2
        for (k, k2) in surv                            # survivor strands
            x1 = xpos(k, n); x2 = xpos(k2, m); c = lc[k]
            println(io, """<path d="M$x1,$y1 C$x1,$ym $x2,$ym $x2,$y2" stroke="$(_RGB[c])" stroke-width="2.4" fill="none" opacity="0.9"/>""")
        end
        for (wi, (s, pl, _, repl)) in enumerate(wins)  # window → node → replacement
            wpos = [mod1(s + j, n) for j in 0:(pl-1)]
            nodex = sum(xpos(k, n) for k in wpos) / length(wpos)
            nodey = y1 + (y2 - y1) * 0.55
            for k in wpos
                x1 = xpos(k, n)
                println(io, """<path d="M$x1,$y1 C$x1,$((y1+nodey)/2) $nodex,$((y1+nodey)/2) $nodex,$nodey" stroke="$(_RGB[lc[k]])" stroke-width="2.4" fill="none"/>""")
            end
            for (j, k2) in enumerate(repl_idx[wi])
                x2 = xpos(k2, m)
                println(io, """<path d="M$nodex,$nodey C$nodex,$((nodey+y2)/2) $x2,$((nodey+y2)/2) $x2,$y2" stroke="$(_RGB[repl[j]])" stroke-width="2.4" fill="none"/>""")
            end
            ncol = _node_rgb(sort(unique(lc[k] for k in wpos)))
            println(io, """<circle cx="$nodex" cy="$nodey" r="6" fill="$ncol"/>""")
        end
        if i == 1                                      # labels on the bottom row
            for k in 1:n
                x = xpos(k, n); c = lc[k]
                println(io, """<circle cx="$x" cy="$y1" r="9" fill="white" stroke="$(_RGB[c])" stroke-width="2"/>""")
                println(io, """<text x="$x" y="$(y1+4)" text-anchor="middle" font-size="11" fill="$(_RGB[c])">$c</text>""")
            end
        end
    end
    if length(steps) < nr                              # stalled: draw last row plain
        lc = rings[nr]; n = length(lc); y = yrow(nr)
        for k in 1:n
            println(io, """<circle cx="$(xpos(k,n))" cy="$y" r="5" fill="$(_RGB[lc[k]])"/>""")
        end
    end
    println(io, "</svg>")
    return String(take!(io))
end

# mean of angles (circular mean, so window arcs near the -π/π seam average right)
function _mean_angle(θs)
    isempty(θs) && return 0.0
    return atan(sum(sin, θs) / length(θs), sum(cos, θs) / length(θs))
end

function _empty_svg(size, label)
    c = size / 2
    return """<svg viewBox="0 0 $size $size" xmlns="http://www.w3.org/2000/svg" font-family="monospace">
<text x="$c" y="$c" text-anchor="middle" font-size="16" fill="#888">$label (no reduction)</text></svg>"""
end

function _empty_svg_wh(width, height, label)
    return """<svg viewBox="0 0 $width $height" xmlns="http://www.w3.org/2000/svg" font-family="monospace">
<text x="$(width/2)" y="$(height/2)" text-anchor="middle" font-size="16" fill="#888">$label</text></svg>"""
end

"""
    save_radial(path, w; kwargs...)

Render `radial_svg(w; kwargs...)` and write it to `path` (a `.svg` file).
Convention: rendered images live under `data/svg/`, not the `data/` top level.
"""
function save_radial(path::AbstractString, w::CircularWord; kwargs...)
    mkpath(dirname(path))
    write(path, radial_svg(w; kwargs...))
    return path
end

"""
    save_linear(path, w; kwargs...)

Render `linear_svg(w; kwargs...)` and write it to `path` (a `.svg` file).
Convention: rendered images live under `data/svg/`, not the `data/` top level.
"""
function save_linear(path::AbstractString, w::CircularWord; kwargs...)
    mkpath(dirname(path))
    write(path, linear_svg(w; kwargs...))
    return path
end
