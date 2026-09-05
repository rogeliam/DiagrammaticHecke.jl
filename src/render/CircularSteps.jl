# render/CircularSteps.jl — the reduction as an IMAGE EQUATION:  D  →  c₁·D₁ + c₂·D₂ + …
#
# ---- WHAT THE TRACER COVERS ----------------------------------------------------
#
# Tracing `CIRCULAR_RULES` alone is not enough: the interesting steps of the
# decorated driver `reduce_to_circular_leave` are exactly
# the ones NOT in `CIRCULAR_RULES` (P1, 2parallel, fusion, D4). A step
# here is therefore generally
#
#     (name::Symbol, f)   with   f(fdm::CircularDecoratedMorphism) -> Nothing | CircularComboR
#
# and [`CIRCULAR_STEPS`](@ref) is the list of these steps IN THE DRIVER'S ORDER. A
# subset is obtained with [`circular_steps`](@ref), e.g. `circular_steps(:circular_rules)` for
# only the base rules; to check a single rule, build a one-element list.
#
# ---- TERMINATION WARNING -------------------------------------------------------
#
# `CIRCULAR_STEPS` contains P1 and 2parallel, which INCREASE `circular_weight` (see the
# headers of CircularParallelMerge.jl / Circular2Parallel.jl). The tracer is therefore
# NOT guaranteed to terminate — it is a display tool, not a driver, and always
# runs against `maxsteps`. For a guaranteed finite sequence, use
# `steps = circular_steps(:circular_rules)` — the `circular_weight` assert applies there.

"""
    CIRCULAR_STEPS :: Vector{Tuple{Symbol, Any}}

The reduction steps in the order `reduce_to_circular_leave` tries them — each as
`(name, f)` with `f(fdm::CircularDecoratedMorphism) -> Union{Nothing, CircularComboR}`:

| name | what |
|---|---|
| `:scalars` | [`circular_extract_scalars`](@ref) — a label at distance 0 moves into the coefficient |
| `:fusion` | [`circular_fusion_step`](@ref) — polynomial over an edge to the marking |
| `:circular_rules` | the first matching rule from [`CIRCULAR_RULES`](@ref) |
| `:parallel_merge` | [`circular_parallel_merge_step`](@ref) — P1, parallel 1/3 edges |
| `Symbol("2parallel")` | [`circular_2parallel_step`](@ref) — 2parallel, parallel same-colour edges |
| `:zamo` | [`circular_zamo_step`](@ref) — Z1/Z2 on the exact morphism key |
| `:d4` | D4 via [`find_circular_d4_any`](@ref), both variants |

THE ORDER IS NOT ARBITRARY: the polynomial steps stand BEFORE `:circular_rules`,
exactly as in `reduce_to_circular_leave`. Otherwise the tracer could try a graph
rule while a polynomial still sits on a boundary-gap-less inner region — the
rule destroys the region and `_circular_transfer_labels` throws. The driver never
hits this because it always runs `circular_extract_scalars`/`circular_fusion_step`
before `reduce_circular`.

Tracer and driver agree on running order: `circular_step_once` runs this list from
the top on EVERY step, checking the polynomial steps each round — and
`reduce_to_circular_leave` does the same.

There are two polynomial steps, `:scalars` and `:fusion`.

A subset (in the same order) is given by [`circular_steps`](@ref). The name of the
2parallel rule is `Symbol("2parallel")` — a literal `:2parallel` wouldn't work
(Julia symbols can't start with a digit); output simply shows `2parallel`.

Meant for [`circular_trace`](@ref) and [`circular_step_once`](@ref). Order and contents
are a DISPLAY of the driver, not a second source of truth — if
`reduce_to_circular_leave` changes, this list must be updated to match.
"""
const CIRCULAR_STEPS = Tuple{Symbol, Any}[
    (:scalars, function (fdm::CircularDecoratedMorphism)
         factor, rest = circular_extract_scalars(fdm)
         isone(factor) && return nothing
         return factor * CircularComboR(_circular_decorated(rest))
     end),
    (:fusion,         circular_fusion_step),
    (:circular_rules, function (fdm::CircularDecoratedMorphism)
         hit = _circular_first_rule_match(_circular_decorated(fdm), CIRCULAR_RULES)
         return hit === nothing ? nothing : hit[2]
     end),
    (:parallel_merge, circular_parallel_merge_step),
    (Symbol("2parallel"), circular_2parallel_step),
    (:zamo, function (fdm::CircularDecoratedMorphism)
         return circular_zamo_step(_circular_decorated(fdm), fdm.m.cut1)
     end),
    (:d4, function (fdm::CircularDecoratedMorphism)
         h = find_circular_d4_any(fdm.m)
         h === nothing && return nothing
         return _apply_circular_d4_hit(_circular_decorated(fdm), h)
     end),
    # The two `:after_d4` stages of `reduce_to_circular_leave`, in its order.
    # Both decline unless the outer label is trivial, exactly as the driver
    # demands, and both follow their switch. The driver additionally suppresses
    # them below a term that has picked up an α (`rex_alpha_blocked`); that is a
    # property of the RECURSION, not of a single step, so a one-step display
    # cannot carry it — a trace may therefore show a step the full driver would
    # have skipped further down.
    (:rex_fusion, function (fdm::CircularDecoratedMorphism)
         (CIRCULAR_REX_FUSION_ENABLED[] && isone(fdm.outer_label)) || return nothing
         r = circular_rex_fusion_step(fdm)
         return r === nothing ? nothing : r.combo
     end),
    (:dot_fusion, function (fdm::CircularDecoratedMorphism)
         (CIRCULAR_DOT_FUSION_ENABLED[] && isone(fdm.outer_label)) || return nothing
         r = circular_dot_fusion_step(fdm)
         return r === nothing ? nothing : r.combo
     end),
]

"""
    circular_steps(names::Symbol...) -> Vector{Tuple{Symbol, Any}}

The sub-list of [`CIRCULAR_STEPS`](@ref) matching the given names, in `CIRCULAR_STEPS`'s
order (not the order of the arguments). Unknown names are an error — a typo
should not silently produce an empty step list.

```julia
circular_trace(g, m; steps = circular_steps(:circular_rules))              # base rules only
circular_trace(g, m; steps = circular_steps(:scalars, :fusion, :d4))  # the D4 path
```
"""
function circular_steps(names::Symbol...)
    known = Set(first.(CIRCULAR_STEPS))
    unknown = setdiff(Set(names), known)
    isempty(unknown) || error("circular_steps: unknown steps $(sort(collect(unknown))) — " *
                                "known are $(sort(collect(known)))")
    return [s for s in CIRCULAR_STEPS if s[1] in names]
end

"""
    circular_step_once(fd::CircularDecorated, m0; steps = CIRCULAR_STEPS)
        -> Union{Nothing, NamedTuple}

A single reduction step on ONE diagram: the first step from `steps` that
applies. Returns `(name, before, after)` with `after::CircularComboR` (the
right-hand side of the image equation, possibly with several terms), or
`nothing` if no step applies.

`m0` only supplies the cuts (marking); the graph comes from `fd`.
"""
function circular_step_once(fd::CircularDecorated, m0; steps = CIRCULAR_STEPS)
    fdm = CircularDecoratedMorphism(CircularMorphismGraph(fd.graph, m0.cut1, m0.cut2),
                               fd.region_labels, fd.outer_label)
    for (name, f) in steps
        res = f(fdm)
        res === nothing && continue
        return (name = name, before = fd, after = res)
    end
    return nothing
end

# ---- Drawing --------------------------------------------------------------------
#
# PNG instead of SVG in a flex box: VS Code skips `image/svg+xml`, and only
# this way do the terms actually stand SIDE BY SIDE instead of stacked.

_circular_tile_png(g::CircularGraph, m0) = Base64.base64encode(_svg_to_png(
    _circular_tutte_morphism(CircularMorphismGraph(g, m0.cut1, m0.cut2);
                        cell_numbers = true, cell_labels = :labels_dist,
                        leaf_numbers = true, slot_numbers = true)))

_circular_tile(g::CircularGraph, m0, width, caption) = string(
    "<div style='margin:0 6px 10px 0;text-align:center;font-family:monospace;",
    "font-size:12px'><img src='data:image/png;base64,", _circular_tile_png(g, m0),
    "' style='width:", width, "px'><div>", caption, "</div></div>")

_circular_glyph(s) = string("<div style='align-self:center;font-family:monospace;",
                       "font-size:26px;margin:0 8px 10px 0'>", s, "</div>")

# Short caption for a term: node count, braids, and (if present) the
# non-trivial region labels — what you actually need when comparing.
function _circular_tile_caption(fd::CircularDecorated, coeff = nothing)
    g = fd.graph
    parts = String[]
    coeff !== nothing && !isone(coeff) && push!(parts, string(coeff, " ·"))
    push!(parts, string(length(g.nodes), " nodes"))
    nb = count(n -> n.kind === :braid && arm_count(n) == 6, g.nodes)
    nb > 0 && push!(parts, string(nb, " Braid", nb == 1 ? "" : "s"))
    lab = [string(x) for x in fd.region_labels if !isone(x)]
    isempty(lab) || push!(parts, join(lab, ", "))
    return join(parts, "  ")
end

"""
    show_circular_step_row(before, after, m0; width = 260, title = "")

Draws the **image equation** `before → c₁·D₁ + c₂·D₂ + …` as a row of
side-by-side diagrams (PNG in an HTML flex box, so VS Code shows them).
`before` is a `CircularDecorated` or `CircularGraph`, `after` a `CircularComboR` (an empty
sum is drawn as `0`). `m0` supplies the cuts.

Below each image: node count, braid count and, if present, the non-trivial
region labels; above the row the `title` (e.g. the rule name). Without a
display stack (script, REPL without display) nothing happens — the same
function works in both contexts.
"""
function show_circular_step_row(before, after::CircularComboR, m0;
                           width::Int = 260, title::AbstractString = "")
    vfd = before isa CircularGraph ? circular_decorated(before) : before
    parts = String[_circular_tile(vfd.graph, m0, width, _circular_tile_caption(vfd)),
                   _circular_glyph("&rarr;")]
    terms = sort(pairs_of(after); by = p -> length(p[1].graph.nodes))
    isempty(terms) && push!(parts, _circular_glyph("0"))
    for (i, (dd, co)) in enumerate(terms)
        i > 1 && push!(parts, _circular_glyph("+"))
        push!(parts, _circular_tile(dd.graph, m0, width, _circular_tile_caption(dd, co)))
    end
    html = string("<div style='font-family:monospace;font-size:13px'>", title,
                  "</div><div style='display:flex;flex-wrap:wrap;",
                  "align-items:flex-start'>", join(parts), "</div>")
    try
        display(MIME("text/html"), html)
    catch
        # no display stack (script/test) — the row is only reported then.
        isempty(title) || println(title)
    end
    return nothing
end

"""
    circular_trace(fd, m0; steps = CIRCULAR_STEPS, maxsteps = 40, draw = true, width = 260)
        -> CircularComboR

Reduces `fd` (a `CircularDecorated`, `CircularGraph` or `CircularComboR`) step by step and
draws **every step as an image equation** `D_vorher → c₁·D₁ + …`
([`show_circular_step_row`](@ref)). Each round takes the FIRST term a step from
`steps` applies to; the remaining terms carry over unchanged. Returns the
combination at the fixed point, or after `maxsteps`.

`draw = false` prints only the rule names — useful for long runs and for tests.

⚠️ `CIRCULAR_STEPS` contains P1/2parallel, which increase `circular_weight`; the tracer
is therefore not guaranteed to terminate and always runs against `maxsteps`
(see the module header). With `steps = CIRCULAR_STEPS[1:1]` you stay within the
base rules, where the `circular_weight` assert applies.

After each step, `euler`/`check_wiring` of all new terms are checked and
violations reported at the end — the house sanity check thus runs automatically.
"""
function circular_trace(fd, m0; steps = CIRCULAR_STEPS, maxsteps::Int = 40,
                   draw::Bool = true, width::Int = 260)
    cur = fd isa CircularComboR ? fd :
          CircularComboR(fd isa CircularGraph ? circular_decorated(fd) : fd)
    hist = Symbol[]
    impure = 0
    for k in 1:maxsteps
        out = CircularComboR()
        fired = nothing
        for (dd, coeff) in pairs_of(cur)
            if fired === nothing
                st = circular_step_once(dd, m0; steps = steps)
                if st === nothing
                    out = out + coeff * CircularComboR(dd)
                else
                    fired = st
                    out = out + coeff * st.after
                end
            else
                out = out + coeff * CircularComboR(dd)
            end
        end
        fired === nothing && break
        push!(hist, fired.name)
        impure += sum(length(check_wiring(dd.graph))
                      for (dd, _) in pairs_of(fired.after); init = 0)
        impure += count(dd -> euler(dd.graph) != 2,
                        [dd for (dd, _) in pairs_of(fired.after)])
        title = string("Step ", k, ":  :", fired.name)
        draw ? show_circular_step_row(fired.before, fired.after, m0;
                                 width = width, title = title) : println(title)
        cur = out
    end
    println("FIXED POINT after ", length(hist), " steps:  ", length(pairs_of(cur)),
            " terms, node counts ",
            sort([length(dd.graph.nodes) for (dd, _) in pairs_of(cur)]))
    println("Step sequence: ", hist)
    println("Violations (wiring/Euler) across all intermediate results: ", impure)
    return cur
end

"""
    show_combo(c::CircularComboR; expected = nothing, m = nothing)

Draws every term of a sum — the shared helper so notebooks don't each define
their own. Per term: a coefficient line (coefficient · degree via
[`circular_term_degree`](@ref), plus a degree check against `expected`, the
non-trivial labels and the node list), then `display` of the term — decorated;
with a marking when `m::CircularMorphismGraph` is passed (`with_marks`). Empty
sum: `0 (empty sum)`.

A single `CircularDecorated`/`CircularGraph` is also accepted (a one-term sum). A
notebook that defines its own `show_combo` simply shadows this
version — both work.
"""
function show_combo(c; expected = nothing, m = nothing)
    combo = c isa CircularComboR ? c :
            CircularComboR(c isa CircularGraph ? circular_decorated(c) : c)
    ps = collect(pairs_of(combo))
    isempty(ps) && (println("0 (empty sum)"); return nothing)
    for (d, co) in ps
        deg  = circular_term_degree(co, d)
        mark = expected === nothing ? "" :
                (deg == expected ? "   ✔ degree ok" : "   ✘ DEGREE $(deg) INSTEAD OF $(expected)")
        lab = [(i, string(x)) for (i, x) in enumerate(d.region_labels) if !isone(x)]
        println("Coefficient ", co, "   degree ", deg, mark)
        println("   Labels: ", isempty(lab) ? "—" : lab,
                "   outer: ", string(d.outer_label),
                "   nodes: ", [(n.kind, n.arms) for n in d.graph.nodes])
        _circular_warn_not_normal(d, m)
        display(m === nothing ? d : with_marks(d, m))
    end
    return nothing
end

"""
    CIRCULAR_WARN_NOT_NORMAL

Switch for the normal-form hint when printing. Default `true`. Turn off if the
warning floods output: `CIRCULAR_WARN_NOT_NORMAL[] = false`.
"""
const CIRCULAR_WARN_NOT_NORMAL = Ref(true)

"""
    _circular_warn_not_normal(d, m) -> nothing

Reports dots whose distance word can end on the dot colour via braid moves
([`circular_dot_reducible`](@ref)) — then the term is NOT in normal form and can be
reduced further.

Needs a marking: without `m` there is no distance and thus no distance word,
so nothing is reported.
"""
function _circular_warn_not_normal(d, m)
    CIRCULAR_WARN_NOT_NORMAL[] || return nothing
    m === nothing && return nothing
    fm = with_marks(d, m)
    mg = fm isa CircularDecoratedMorphism ? fm.m :
         (fm isa CircularMorphismGraph ? fm : nothing)
    mg === nothing && return nothing
    bad = try
        circular_non_normal_dots(mg)
    catch
        return nothing          # never let this hint fail the caller
    end
    isempty(bad) && return nothing
    for v in bad
        w = circular_distance_word(mg, v)
        s = arm_colour(mg.graph.nodes[v], 1)
        rotations = [u for u in circular_braid_class(w) if !isempty(u) && u[end] == s]
        println("   ⚠ NOT IN NORMAL FORM: dot ", v, " (colour ", s,
                "), distance word ", join(w),
                " can end on ", s, " (e.g. ", join(first(sort(rotations))), ")")
    end
    return nothing
end

