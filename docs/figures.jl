# docs/figures.jl — render every figure of the manual to `docs/src/assets/figures/`.
#
#     julia --project=docs docs/figures.jl
#
# The manual is static HTML, so the pictures are static too: this script draws
# them once with the package itself and writes plain SVG files, plus the two
# generated pages `rules-gallery.md` and `gallery.md`. Run it before
# `docs/make.jl` whenever a rule, a fixture or a renderer changes.
using DiagrammaticHecke
const DH = DiagrammaticHecke

const OUT = joinpath(@__DIR__, "src", "assets", "figures")
mkpath(OUT)

debug_circular(true)          # arm, region and leaf numbers in every picture

# ---- 1. Composing SVGs -----------------------------------------------------
#
# A figure is a row: diagrams side by side with `→`, `+` or `=` between them.
# The tiles are nested `<svg>` elements, which keep their own viewBox, so no
# coordinate in the renderers has to be touched.

const TILE = 300          # px per diagram in a row
const GLYPH = 46          # px per separator column
const CAPTION = 20        # px of caption strip under a tile

"Strip the size attributes of the outermost `<svg>` and place it at `(x, y)`."
function _nest(svg::AbstractString, x::Real, y::Real, w::Real)
    head_end = findfirst('>', svg)
    head, body = svg[1:head_end], svg[head_end+1:end]
    head = replace(head, r"\s(width|height)=\"[^\"]*\"" => "")
    head = replace(head, "<svg" => "<svg x=\"$x\" y=\"$y\" width=\"$w\" height=\"$w\"")
    return head * body
end

_text(x, y, s; size = 14, anchor = "middle") = string(
    "<text x=\"$x\" y=\"$y\" text-anchor=\"$anchor\" font-size=\"$size\" ",
    "font-family=\"monospace\" fill=\"#888\">", s, "</text>")

"""
    svg_row(items; file, tile = TILE) -> String

`items` is a list of `svg::String`, `(svg, caption)` or `glyph::Symbol`
(`:arrow`, `:plus`, `:eq`, `:zero`). Writes `<file>.svg` under the figure
directory and returns the path, relative to `docs/src`.
"""
function svg_row(items::Vector; file::AbstractString, tile::Int = TILE)
    xs, w = Float64[], 0.0
    for it in items
        push!(xs, w)
        w += it isa Symbol ? GLYPH : tile
    end
    h = tile + CAPTION
    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $w $h\" ",
              "width=\"$(round(Int, w))\" height=\"$(round(Int, h))\">")
    for (x, it) in zip(xs, items)
        if it isa Symbol
            g = it === :arrow ? "&#8594;" : it === :plus ? "+" : it === :eq ? "=" : "0"
            print(io, _text(x + GLYPH/2, tile/2 + 8, g; size = 24))
        else
            svg, cap = it isa Tuple ? it : (it, "")
            print(io, _nest(svg, x, 0, tile))
            isempty(cap) || print(io, _text(x + tile/2, tile + 14, cap; size = 12))
        end
    end
    print(io, "</svg>")
    path = joinpath(OUT, file * ".svg")
    tmp = path * ".tmp"
    write(tmp, String(take!(io)))
    mv(tmp, path; force = true)
    return "assets/figures/" * file * ".svg"
end

svg_one(svg; file, tile = TILE) = svg_row(Any[svg]; file = file, tile = tile)

# ---- 2. Tiles --------------------------------------------------------------

circ_tile(g::CircularGraph, m0) = DH._circular_tutte_morphism(
    CircularMorphismGraph(g, m0.cut1, m0.cut2);
    cell_numbers = true, cell_labels = :labels_dist,
    leaf_numbers = true, slot_numbers = true)
circ_tile(g::CircularGraph) = circ_tile(g, CircularMorphismGraph(g, 0, 0))
circ_tile(fm::CircularMorphismGraph) = circ_tile(fm.graph, fm)
circ_tile(fd::CircularDecorated, m0) = circ_tile(fd.graph, m0)

word_tile(g::WordGraph) = tutte_svg(g)
word_tile(m::MorphismGraph) = DH._tutte_morphism(m)

caption(fd::CircularDecorated, coeff = nothing) = DH._circular_tile_caption(fd, coeff)

"""
    step_row(before, after, m0; file, reversed = false)

The image equation `before → c₁·D₁ + …` as one SVG file. `reversed` draws the
same equation the other way round, `c₁·D₁ + … → before`, for the rules whose
natural reading is the one the driver applies backwards.
"""
function step_row(before, after::CircularComboR, m0; file, reversed::Bool = false)
    fd = before isa CircularGraph ? circular_decorated(before) : before
    lhs = Any[(circ_tile(fd.graph, m0), caption(fd))]
    rhs = Any[]
    terms = sort(pairs_of(after); by = p -> length(p[1].graph.nodes))
    isempty(terms) && push!(rhs, :zero)
    for (i, (dd, co)) in enumerate(terms)
        i > 1 && push!(rhs, :plus)
        push!(rhs, (circ_tile(dd.graph, m0), caption(dd, co)))
    end
    items = reversed ? vcat(rhs, Any[:arrow], lhs) : vcat(lhs, Any[:arrow], rhs)
    return svg_row(items; file = file)
end

# ---- 3. The figures of the guide pages -------------------------------------

println("figures: guide"); flush(stdout)

const ll = light_leaf([1, 2, 1], [1, 0, 1])
const fll = circular_morphism(ll)

svg_one(DH._tutte_morphism(ll; leaf_numbers = true); file = "light-leaf-plain")
svg_one(circ_tile(fll); file = "light-leaf-circular")
svg_row(Any[(DH._tutte_morphism(ll; leaf_numbers = true), "MorphismGraph"), :eq,
            (circ_tile(fll), "CircularGraph")]; file = "light-leaf-both")

# The region word and the tree it is read off. `:words_num` prints `region|word`.
svg_one(DH._circular_tutte_morphism(fll; cell_numbers = true, cell_labels = :words_num,
                                    leaf_numbers = true, region_tree = true);
        file = "region-words")

# Named plain fixtures — the diagrams the dot rules need.
svg_row(Any[(tutte_svg(unit_graph(); leaf_numbers = true), "unit_graph"),
            (tutte_svg(needle_a_graph(); leaf_numbers = true), "needle_a_graph"),
            (tutte_svg(needle_b_graph(); leaf_numbers = true), "needle_b_graph"),
            (tutte_svg(barbell_graph(); leaf_numbers = true), "barbell_graph")];
        file = "fixtures")

# Double leaves of low degree over 121 -> 121.
let ds = [d for d in double_leaves([1, 2, 1], [1, 2, 1]) if d.degree <= 1]
    items = Any[]
    for d in ds[1:min(4, length(ds))]
        isempty(items) || push!(items, :plus)
        push!(items, (DH._tutte_morphism(d.morphism; leaf_numbers = true),
                      "DL($(join(d.e))|$(join(d.f)))  deg $(d.degree)"))
    end
    svg_row(items; file = "double-leaves")
end

# ---- 4. A reduction, step by step ------------------------------------------

println("figures: trace"); flush(stdout)

"""
    trace_figures(fm, stem; maxsteps) -> Vector{Tuple{String,String}}

Walks the reduction of `fm` with `circular_step_once`, writes one row per step
and returns `(rule name, figure path)`. A multi-term step ends the walk: in a
linear trace the picture of a sum is the last thing there is to show.
"""
function trace_figures(fm, stem::AbstractString; maxsteps::Int = 8)
    rows = Tuple{String,String}[]
    cur = circular_decorated(fm.graph)
    for k in 1:maxsteps
        st = circular_step_once(cur, fm)
        st === nothing && break
        push!(rows, (string(st.name),
                     step_row(st.before, st.after, fm; file = "$stem-$k")))
        terms = collect(pairs_of(st.after))
        length(terms) == 1 || break
        cur = terms[1][1]
    end
    return rows
end

"""
    longest_trace(pairs; least) -> (rows, entry)

The double leaf, over the given word pairs, whose reduction takes the most
steps. The search stops at `least` steps: a longer trace makes a longer page,
not a better example. Degree is capped so the walk stays over short diagrams.
"""
function longest_trace(pairs; least::Int = 3, degree_max::Int = 2)
    best, best_rows, best_d = 0, Tuple{String,String}[], nothing
    for (x, y) in pairs, d in double_leaves(x, y)
        d.degree <= degree_max || continue
        rows = trace_figures(circular_morphism(d.morphism), "trace")
        length(rows) > best && ((best, best_rows, best_d) = (length(rows), rows, d))
        best >= least && break
    end
    best_d === nothing && return best_rows, best_d
    println("  trace: DL(", join(best_d.e), "|", join(best_d.f), ") in ", best, " step(s)")
    # The search overwrites `trace-k.svg` for every candidate it walks, so the
    # winner is drawn once more, last.
    return trace_figures(circular_morphism(best_d.morphism), "trace"), best_d
end

const TRACE, TRACE_DL = longest_trace([([1, 2, 1], [1, 2, 1]), ([1, 2, 1, 2], [1, 2, 1, 2]),
                                       ([2, 1, 3, 2], [2, 1, 3, 2])]; least = 4)

# ---- 5. Zamolodchikov ------------------------------------------------------

println("figures: zamolodchikov"); flush(stdout)

const ZA = Zamo(1, 8)
const ZB = flip(Zamo(8, 1))
svg_row(Any[(DH._tutte_morphism(ZA; leaf_numbers = true), "Zamo(1,8)"), :eq,
            (DH._tutte_morphism(ZB; leaf_numbers = true), "flip(Zamo(8,1))")];
        file = "zamolodchikov")
const ZTRACE = trace_figures(circular_morphism(ZB), "zamo-trace"; maxsteps = 4)

# ---- 6. Every rule with one drawn example ----------------------------------

println("figures: rules"); flush(stdout)

# Harvesting the rule examples is the expensive part of this script. With
# `DH_FIGURES=nourules` it is skipped and the rows are read back off disk, which
# is what one wants while editing the prose of the pages.
const SKIP_RULES = get(ENV, "DH_FIGURES", "") == "nourules"

"First diagram each rule of `rules` fires on, harvested from double leaves."
function rule_examples(pairs; rules = CIRCULAR_RULES, degree_max = 2, extra = [])
    found = Dict{Symbol,Tuple{CircularDecorated,CircularComboR,Any}}()
    inputs = Any[circular_morphism(d.morphism) for (x, y) in pairs
                 for d in double_leaves(x, y) if d.degree <= degree_max]
    append!(inputs, extra)
    for fm in inputs
        combo = CircularComboR(circular_decorated(fm.graph))
        for _ in 1:200
            next = CircularComboR(); fired = false
            for (dd, coeff) in pairs_of(combo)
                hit = nothing
                for r in rules
                    res = r.apply(dd)
                    res === nothing && continue
                    haskey(found, r.name) || (found[r.name] = (dd, res, fm))
                    hit = res; break
                end
                next = next + coeff * (hit === nothing ? CircularComboR(dd) : hit)
                hit === nothing || (fired = true)
            end
            combo = next
            fired || break
        end
    end
    return found
end

"The graph `g` with a dot on boundary leaf `k`."
function with_dot(g::CircularGraph, k::Int)
    lets = letters(g.word)
    nodes = vcat(copy(g.nodes), [circular_node([lets[k]])])
    dot = length(nodes)
    shift(p) = p isa Leaf ? (p.k == k ? NodePort(dot, 1) : Leaf(p.k > k ? p.k - 1 : p.k)) : p
    CircularGraph(CircularWord(vcat(lets[1:k-1], lets[k+1:end])), nodes,
                  Edge[Edge(e.colour, shift(e.a), shift(e.b)) for e in g.edges])
end

const PAIRS = [([1, 2, 1], [1, 2, 1]), ([1, 2, 1, 3], [1, 2, 1, 3]),
               ([2, 1, 3, 2], [2, 1, 3, 2]), ([1, 2, 1, 3, 2], [1, 2, 1, 3, 2]),
               ([1, 2, 1, 3, 2, 1], [1, 2, 1, 3, 2, 1])]

examples = SKIP_RULES ? Dict{Symbol,Tuple{CircularDecorated,CircularComboR,Any}}() :
           rule_examples(PAIRS)

# The dot rules need a dot next to a trivalent or on a bigon: named fixtures.
SKIP_RULES || merge!(examples, filter(p -> !haskey(examples, p.first),
    rule_examples([]; extra = [CircularMorphismGraph(circular(g), 0, 0)
                               for g in (unit_graph(), needle_a_graph(), needle_b_graph(),
                                         barbell_graph(), r3_braid_dot(), r4_braid_dot())])))

# The big-node family starts from two braid nodes sharing three edges.
SKIP_RULES || let cluster = general_braid_cluster(2)
    c13, _ = reduce_circular(cluster)
    g8 = first(pairs_of(c13))[1].graph                    # the 8-armed node
    merge!(examples, filter(p -> !haskey(examples, p.first),
        rule_examples([]; extra = [CircularMorphismGraph(g, 0, 0)
                                   for g in (cluster, with_dot(with_dot(g8, 1), 1))])))
    # A single dot on an 8-armed node is left alone by the dot policy, so C14
    # has to be called directly.
    if !haskey(examples, :dot_on_gen12)
        gd = with_dot(g8, 1)
        f = DH._lift(g -> DH._fr_dot_on_gen12(g; policy = false))
        examples[:dot_on_gen12] = (circular_decorated(gd), f(circular_decorated(gd)),
                                   CircularMorphismGraph(gd, 0, 0))
    end
end

# Seven rules need a many-armed node, a particular double connection or an
# unfolded preimage; `CircularRuleFixtures.jl` builds a witness for each.
for (nm, g) in (SKIP_RULES ? [] : [(:bead,                circular_bead_graph()),
                (:commutation_merge,   circular_counter_joined([1,3,1,3,1,3], [3,1,3,1,3,1])),
                (:braid_on_gen12,      circular_glued_braids(8, 6, 3)),
                (:gen12_on_gen12,      circular_glued_braids(8, 8, 3)),
                (:gen12_on_gen12_null, circular_glued_braids(8, 6, 4)),
                (:gen12_two_edges,     circular_two_edge_fixture()),
                (:mixed_into_gen12,    circular_mixed_at_braid(3, (1,2), [1,3,1,1,3]))])
    haskey(examples, nm) && continue
    g === nothing && continue
    r = only(r for r in CIRCULAR_RULES if r.name === nm)
    before = circular_decorated(g)
    after = r.apply(before)
    after === nothing && continue
    examples[nm] = (before, after, CircularMorphismGraph(g, 0, 0))
end

# ---- 6b. Better pictures for a few rules -----------------------------------
#
# The general harvest keeps the FIRST diagram a rule fires on, which is not
# always the clearest one: it can carry a third colour the rule has nothing to
# do with, or twice the nodes needed. For the rules below the search is run
# again, for that rule alone, over inputs chosen for the picture, and the
# SMALLEST match wins.

"Every diagram along the reductions of `inputs` on which `r` fires."
function rule_hits(r, inputs)
    out = Tuple{CircularDecorated,CircularComboR,Any}[]
    for fm in inputs
        combo = CircularComboR(circular_decorated(fm.graph))
        for _ in 1:200
            next = CircularComboR(); fired = false
            for (dd, coeff) in pairs_of(combo)
                res = r.apply(dd)
                res === nothing || push!(out, (dd, res, fm))
                hit = res
                if hit === nothing
                    for rr in CIRCULAR_RULES
                        z = rr.apply(dd)
                        z === nothing || (hit = z; break)
                    end
                end
                next = next + coeff * (hit === nothing ? CircularComboR(dd) : hit)
                hit === nothing || (fired = true)
            end
            combo = next
            fired || break
        end
    end
    return out
end

morphisms_of(pairs; degree_max = 2) =
    Any[circular_morphism(d.morphism) for (x, y) in pairs
        for d in double_leaves(x, y) if d.degree <= degree_max]

# `121` and `1212` are the two-colour words: no diagram over them can contain a
# 3, which is what the braid rules should be shown on.
const TWO_COLOUR = [([1, 2, 1], [1, 2, 1]), ([1, 2, 1, 2], [1, 2, 1, 2])]

const PREFERRED = Dict{Symbol,Any}(
    :pitchfork      => () -> morphisms_of(TWO_COLOUR),
    :braid_back     => () -> morphisms_of(TWO_COLOUR),
    :braid_relation => () -> morphisms_of(TWO_COLOUR),
    # The needle IS the named fixture: a strand that runs back into its own
    # node. Nothing shows the self-connection more plainly.
    :needle         => () -> [CircularMorphismGraph(circular(g), 0, 0)
                              for g in (needle_b_graph(), needle_a_graph())],
    :gen12_two_edges => () -> [CircularMorphismGraph(g, 0, 0)
                               for g in (circular_two_edge_fixture(),
                                         circular_glued_braids(8, 8, 3))
                               if g !== nothing],
)

# The dot on a big node reads naturally as the expansion, which is the
# direction opposite to the one the driver applies.
const REVERSED = Set([:dot_on_gen12])

SKIP_RULES || for (nm, inputs) in PREFERRED
    r = only(rr for rr in CIRCULAR_RULES if rr.name === nm)
    hits = rule_hits(r, inputs())
    isempty(hits) && continue
    examples[nm] = argmin(h -> length(h[1].graph.nodes), hits)
end

"One row per rule of `registry` that has an example; returns `(name, path)`."
function rule_rows(registry, found, stem)
    rows = Tuple{String,String}[]
    for r in registry
        haskey(found, r.name) || continue
        before, after, fm = found[r.name]
        push!(rows, (string(r.name),
                     step_row(before, after, fm; file = "$stem-$(r.name)",
                              reversed = r.name in REVERSED)))
    end
    return rows
end

"The rows of `registry` whose figure is already on disk."
on_disk(registry, stem) =
    [(string(r.name), "assets/figures/$stem-$(r.name).svg") for r in registry
     if isfile(joinpath(OUT, "$stem-$(r.name).svg"))]

const RULE_ROWS = SKIP_RULES ? on_disk(CIRCULAR_RULES, "rule") :
                  rule_rows(CIRCULAR_RULES, examples, "rule")
const MISSING_RULES = [string(r.name) for r in CIRCULAR_RULES
                       if !(string(r.name) in first.(RULE_ROWS))]

# ---- 7. Double leaves to circular leaves -----------------------------------

println("figures: DL to CL"); flush(stdout)

# `repr` on a diagram gives exactly the picture a notebook shows, whatever the
# type is: the objects carry their own `show(::MIME"image/svg+xml", …)`.
own_tile(x) = repr(MIME("image/svg+xml"), x)

# One double leaf of degree 2 over 121, and the circular leaves it decomposes
# into. This is the computation the whole package exists for.
const DLCL = let en = dl([1, 2, 1], [1, 2, 1])[2][1]
    c = cl_reduce(en.morphism)
    items = Any[(own_tile(en.morphism), "DL($(join(en.e))|$(join(en.f)))  deg $(en.degree)"), :eq]
    for (i, (d, co)) in enumerate(pairs_of(c))
        i > 1 && push!(items, :plus)
        push!(items, (own_tile(d), string(co, " ·  CL")))
    end
    svg_row(items; file = "dl-to-cl")
    (e = en.e, f = en.f, degree = en.degree, terms = length(pairs_of(c)))
end

# The degree-0 circular leaves of 121: the basis the decomposition lands in.
const CLBASIS = let b = cl_basis([1, 2, 1], [1, 2, 1])
    items = Any[]
    for (i, d) in enumerate(b.leaves)
        b.degrees[i] == 0 || continue
        isempty(items) || push!(items, :plus)
        push!(items, (own_tile(d), "CL$i,  degree 0"))
    end
    isempty(items) || svg_row(items; file = "cl-basis")
    (n = length(b.leaves), zero = count(iszero, b.degrees))
end

# ---- 7. The generated pages ------------------------------------------------

println("pages"); flush(stdout)

"Write `docs/src/<name>.md`, through a temporary file."
function write_page(name::AbstractString, text::AbstractString)
    path = joinpath(@__DIR__, "src", name)
    tmp = path * ".tmp"
    write(tmp, text)
    mv(tmp, path; force = true)
    println("  ", name)
end

fig(path, alt) = "![$alt]($path)\n"

# --- rules-gallery.md
io = IOBuffer()
print(io, """
# Every rule, drawn

Each row is one rewrite rule at work: the diagram it matches on the left, the
sum it is replaced by on the right, `0` if the right-hand side is empty. The
captions carry the node count, the number of braid nodes and any region label
that is not `1`. Numbers in the pictures are the debug decorations — arm,
region and leaf numbers.

The examples are not hand-made. Reducing the double leaves of a few word pairs
produces every intermediate diagram, and for each rule the first diagram it
fires on is kept. The rules whose pattern needs a many-armed node or a
particular double connection get a witness built by
`circular/rules/CircularRuleFixtures.jl`; the planar embedding there is
searched, not guessed.

This page is generated by `docs/figures.jl`. The executable version, where the
diagrams can be poked at, is
[`examples/07-all-rules.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/07-all-rules.ipynb).

## The structural rules

`CIRCULAR_RULES`, in registry order — the order in which
[`reduce_circular`](@ref) tries them, which is part of the definition. The
region rules of `CIRCULAR_REGION_RULES` are not drawn separately: they match the
same relations, read off the boundary word of a region instead of off node
shapes, so their pictures would only repeat the ones below.

""")
for (nm, path) in RULE_ROWS
    print(io, "### `", nm, "`\n\n", fig(path, nm), "\n")
end
isempty(MISSING_RULES) ||
    print(io, "\nWithout a drawn example in this build: `",
          join(MISSING_RULES, "`, `"), "`.\n")
write_page("rules-gallery.md", String(take!(io)))

# --- gallery.md
io = IOBuffer()
print(io, """
# A gallery

The notebooks in `examples/` are the executable tutorial. This page is a
still-image tour of the same material, for reading without a Julia kernel.

## A light leaf, two ways

`light_leaf([1,2,1], [1,0,1])` as a `MorphismGraph`, and the same diagram as a
`CircularGraph`: the bottom and top boundaries become one circle, cut at the
marking.

""")
print(io, fig("assets/figures/light-leaf-both.svg", "a light leaf, plain and circular"), "\n")
print(io, """
## The region word

Every region gets a word: the letters crossed on the way in from the marked
outer region. The arrows are the region tree the walk follows, the labels are
`region|word`. Several rules trigger on this word rather than on node shapes.

""")
print(io, fig("assets/figures/region-words.svg", "the region words of a light leaf"), "\n")
print(io, """
## Named diagrams

The fixtures of `diagram/Fixtures.jl`: the unit, the two needles, the barbell.
They are the smallest diagrams on which the dot rules fire.

""")
print(io, fig("assets/figures/fixtures.svg", "unit, needles, barbell"), "\n")
print(io, """
## Double leaves

The double leaves of `121 -> 121` in degree at most 1. `DL(e|f)` names the pair
of subexpressions the leaf is built from.

""")
print(io, fig("assets/figures/double-leaves.svg", "double leaves of 121"), "\n")
TRACE_DL === nothing || print(io, """
## A reduction step

Reduction branches. Most steps replace a diagram by a SUM, so a trace is a tree
and not a line, and the walk below stops at the first branching. Each row is one
step of [`reduce_to_circular_leave`](@ref) on
`DL($(join(TRACE_DL.e))|$(join(TRACE_DL.f)))`, labelled with the rule that
fired. [`circular_trace`](@ref) walks the whole thing interactively.

""")
for (nm, path) in TRACE
    print(io, "**`", nm, "`**\n\n", fig(path, nm), "\n")
end
print(io, """
## Zamolodchikov

The Zamolodchikov axiom: the two sides are the two ways of rewriting the same
reduced expression, and they have to be equal in the category.

""")
print(io, fig("assets/figures/zamolodchikov.svg", "the two sides of the Zamolodchikov axiom"), "\n")
if !isempty(ZTRACE)
    print(io, "The Zamolodchikov stage at work on the right-hand side:\n\n")
    for (nm, path) in ZTRACE
        print(io, "**`", nm, "`**\n\n", fig(path, nm), "\n")
    end
end
print(io, """
## From double leaves to circular leaves

The point of the machine. A double leaf of degree $(DLCL.degree) over
`121 -> 121` on the left, and on the right what [`cl_reduce`](@ref) turns it
into: $(DLCL.terms) circular leaf$(DLCL.terms == 1 ? "" : "s") with a
coefficient in `R`. That is one row of the change-of-basis matrix, drawn.

$(fig("assets/figures/dl-to-cl.svg", "a double leaf and its circular leaves"))
""")
CLBASIS.zero == 0 || print(io, """
The circular leaves of degree 0 over `121 -> 121`, $(CLBASIS.zero) of the
$(CLBASIS.n) leaves of the basis:

$(fig("assets/figures/cl-basis.svg", "the degree-0 circular leaves of 121"))
""")
print(io, """
Every rule of the machine is drawn on its own page,
[Every rule, drawn](rules-gallery.md).
""")
write_page("gallery.md", String(take!(io)))

println("done: ", length(readdir(OUT)), " figures in docs/src/assets/figures")
