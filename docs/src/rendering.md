# Rendering

`src/render/`. Every diagram type draws itself. In a notebook `display(g)` is
enough; in a script, register the display sink first.

```julia
notebook_display_sink!()    # a display that accepts image/svg+xml
```

`docs/figures.jl` is the third context: it renders the figures of this manual
to plain SVG files with the same functions, so the pictures here are the ones
the package actually draws.

## Layouts

There are three layouts, for three different questions.

| function | layout | for |
|---|---|---|
| [`tutte_svg`](@ref) | Tutte embedding, boundary on a circle | a `WordGraph`, the default |
| [`circular_tutte_svg`](@ref) | the same for a `CircularGraph` | many-armed nodes |
| [`morphism_rect_svg`](@ref) | a rectangle, bottom to top | a `MorphismGraph` read as a morphism |
| [`radial_svg`](@ref) / [`linear_svg`](@ref) | a word reduction, outside-in or in rows | a reduction path |

Each has a `display_*` companion ([`display_tutte`](@ref),
[`display_circular_tutte`](@ref), [`display_rect`](@ref)) and a file writer
([`save_diagram`](@ref), [`save_radial`](@ref), [`save_linear`](@ref)).
[`diagram_png`](@ref) and [`morphism_rect_png`](@ref) rasterize through Cairo
for front-ends that do not render SVG.

## Labels on regions

[`decorated_svg`](@ref) and [`circular_decorated_svg`](@ref) draw the
polynomial of each region into the region, placed at a repaired centroid.
[`display_decorated`](@ref) and [`display_circular_decorated`](@ref) are the
notebook forms.

## Debug overlays

```julia
debug_labels(true)     # arm and leaf numbers in every WordGraph picture
debug_circular(true)   # arm, region and leaf numbers in every circular picture
display_cell_number(g) # cell ids
display_cell_distance(m)
```

Every notebook opens with one of the two switches on, because a picture without
numbers cannot be checked against a printed region index.

## Reductions as pictures

[`circular_trace`](@ref) draws a whole reduction as a stack of image equations
`D → c₁·D₁ + …`, one row per step, with the name of the step that fired.
[`show_circular_step_row`](@ref) draws a single before/after row, which is how
the rule inventory notebook shows one example per rule.

## Sums

[`show_combo`](@ref) prints or draws a `CircularComboR`; passing `m = fm` keeps
the marking of the original morphism so the terms are drawn the same way up.
