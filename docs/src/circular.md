# Circular diagrams

`src/circular/`. This is the computing core. A [`CircularGraph`](@ref) is the
same picture as a [`WordGraph`](@ref), with one difference that decides
everything: a [`CircularNode`](@ref) carries an **open arm count**, so two
nodes can merge into one bigger node instead of producing a new diagram shape
for every case.

```julia
m  = light_leaf([1, 2, 1], [1, 0, 1])
fm = circular_morphism(m)          # keeps the boundary marking
display(fm)
arm_count(fm.graph.nodes[1]), arms(fm.graph.nodes[1])
circular_degree(fm.graph)
```

[`circular(g)`](@ref circular) converts a `WordGraph`; `circular_morphism` does
the same for a morphism and keeps its two cuts.

The same light leaf as a `MorphismGraph` and as a `CircularGraph`:

![a light leaf, plain and circular](assets/figures/light-leaf-both.svg)

## Merging

[`merge_at_edge`](@ref), [`merge_nodes`](@ref) and [`merge_all`](@ref) are the
single operation the rules are written in terms of. Two braid-like nodes joined
by an edge become one node whose arm list is the two arm lists spliced at the
shared edge. A `2k`-armed node stands for the generalised braid generator on
`k` strands, so the rules never enumerate cases by node type.

## Keys and combinations

[`circular_canonical_key`](@ref) canonicalises a node by rotation and the whole
graph up to relabelling — but **not** up to rotating the boundary word. The
start of the circular word is part of the identity of a circular diagram,
because the left marking of a morphism is. That is the one place where the two
canonical keys of the package differ, and the reason `CircularCombo` and
`DiagramCombo` are separate types.

[`CircularCombo`](@ref) and `CircularComboR` are the `ℤ`- and `R`-linear
combinations keyed that way. [`pairs_of`](@ref) iterates the terms;
[`show_combo`](@ref) prints or draws the sum.

## Regions and their words

A region is a connected piece of the complement of the diagram. Regions carry
the polynomial labels, and their combinatorics drives several rules.

```julia
regions(g)
circular_region_distances(m)         # node-weighted distance from the marking
circular_region_tree_words(m)        # the colours read along the tree path
circular_region_words(m)             # the other variant, a different alphabet
circular_region_boundary_word(m, r)  # invariant under merging
```

There are two distances and two word variants, and they answer different
questions. The tree word says *how one gets to a region from the marking*; it
therefore depends on the marking, and moving the cut changes it. That is not a
defect: the marking is part of the question.

[`circular_unreduced_transition`](@ref) is the criterion the fusion rules
trigger on. It does not look for an unreduced region word. It looks for a
**transition**: a region whose word plus the next colour fails to be reduced,
together with the pair of edges that causes it. The result names the region,
the word, the colour, the `kind` (`:edge` or `:dot`) and the `case`
(`:direct` or `:dihedral`).

## Decorated circular diagrams

[`CircularDecorated`](@ref) is a `CircularGraph` with one `SoergelPoly` per
region, and [`CircularDecoratedMorphism`](@ref) adds the cuts. This is the type
the reduction actually runs on.

```julia
fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
combo = reduce_to_circular_leave(fdm)
show_combo(combo; m = fm)
```

## Components and blow-ups

[`circular_connected_components`](@ref) and
[`circular_extract_floating_components`](@ref) pull a floating component out of
a diagram and evaluate it locally; [`circular_extract_scalars`](@ref) turns a
closed component into its coefficient.

When a reduction grows instead of shrinking, `CIRCULAR_EXPLOSION_ENABLED[]`
keeps the deepest diagram so [`circular_save_explosion`](@ref) can write it to
disk for inspection.

See [`examples/05-circular-diagrams-and-rules.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/05-circular-diagrams-and-rules.ipynb)
and [`examples/06-the-region-word.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/06-the-region-word.ipynb).
