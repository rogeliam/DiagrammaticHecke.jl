# The rewrite rules

`src/circular/rules/`. The rules are the machine. Each one is a local pattern
plus a replacement, and the two drivers apply them until nothing matches.

## Naming

Every rule of both registries is drawn, one image equation each, on
[Every rule, drawn](rules-gallery.md).

| prefix | what it names |
|---|---|
| `Rn` | a relation of the diagrammatic category on plain diagrams (R4 dot into braid, R5 braid there and back, R8 trivalent past a crossing, R9 the braid relation, R11 the bigon, R12 commutation) |
| `Cn` | the circular rule implementing a relation on `CircularGraph` |
| `Pn` | the parallel rules |
| `D4` | dot plus connection |
| `Z1`, `Z2` | the Zamolodchikov rules |

## The signature

Every rule has the same shape:

```julia
_fr_<name>(g::CircularGraph) -> Union{Nothing, CircularComboR}
```

`nothing` means no match. An **empty** `CircularComboR` means the diagram is
zero, which is a match, not a failure. [`CIRCULAR_RULES`](@ref) is the registry
of structural rules, [`CIRCULAR_REGION_RULES`](@ref) the ones that match on a
region boundary word, and [`CIRCULAR_STEPS`](@ref) the driver stages.

```julia
[r.name for r in CIRCULAR_RULES]
[r.name for r in CIRCULAR_REGION_RULES]
first.(CIRCULAR_STEPS)
```

## Order is substantive

The registry order is part of the definition, not a cosmetic choice. Several
pairs of rules are complementary at the same site, and the cheaper or more
specific one has to fire first — a same-coloured double connection is `0` and
must be seen before the merge rule reads it as a general double connection, for
instance. The file headers state each ordering constraint where it arises.

## Termination

[`circular_weight`](@ref) is the measure. It is lexicographic and falls
strictly under every structural rule. A pure node count would not do: R4 and
R5 raise the node count, so the first component counts **arms of braid-like
nodes** instead. [`CIRCULAR_WEIGHT_MODE`](@ref) selects between the three
readings [`circular_arm_weight`](@ref), [`circular_sep_weight`](@ref) and
[`circular_chain_weight`](@ref); `:auto` follows the direction switches.

## The two drivers

[`reduce_circular`](@ref) runs the structural rules on a `CircularComboR` and
returns the normal form together with the names of the rules that fired.

```julia
combo, history = reduce_circular(circular(unit_graph()))
history
show_combo(combo)
```

[`reduce_to_circular_leave`](@ref) is the full reduction on decorated
diagrams. It runs in rounds: structural rules, then Zamolodchikov, then
2parallel and D4, then rex fusion, then dot fusion, each success going back to
the head of the round.

```julia
fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly), region_count(fm.graph)))
reduce_to_circular_leave(fdm)
```

[`cl_reduce`](@ref) is the convenience wrapper that takes any morphism.

## Switches

These change results, not output. All default to `true` except the dot slide.

| switch | what it does |
|---|---|
| [`CIRCULAR_ZAMO_ENABLED`](@ref) | the Zamolodchikov replacement stage |
| [`CIRCULAR_DOT_FUSION_ENABLED`](@ref) | the general D4, the dot analogue of rex fusion |
| [`CIRCULAR_REX_FUSION_ENABLED`](@ref) | the rex-cycle relation applied locally |
| [`CIRCULAR_DOT_SLIDE_ENABLED`](@ref) | the hand-built dot-slide surgery (off) |
| [`CIRCULAR_DOT_POLICY`](@ref) | the condition under which a dot moves |
| [`CIRCULAR_EXPLOSION_ENABLED`](@ref) | keep the deepest diagram of a blow-up |

The dot policy is worth spelling out: a dot on a **minimal** node of its colour
pair (6 arms for a 12- or 23-node, 4 for a 13-node), or on any non-braid node,
always moves. On a bigger braid node it moves only if two dotted arms are
adjacent, or if fewer than the Coxeter number of undotted arms of one colour
remain ([`circular_dots_reducible`](@ref)).

A `false` default marks a rule whose place is not settled, not one that is
wrong.

## Watching a reduction

[`circular_trace`](@ref) draws the whole reduction as an image equation,
naming the step that fired at each line; [`circular_step_once`](@ref) does one
step; [`circular_steps`](@ref) selects which stages take part.

```julia
circular_trace(fm.graph, fm)
circular_trace(fz.graph, fz; steps = circular_steps(:zamo))
```

See [`examples/07-all-rules.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/07-all-rules.ipynb),
which shows every rule in the registry with one drawn example.
