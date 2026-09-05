# Words and rules

`src/words/`. This layer is pure combinatorics: no diagrams, no polynomials.
It is the boundary of everything drawn later.

## Circular words

A [`CircularWord`](@ref) is a word over `{1,2,3}` up to **cyclic rotation**.
Two words that differ by a rotation are the same element.

```julia
w = CircularWord([1, 3, 1, 2])
CircularWord([3, 1, 2, 1]) == w    # true, a rotation
letters(w), compact(w), length(w)
rotations(w)
```

[`EMPTY`](@ref) is the empty word. The reason for the quotient is geometric:
the word is read off the boundary circle of a diagram, and the circle has no
starting point of its own.

## The rewrite rules

[`BASE_RULES`](@ref) lists the rules, each usable in both directions and
applied **cyclically** — a pattern may wrap around the end of the word.

```julia
[m for m in moves(w) if m.rule == :braid_12]   # every site, wrapping included
neighbors(w)                                   # distinct results, with labels
```

In `1312` the braid pattern `121` sits at the wrap. A linear matcher misses it;
that is the whole point of the type.

[`moves`](@ref) enumerates sites, [`neighbors`](@ref) deduplicates them to the
distinct resulting words paired with the length change.

## Directed moves

[`EXTRA_REDUCERS`](@ref) are strictly reducing shortcuts derivable from the
base rules. [`down_neighbors`](@ref) and [`up_neighbors`](@ref) split the moves
by whether they shorten the word.

```julia
down_neighbors(CircularWord([1, 2, 1, 2, 1, 2]))   # reaches ε directly
up_neighbors(CircularWord([1]))
```

## Braid expanders

Some words shorten only after they first grow. The eight **braid expanders**
[`BRAID_EXPANDERS`](@ref) tile a word into blocks that each expand.

```julia
v = CircularWord([1, 2, 2, 1])
has_full_expansion(v)
ex = a_full_expansion(v)      # (rotation, blocks, result)
increasing_expansions(v)
```

[`ColoredWord`](@ref) records which block a letter came from, so that "no
simplification inside one expanded block" can be enforced;
[`colored_down_path`](@ref) and [`forget_color`](@ref) run and drop that
bookkeeping.

## The reduction graph

[`build_reduction_graph`](@ref) builds the graph of all circular words up to a
length bound, with one edge per move. It is the search space for the question
"which words reduce to ε".

```julia
g  = build_reduction_graph(maxlen = 10)
dg = non_increasing_subgraph(g)
reaches_empty_down(down_graph(collect(nodes(g))))
counts_by_length(reducible_words(g))
```

The build is long enough to want a checkpoint. [`build_or_resume`](@ref)
carries a [`BuildState`](@ref) across sessions through the `.cwg` format
([`save_build`](@ref), [`load_build`](@ref)), and [`cached_graph`](@ref) keeps
one graph per bound on disk.

## What a word draws

Every rule above is a relation between planar diagrams, and every circular word
is the boundary of one. The rest of the manual makes that literal.

```julia
b = braid(1, 2; m = 3)
compact(b.word)
display(b)
```

See [`examples/01-circular-words-and-rules.ipynb`](https://github.com/<user>/DiagrammaticHecke.jl/blob/main/examples/01-circular-words-and-rules.ipynb).
