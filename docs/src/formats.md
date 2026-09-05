# File formats

Every object that takes time to build has a text format. All of them are
line-based, versioned and git-diffable: a diagram in a repository should be
readable and reviewable, not a blob.

## The graph family

There are two graph record formats, and they are the same format with one line
swapped.

| suffix | type | writer / reader |
|---|---|---|
| `.wg` | [`WordGraph`](@ref) | [`wordgraph_to_string`](@ref), [`wordgraph_from_string`](@ref) |
| `.wgm` | [`MorphismGraph`](@ref) | [`morphismgraph_to_string`](@ref), [`morphismgraph_from_string`](@ref) |
| `.fwg` | [`CircularGraph`](@ref) | [`circulargraph_to_string`](@ref), [`circulargraph_from_string`](@ref) |
| `.fwgm` | [`CircularMorphismGraph`](@ref) | [`circularmorphismgraph_to_string`](@ref), [`circularmorphismgraph_from_string`](@ref) |

The layout, with `<tag>` being `wg` or `fwg`:

```
#wg 2                     format version
word 1212121              boundary word, "eps" for ε
# nodes                   one per line, in index order
n braid 1 2 3
# edges                   one per line: colour  portA  portB
e 1 L1 n2.3               colour 1, leaf 1 to node 2 slot 3
cuts 0 3                  optional, the morphism cuts
```

A port is `L<k>` for boundary leaf `k`, `n<node>.<slot>` for a node connector,
or `C<colour>` for a free circle. Only the node line differs between the two
families, because a [`Node`](@ref) carries colours and `m` while a
[`CircularNode`](@ref) carries an arbitrary arm vector. Several records go in
one file separated by `---`, each with a free-form metadata comment;
[`save_morphisms`](@ref), [`load_morphisms`](@ref),
[`save_circular_morphisms`](@ref) and [`load_circular_morphisms`](@ref) handle
those containers.

## Decorated diagrams and bases

| suffix | content |
|---|---|
| `.cd` | a [`CircularDecorated`](@ref): a label header wrapping an `.fwg` body |
| `.cdb` | a basis of decorated diagrams |
| `.cws` | a word set, from the reduction-graph build |
| `.cwg` | a [`BuildState`](@ref): a long build, resumable |

```julia
save_circulardecorated(path, d);  load_circulardecorated(path)
save_cl_basis(path, b);           load_cl_basis(path; zamo = nothing)
save_wordset(path, ws; maxlen);   load_wordset(path)
save_build(path, st);             load_build(path)
```

Polynomials go through [`soergelpoly_to_string`](@ref) and
[`soergelpoly_from_string`](@ref), so a label is readable in the file.

## Caches

Two computations are expensive enough to cache on disk: the circular-leaf basis
of a word pair ([`cl_cache_path`](@ref)) and the rex-cycle relation
([`circular_rex_cache_path`](@ref)). Both live under `data/cache/`.
[`clear_dl_cache!`](@ref) and [`clear_graph_cache!`](@ref) drop the in-memory
ones.

## Blow-ups

When a reduction grows without terminating, the deepest diagram is written out
by [`circular_save_explosion`](@ref) into `CIRCULAR_EXPLOSION_DIR[]`, together
with the weight, the package version and the state of every switch — enough to
reproduce the case. [`circular_explosions`](@ref) lists what has been captured,
[`circular_reset_explosion`](@ref) clears it.
