# DiagrammaticHecke.jl

Computer algebra for the **diagrammatic Hecke category** in type A₃: circular
words over the alphabet `{1,2,3}`, planar Soergel diagrams (`WordGraph` /
`CircularGraph`), local rewrite rules, and the light-leaves / double-leaves
bases — together with SVG rendering of every diagram.

## Installation

The package is not registered. From a Julia session:

```julia
using Pkg
Pkg.add(url = "https://github.com/<user>/DiagrammaticHecke.jl")
```

or, from a clone,

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Julia 1.11 or newer. The dependencies are `AbstractAlgebra.jl` for the
coefficient ring `R = ℚ[α₁,α₂,α₃]` and the exact linear algebra over it,
[`CoxeterGroups.jl`](https://github.com/ulthiel/CoxeterGroups.jl) for the
Coxeter group of type A₃ (installed from GitHub through the `[sources]` entry
of `Project.toml`), and `Cairo.jl` + `Rsvg.jl` to rasterize diagram SVGs for
front-ends that do not render `image/svg+xml`. Everything else — words,
diagrams, rewrite rules — is plain Julia.

## Example

```julia
using DiagrammaticHecke

m  = light_leaf([1, 2, 1], [1, 0, 1])   # a light leaf BS(121) -> BS(ε)
fm = circular_morphism(m)               # the same diagram, circular model
display(fm)                             # draws it

labels = fill(one(SoergelPoly), region_count(fm.graph))
combo  = reduce_to_circular_leave(CircularDecoratedMorphism(fm, labels))
show_combo(combo; m = fm)               # the normal form: Σ cᵢ·Dᵢ over R
```

## Documentation

The manual lives in `docs/`: a guide that follows the layers of the package and
an API reference generated from the docstrings. To build it locally,

```
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=. docs/figures.jl
julia --project=docs docs/make.jl
```

and open `docs/build/index.html`. Set `DOCS_REPO` to the repository once it has
a remote; the placeholder in `docs/make.jl` is what stands in until then.

`docs/figures.jl` is the middle step. It draws every picture of the manual with
the package itself, writes them to `docs/src/assets/figures/` and generates the
two picture pages: a gallery of the notebook material, and one image equation
per rewrite rule. Run it again whenever a rule, a fixture or a renderer
changes.

## Notebooks

Ten notebooks in `examples/`, meant to be read top to bottom:

1. `01-circular-words-and-rules.ipynb` — circular words and the cyclic rewrite rules
2. `02-the-soergel-ring.ipynb` — `R = ℚ[α₁,α₂,α₃]`, the A₃ action, the pairing
3. `03-word-diagrams.ipynb` — `WordGraph`, faces, regions, saving and loading
4. `04-light-leaves.ipynb` — subexpressions, light and double leaves, Zamolodchikov
5. `05-circular-diagrams-and-rules.ipynb` — `CircularGraph`, the rule driver, the circular file formats, the debug options one by one
6. `06-the-region-word.ipynb` — how the word is computed, the marking it depends on, the two variants and their two distances, and the transition the rules trigger on
7. `07-all-rules.ipynb` — every rewrite rule with one drawn example
8. `08-the-double-leaves-basis.ipynb` — the DL basis by degree, the coordinates of an arbitrary morphism, the change of basis, composites and flips
9. `09-zamolodchikov.ipynb` — the Zamolodchikov axiom and the derivation of the second Zamolodchikov rule
10. `10-double-leaves-to-circular-leaves.ipynb` — the decomposition `DL ↦ CL`, and what the coefficients look like

The region word comes before the rules on purpose: several of them trigger on it.

They need this project's environment. Install a kernel bound to it once:

```julia
julia --project=. -e 'using IJulia;
                      installkernel("DiagrammaticHecke", "--project=" * pwd())'
```

then open `examples/` in Jupyter and pick the **DiagrammaticHecke** kernel. To
run notebook-style code in a plain script instead, call
`notebook_display_sink!()` first: it registers a display that accepts
`image/svg+xml` and friends, which a bare `julia script.jl` lacks.

## Layout

- `src/` — the package, grouped by layer (`words/`, `algebra/`, `diagram/`,
  `morphism/`, `circular/`, `render/`). The include order in
  `src/DiagrammaticHecke.jl` is the dependency order.
- `test/` — one file per layer.
- `examples/` — the notebooks above.

See [PLAN.md](PLAN.md) for the layer structure and the house rules.

## Tests

```
julia --project=. -e "using Pkg; Pkg.test()"                 # everything
DH_TESTS=words,circular julia --project=. test/runtests.jl   # selected groups
julia --project=. test/runtests.jl --list                    # group names
```

## TODO

The four things this package is missing, in the order they are being worked on.

1. **No 13-node.** The author recently found that merging a 1 and a 3 into one
   node is not correct. There is no 13-node mathematically, and building one
   makes the diagram non-planar. Most of what this package computes is
   unaffected and correct as it stands. The rework of the merge family, and of
   everything downstream that reads such a node, comes in a future version.
2. **Every Coxeter group.** The package is type A₃ throughout: the alphabet, the
   ring, the group, the braid lengths, the fixtures. The goal is to take any
   Coxeter group and build the package for it, with the Coxeter matrix as the
   input the rules read their braid orders off.
3. **The Zamolodchikov relation in H₃.** In A₃ it is an axiom and the rules Z1
   and Z2 encode it. In H₃ it has to be computed first.
4. **Idempotents and Soergel bimodules.** The bases are here, the idempotents
   that cut out the indecomposables are not, and neither is the bridge back to
   bimodules.

## Authorship

The mathematics is mine, and it will be presented in an upcoming paper,
*Circular leaves for the diagrammatic Hecke category*
([`Ro-circular-leaves`](CITATION.bib), in preparation).

The code is largely written by large language models: Claude for most of it,
Kimi for a smaller part, and GitHub Copilot alongside. They did the
implementation work under my direction, from my ideas and my specifications. I
tested it a lot. The responsibility for what the package claims is mine.

## Citing

If this package is useful in your work, please cite it, and the paper it belongs
to once that is out. Both entries are in [CITATION.bib](CITATION.bib):

```bibtex
@software{DiagrammaticHecke,
  author  = {Rogel, Liam},
  title   = {{DiagrammaticHecke.jl}: a {J}ulia package for the diagrammatic
             {H}ecke category},
  year    = {2026},
  note    = {In preparation},
}

@unpublished{Ro-circular-leaves,
  author  = {Rogel, Liam},
  title   = {Circular leaves for the diagrammatic {H}ecke category},
  year    = {2026},
  note    = {In preparation},
}
```

## License

GNU General Public License, version 3 or later. The full text is in
[LICENSE](LICENSE).

    DiagrammaticHecke.jl, computer algebra for the diagrammatic Hecke category.
    Copyright (C) 2026 Liam Rogel

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version. It is distributed WITHOUT ANY WARRANTY; see the licence
    for details.


