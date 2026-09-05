# PLAN — structure and house rules

## Layers

The include order in `src/DiagrammaticHecke.jl` **is** the dependency order:
each file may use everything included above it, and a shared core has to come
before every file that uses it (`diagram/IO.jl` before `GraphIO.jl`, for
instance). That is the one structural rule to keep.

```
src/DiagrammaticHecke.jl   include order = dependency order
src/words/       CircularWord, Rules (cyclic rewriting), ColoredReduction,
                 SequenceGraph (backward BFS from ε), Checkpoint
src/algebra/     Coxeter (the group of type A₃ from CoxeterGroups.jl), Ring (R = ℚ[α₁,α₂,α₃],
                 Demazure), SoergelMatrix, Pairing
src/diagram/     WordGraph core: Graph, IO, CanonicalKey, Combo, Helpers,
                 Faces, WiringCheck, Decorated, Join
src/morphism/    MorphismGraph, GraphToPath/PathToGraph, BraidMoves,
                 LightLeaves, Zamolodchikov, and the bases
                 (CircularLightLeaves, DLBasis, CircularPairing, CLBasis)
src/circular/    CircularGraph and everything on it (merge, faces, regions,
                 decorated, components) + rules/ (the driver and rule files)
src/render/      Tutte, Rect, Cells, the circular renderers, steps
test/            one file per layer, group selection via DH_TESTS
examples/        one explanatory notebook per layer
```

Two placements worth stating:

* the bases live under `morphism/`, not in a `basis/` folder of their own;
* the named example graphs live in `diagram/Fixtures.jl`, not in the renderer.

## Two canonical keys, two combo types

`Combo`/`CircularCombo` and `Decorated`/`CircularDecorated` are separate types
because the two canonical keys are different notions of "the same diagram".

* `circular_canonical_key` — the START of the circular word is part of the
  identity. It deliberately does not minimise over
  `_circular_leaf_rotations(g.word)`, because the left marking belongs to the
  diagram. Rotating the boundary word gives a DIFFERENT key.
* `canonical_key` — no such thing; orientation does not matter.

Some rules can be made orientation-free and some cannot, so the distinction has
to stay visible in the types. Two combo types keep it visible; one parametric
type would blur it into a dispatch detail.

The `Any` representative in `CircularCombo` is bimodal, not heterogeneous:
`_lift` (`circular/rules/CircularDriver.jl`) is the boundary, and it CONVERTS by
building a fresh combo rather than mixing types within one.

## Rule switches

`CIRCULAR_REX_FUSION_ENABLED` defaults to `true`: the trigger is the region word,
and the driver calls the step in the `:after_d4` group — structural rules, then
Zamo, then 2parallel and D4, then rex fusion, then dot fusion, each going back to
the head of the round.

`CIRCULAR_ZAMO_ENABLED`, `CIRCULAR_DOT_FUSION_ENABLED`,
`CIRCULAR_EXPLOSION_ENABLED` and the dot policy `CIRCULAR_DOT_POLICY` default
to `true`. The dot policy is the condition under which a dot moves: a dot on a
minimal node of its colour pair (6 arms for a 12- or 23-node, 4 for a 13-node)
or on any non-braid node always; on a bigger braid node only if two dotted arms
are adjacent, or if fewer than the Coxeter number of undotted arms of one colour
remain (`circular_dots_reducible`).

## Names of relations and rules

`Rn` names a relation of the diagrammatic category on plain diagrams (R4 dot
into braid, R5 braid there and back, R8 trivalent past a crossing, R9 the braid
relation, R11 the bigon, R12 commutation). `Cn` names the circular rule that
implements a relation on `CircularGraph`; `Pn` the parallel rules, `D4` the
dot-plus-connection rule, `Z1`/`Z2` the Zamolodchikov rules.

## House rules

- Comments and docs **English only**.
- Docstrings say what a thing IS and WHY it is built that way, in the present
  tense. No chronicle: no "used to", no "renamed from", no plan or option
  references, no dates on a reason, no pointers to files outside this
  repository.
- No file over ~700 lines; split by topic, not by history.
- Tests stay minimal and focused (a handful of assertions per feature); test
  data fixtures go in `test/fixtures_*.jl`.
- Every layer gets ONE explanatory notebook in `examples/` that shows the
  concepts with drawn diagrams (debug mode on, `display` only).
- The manual is a picture book too. `docs/figures.jl` renders every figure with
  the package itself and generates `docs/src/gallery.md` and
  `docs/src/rules-gallery.md`; a new rule shows up there by itself, and a rule
  without an example is reported by name at the end of the page.
- Verify a "green" suite by the **test count**, not the exit status.
