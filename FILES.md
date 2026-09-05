# Files in dependency order

The include order in `src/DiagrammaticHecke.jl` is the dependency order: a
file may use everything listed above it and nothing below. Cross-layer jumps
are marked with "needs". Line counts are approximate; the house limit is 700.

## 1. Words and algebra — signed off

The user has reviewed and approved the code of this block. No further
changes are planned here; anything that would touch these files has to be
raised first.

| file | lines | what it is |
|---|---|---|
| `algebra/Coxeter.jl` | 22 | the Coxeter group A₃ from CoxeterGroups.jl: length, descents, short-lex, enumeration |
| `words/CircularWord.jl` | 136 | `CircularWord`: a word over `{1,2,3}` up to rotation |
| `algebra/Ring.jl` | 188 | `SoergelPoly`: the ring R = ℚ[α₁,α₂,α₃], the A₃ action, Demazure operators |
| `algebra/SoergelMatrix.jl` | 202 | determinants and invertibility of matrices over R |
| `algebra/Pairing.jl` | 140 | the Bott–Samelson pairing computed in the Hecke algebra (needs Coxeter) |
| `words/Rules.jl` | 470 | the cyclic rewrite rules on circular words, expansions, move enumeration |
| `words/ColoredReduction.jl` | 188 | coloured words: "no simplification inside one expanded block" |

## 2. Diagrams (`WordGraph`) — signed off

Reviewed and approved like block 1; raise anything that would change these
files before touching them.

| file | lines | what it is |
|---|---|---|
| `diagram/Graph.jl` | 342 | `WordGraph`, `Node`, `Edge`, ports; the generators dot/trivalent/braid |
| `diagram/Fixtures.jl` | 72 | the named example graphs used by tests and notebooks |
| `diagram/IO.jl` | 237 | the shared line-based text format core (`.wg`/`.fwg` families) |
| `diagram/GraphIO.jl` | 76 | `.wg`: save/load a `WordGraph` (needs IO) |
| `diagram/CanonicalKey.jl` | 233 | the one Weisfeiler–Leman canonical-key core for both graph types |
| `diagram/Combo.jl` | 290 | `DiagramCombo`: formal linear combinations keyed by `canonical_key` |
| `reduce/PathDiagram.jl` | 134 | the layer/vertex model of a reduction path (for drawing and printing) |
| `diagram/Helpers.jl` | 41 | incidence and rewiring helpers on a `WordGraph` |

## 3. Morphisms (bottom → top view) — signed off

Reviewed and approved like blocks 1 and 2.

| file | lines | what it is |
|---|---|---|
| `morphism/MorphismGraph.jl` | 154 | `MorphismGraph`: a `WordGraph` with two cuts |
| `morphism/GraphToPath.jl` | 238 | peel a graph outside-in into a reduction path (needs PathDiagram) |
| `morphism/PathToGraph.jl` | 178 | wire a path back into a graph; `path_from_words` |
| `diagram/MorphismIO.jl` | 96 | `.wgm`: save/load a `MorphismGraph` (needs IO, MorphismGraph) |
| `diagram/Faces.jl` | 501 | face tracer, cells, regions; the generic cores the circular side reuses |
| `diagram/WiringCheck.jl` | 430 | Euler test, wiring-convention check, slot-colour audit (needs Faces) |
| `diagram/Decorated.jl` | 241 | `DecoratedDiagram`: a polynomial per inner cell (needs Faces, Ring) |

## 4. Circular diagrams (`CircularGraph`) — signed off

Reviewed and approved.

| file | lines | what it is |
|---|---|---|
| `circular/CircularGraph.jl` | 449 | `CircularNode`/`CircularGraph`: nodes with open arm count; accessors; the dot policy |
| `circular/CircularGraphKey.jl` | 201 | rotation canonicalisation of a node, `circular_canonical_key`, `==`/`hash` |
| `circular/CircularGraphDisplay.jl` | 142 | `is_wired`, console display, the `4n`-cluster fixture |
| `circular/CircularMerge.jl` | 433 | the merge operation on nodes and edges |
| `circular/CircularFaces.jl` | 54 | cells for `CircularGraph` as thin wrappers over `Faces.jl` |
| `circular/CircularConvert.jl` | 186 | `WordGraph` → `CircularGraph`, slot by slot (needs MorphismGraph) |
| `circular/CircularGraphIO.jl` | 149 | `.fwg`/`.fwgm`: save/load (needs IO) |
| `circular/CircularRegion.jl` | 560 | region adjacency, distances, distance words, braid classes (needs BraidMoves later only through `reduced_words`, see note) |
| `circular/CircularBoundaryWord.jl` | 151 | the boundary word of a region, invariant under merging |
| `circular/CircularCombo.jl` | 143 | `CircularCombo`: R-linear combinations keyed by `circular_canonical_key` |
| `circular/CircularDecorated.jl` | 461 | `CircularDecorated`: a polynomial per region |
| `circular/CircularDecoratedIO.jl` | 283 | `.cd`/`.cdb`: decorated graphs and basis files |
| `diagram/Join.jl` | 418 | `join_graphs`, tensor, compose, flips of a `MorphismGraph` |
| `diagram/JoinShift.jl` | 145 | `shift_right`/`shift_left`: moving a boundary letter between the cuts |
| `diagram/JoinCircularMirror.jl` | 150 | `hflip`/`vflip` of a `CircularMorphismGraph` (needs the circular types) |

Note: `CircularRegion.jl` calls `reduced_words` from `morphism/BraidMoves.jl`,
which is included later. This works because the call sits inside a function
body, but it is the one place where the include order is not the call order.

## 5. Light leaves and bases — signed off

Reviewed and approved.

| file | lines | what it is |
|---|---|---|
| `morphism/BraidMoves.jl` | 260 | braid moves on reduced words, `reduced_words`, `canonical_word` (needs Coxeter) |
| `morphism/LightLeaves.jl` | 546 | light leaves and double leaves as `MorphismGraph`s (needs Join, BraidMoves) |
| `morphism/Zamolodchikov.jl` | 167 | the Zamolodchikov relation as a 14-step braid cycle |
| `morphism/CircularLightLeaves.jl` | 128 | the same leaves converted to circular morphisms |
| `morphism/DLBasis.jl` | 625 | the double-leaves basis tables, ranks, and their circular counterparts |
| `circular/CircularExplosion.jl` | 217 | save a diagram to disk when a computation blows up (needs DecoratedIO) |

## 6. Circular rules and the driver

| file | lines | what it is |
|---|---|---|
| `circular/rules/CircularRules.jl` | 493 | rule signature, shared helpers, C1–C5 and the rule list |
| `circular/rules/CircularMergeRules.jl` | 547 | the merge family: mono double, C10, C12, the merge lock |
| `circular/rules/CircularBraidRules.jl` | 389 | C6–C8: the two-colour rules |
| `circular/rules/CircularBraidChannels.jl` | 261 | channel normalisation at braid pairs |
| `circular/rules/CircularGen12Merge.jl` | 293 | C13: two braid nodes plus two trivalents → one gbraid |
| `circular/rules/CircularDotOnGen12.jl` | 135 | C14: a dot on a gbraid arm |
| `circular/rules/CircularGen12Expand.jl` | 323 | C13⁻¹: the hard-coded preimages of an 8-armed gbraid |
| `circular/rules/CircularViaPreimage.jl` | 118 | the combinator "expand, reduce, compare" used by C18 |
| `circular/rules/CircularBraidOnGen12.jl` | 276 | C15: a braid hanging off a gbraid by three edges; the shared glue-3 surgery |
| `circular/rules/CircularGen12OnGen12.jl` | 91 | C16: two gbraids joined by three edges |
| `circular/rules/CircularGen12Null.jl` | 175 | C22: two braid-like nodes joined by four or more edges are 0 |
| `circular/rules/CircularGen12TwoEdges.jl` | 252 | C18: gbraid and braid joined by two edges or a channel |
| `circular/rules/CircularDotMerge.jl` | 156 | C19/C20: a dot on a braid-like node |
| `circular/rules/CircularWeight.jl` | 113 | the termination weight of the driver |
| `circular/rules/CircularDriver.jl` | 479 | `reduce_circular`: the rule driver on `CircularComboR` |
| `circular/rules/CircularRegionRules.jl` | 192 | rules matching on region boundary words |
| `circular/rules/ZamoRules.jl` | 385 | Zamolodchikov fixtures and anchors on `WordGraph` |
| `circular/rules/CircularDecoratedRules.jl` | 700 | D4 (dot plus connection) on decorated circular graphs |
| `circular/rules/CircularLeaveDriver.jl` | 691 | `reduce_to_circular_leave`: the circular-leaf reduction driver |
| `circular/CircularComponents.jl` | 238 | floating connected components and their local evaluation |
| `circular/rules/CircularD4Node.jl` | 377 | the node variant of D4; scalar extraction; fusion step |
| `morphism/DLBasisWriter.jl` | 449 | writing a morphism in the DL basis; structure constants (needs LeaveDriver) |
| `morphism/CircularPairing.jl` | 142 | the pairing ⟨DLᵢ, DLⱼ⟩ as a function |
| `circular/rules/ZamoTermRules.jl` | 293 | the two Zamolodchikov rules as whole-term replacements |
| `circular/rules/CircularZamoRegion.jl` | 419 | the Zamolodchikov rule as a multi-region pattern: pattern and matcher |
| `circular/rules/CircularZamoRegionStep.jl` | 194 | the step: delete the cluster, splice in the right-hand side |
| `circular/rules/CircularZamoRegionInverse.jl` | 299 | Zamo triangles and regions, the backward direction and its guard |
| `circular/rules/CircularParallelMerge.jl` | 255 | P1: merge parallel 1/3 edges |
| `circular/rules/Circular2Parallel.jl` | 663 | 2parallel: the conditions and `find_circular_2parallel` |
| `circular/rules/Circular2ParallelApply.jl` | 84 | the 2parallel surgery and the rule step |
| `circular/rules/CircularRegionWord.jl` | 694 | the region-word criterion that triggers 2parallel |
| `circular/rules/CircularDoubleBraid.jl` | 293 | P3: two braids doubly connected (disabled) |
| `circular/rules/CircularLeafDegree.jl` | 161 | the degree invariant and the per-degree count test |
| `circular/rules/CircularInverseRules.jl` | 442 | the circular rules run backwards |
| `circular/rules/CircularRuleFixtures.jl` | 193 | named example diagrams for the rules no short double leaf produces |
| `circular/rules/CircularDotSlide.jl` | 318 | the hand-built dot-slide surgery behind the dihedral region-word step |
| `morphism/CLBasis.jl` | 219 | the circular-leaf basis of a word pair, computed and cached |
| `circular/rules/CircularRexRelation.jl` | 250 | the rex-cycle relation, derived and cached (needs CLBasis) |
| `circular/rules/CircularSplice.jl` | 122 | splice a morphism into a cut curve |
| `circular/rules/CircularRexFusion.jl` | 232 | the rex-cycle relation applied locally (`CIRCULAR_REX_FUSION_ENABLED`, off) |
| `circular/rules/CircularDotFusion.jl` | 137 | the general D4 rule, dot analogue of rex fusion |

## 7. Rendering

| file | lines | what it is |
|---|---|---|
| `render/Morphism.jl` | 493 | radial and linear SVG of a word reduction |
| `render/Graph.jl` | 294 | SVG/PNG output of a `WordGraph`; the notebook display sink |
| `render/Rect.jl` | 286 | a `MorphismGraph` drawn as a rectangle |
| `render/tutte/Avoid.jl` | 193 | geometry helpers shared by both Tutte renderers |
| `render/tutte/Markers.jl` | 250 | boundary markers and gap angles |
| `render/tutte/Tutte.jl` | 166 | node neighbours and the braid-node leg angles |
| `render/tutte/TutteSVG.jl` | 534 | `tutte_svg`: the Tutte embedding of a `WordGraph` as SVG |
| `render/tutte/Wiring.jl` | 163 | debug labels and the wiring table |
| `render/tutte/Display.jl` | 157 | the `display_*` helpers |
| `render/CircularTutte.jl` | 492 | the Tutte-style layout for `CircularGraph` |
| `render/CircularTutteSVG.jl` | 661 | the SVG emitter for that layout |
| `render/Cells.jl` | 552 | cell anchors: centroid placement and the repair passes |
| `render/CellLabels.jl` | 132 | label text and SVG text elements for cells and regions |
| `render/CellDistance.jl` | 129 | cell numbers, cell adjacency and distances, their SVG |
| `render/CircularCells.jl` | 268 | cell labels for circular graphs (shares the anchor core) |
| `render/CircularRegionTree.jl` | 120 | the dual-graph overlay: the region-word tree |
| `render/Decorated.jl` | 90 | polynomial labels on the Tutte SVG |
| `render/CircularSteps.jl` | 331 | the reduction drawn as an image equation; collects all step functions |

## 8. Word exploration (independent of diagrams)

| file | lines | what it is |
|---|---|---|
| `words/SequenceGraph.jl` | 642 | the reduction graph of circular words, reachability of ε |
| `words/Checkpoint.jl` | 272 | save/resume a long build; the `.cwg` format |

## Tests and notebooks

`test/runtests.jl` runs the groups `words, algebra, diagram, faces, morphism,
circular, decorated, pairing, dlbasis, circularpairing, clbasis,
circularcomponents, zamo` in that order; `test/fixtures_circular.jl` is loaded
before the circular files, which follow the include order of `src/circular/`.
The notebooks `examples/01`–`05` follow layers 1 to 6, `03` and `05` closing
with the file formats of their layer; `06` is the region word, which comes
BEFORE `07`, the rule inventory, because several rules trigger on it; `08` is
the double-leaves basis together with the coordinates of an arbitrary morphism
and the change of basis, `09` the Zamolodchikov rules, and `10` the reduction
of double leaves to circular leaves.

## The manual

`docs/make.jl` builds the manual with Documenter; `docs/figures.jl` builds its
pictures. The figure script renders every diagram of the manual to
`docs/src/assets/figures/` and writes two generated pages, `docs/src/gallery.md`
(the notebook material as still images) and `docs/src/rules-gallery.md` (one
image equation per rewrite rule, harvested the same way `examples/07` harvests
them). Both are generated, so they are edited by editing `docs/figures.jl`.
