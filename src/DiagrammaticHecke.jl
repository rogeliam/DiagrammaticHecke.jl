module DiagrammaticHecke

# Cairo + Rsvg rasterize the diagram SVG to PNG for front-ends that do not
# render image/svg+xml (notably VS Code notebooks). See render/Graph.jl.
using Cairo, Rsvg
using Base64
using PrecompileTools

# The diagrammatic Hecke category in type A3: circular words, planar Soergel
# diagrams and their local rewrite rules, light/double leaves, rendering.
#
# The include order below IS the dependency order: each file may use everything
# included above it. Layers:
#
#   words/     circular words over {1,2,3} and the cyclic rewrite rules
#   algebra/   the Coxeter group of type A₃, the coefficient ring
#              R = QQ[α₁,α₂,α₃], linear algebra over it, the Soergel pairing
#   diagram/   the WordGraph diagram type
#   morphism/  the bottom→top view, light/double leaves and the bases
#   circular/  CircularGraph — the computing core + rule driver
#   render/    SVG rendering

include("algebra/Coxeter.jl")

include("words/CircularWord.jl")
export CircularWord, EMPTY, letters, rotations, compact

include("algebra/Ring.jl")
export SoergelPoly, R, alpha, act, demazure, δ, degree
export soergel_monomial, poly_terms, poly_from_terms, constant_term, soergel_str

include("algebra/SoergelMatrix.jl")
export soergel_det, soergel_det_peeled, is_invertible_over_frac_ring
export is_permutation_triangular

include("algebra/Pairing.jl")
export bs_pairing, pairing_table, expected_count

include("words/Rules.jl")
export Rule, Move, BASE_RULES, EXTRA_REDUCERS, EXTENDED_RULES, UP_RULES
export moves, neighbors, down_neighbors, up_neighbors, sign_of
export increasing_expansions, has_full_expansion, a_full_expansion, all_full_expansions
export uses_all_letters
export BRAID_EXPANDERS, BRAID_EXPANDERS_2, BRAID_EXPANDERS_1

include("words/ColoredReduction.jl")
export ColoredWord, forget_color, colored_from_expansion
export down_neighbors_colored, colored_down_path

include("diagram/Graph.jl")
export WordGraph, Node, Edge, Port, Leaf, NodePort, Circle, boundary
export dot, trivalent, braid, join_trivalents, show_wordgraph
export opposite_slot, is_planar_braid
const NP = NodePort
export NP

include("diagram/Fixtures.jl")
export unit_graph, needle_a_graph, needle_b_graph, r3_braid_dot, r4_braid_dot, barbell_graph

# diagram/IO.jl — the line-based text-format core shared by the `.wg`/`.wgm`
# (GraphIO/MorphismIO, right below) and `.fwg`/`.fwgm` (circular/CircularGraphIO.jl,
# further down) formats. Each of those supplies only its node line.
include("diagram/IO.jl")
export GraphTextFormat, graph_to_text, graph_from_text

include("diagram/GraphIO.jl")
export wordgraph_to_string, wordgraph_from_string, save_wordgraph, load_wordgraph

# diagram/CanonicalKey.jl — the shared Weisfeiler–Leman canonical-key core, used
# by both `canonical_key` (Combo.jl, right below) and `circular_canonical_key`
# (circular/CircularGraph.jl, further down).
include("diagram/CanonicalKey.jl")

include("diagram/Combo.jl")
export DiagramCombo, DiagramRule, canonical_key, coefficient, pairs_of, simplify
export slot_rotation_keys

include("reduce/PathDiagram.jl")
export Vertex, VKind, VDot, VTrivalent, VBraid, VComm
export PathDiagram, path_diagram
export reduction_path_string, show_reduction_path

# diagram/Helpers.jl — incidence / rewiring helpers on a WordGraph.
include("diagram/Helpers.jl")

include("morphism/MorphismGraph.jl")
include("morphism/GraphToPath.jl")
include("morphism/PathToGraph.jl")
export MorphismGraph, morphism_graph, bottom, top, leaf_colour, show_morphism
export graph_to_path, path_to_graph
export path_from_words
export dot_count, trivalent_count

include("diagram/MorphismIO.jl")
export morphismgraph_to_string, morphismgraph_from_string, save_morphism, load_morphism
export MorphismEntry, save_morphisms, load_morphisms

include("diagram/Faces.jl")
export Cells, Cell, face_count, inner_faces, face_of_port, gap_cell, sector_cell
export show_faces
export region_count, boundary_regions, regions, Region

# diagram/WiringCheck.jl — the package's sanity check: the Euler test and the
# direct wiring-convention test. Needs the tracer/Cells from Faces.jl above.
include("diagram/WiringCheck.jl")
export euler, is_planar_embedding
export WiringViolation, check_wiring, check_wiring_table, show_wiring_check
export slot_colour_audit, slot_colour_audit_log

include("diagram/Decorated.jl")
export DecoratedDiagram, DiagramComboR, decorated, face_label, outer_label, poly_type

# circular/ — the diagram type with open arm count: a CircularNode carries only
# its cyclic arm-colour sequence. Needs Faces (tracer), Combo (WL helpers) and
# MorphismGraph (circular(m)), all above.
include("circular/CircularGraph.jl")
include("circular/CircularGraphKey.jl")
include("circular/CircularGraphDisplay.jl")
include("circular/CircularMerge.jl")
include("circular/CircularFaces.jl")
include("circular/CircularConvert.jl")
export CircularNode, CircularGraph, circular, circular_node, arms, arm_colour, arm_count, circular_degree, rotate_arms
export circular_dot_may_pass, CIRCULAR_DOT_MIN_ARMS
export circular_dots_reducible
export CIRCULAR_DOT_POLICY
export general_braid_cluster, general_braid_cluster_swapped
export circular_key, circular_canonical_key, circular_boundary, circular_opposite_slot
export is_wired, merge_nodes, merge_at_edge, merge_all, show_circulargraph
export CircularMorphismGraph, circular_morphism_graph, circular_morphism

include("circular/CircularGraphIO.jl")
export circulargraph_to_string, circulargraph_from_string, save_circulargraph, load_circulargraph
export circularmorphismgraph_to_string, circularmorphismgraph_from_string
export save_circular_morphism, load_circular_morphism
export CircularMorphismEntry, save_circular_morphisms, load_circular_morphisms
export save_zamo_circular_morphisms

include("circular/CircularRegion.jl")
export circular_region_adjacency, circular_region_distances, circular_boundary_edge, circular_region_of_dot
export circular_region_distances_edges_only
export circular_region_tree_words
export circular_braid_class, circular_distance_word, circular_dot_reducible, circular_non_normal_dots
export circular_region_tree_parents
export CIRCULAR_WARN_NOT_NORMAL
export circular_node_sector_regions
export node_class, node_cap

include("circular/CircularBoundaryWord.jl")
export CircularBoundaryWord, CircularBoundaryLetter, CircularCorner
export circular_boundary_words, circular_region_boundary_word, is_interior

include("circular/CircularCombo.jl")
export CircularCombo, CircularComboR

include("circular/CircularDecorated.jl")
export CircularDecorated, circular_decorated, region_label

include("circular/CircularDecoratedIO.jl")
export soergelpoly_to_string, soergelpoly_from_string
export circulardecorated_to_string, circulardecorated_from_string
export save_circulardecorated, load_circulardecorated
export CLBasis, cl_basis_index, cl_basis_to_string, cl_basis_from_string
export save_cl_basis, load_cl_basis

# diagram/Join.jl — compose/tensor/flip and the mirrors, for MorphismGraph AND
# CircularMorphismGraph, hence after the circular types.
include("diagram/Join.jl")
include("diagram/JoinShift.jl")
include("diagram/JoinCircularMirror.jl")
export join_graphs, tensor, compose, flip, hflip, shift_right, shift_left
export vflip

include("morphism/BraidMoves.jl")
export apply_braid_move, replay_braid_moves, braid_to_end_with, braid_to_word
export same_element, is_reduced, canonical_word
export braid_move_sites, reduced_words, reex_graph

include("morphism/LightLeaves.jl")
export identity_strand, dot_morphism, light_leaf_up
export braid_move_morphism, braid_top, braid_top_to
export identity_morphism, merge_morphism, split_morphism, cap_morphism
export subexpressions, decorations, defect, expressed_word
export light_leaf, double_leaf, light_leaves, double_leaves
export a3_words

include("morphism/Zamolodchikov.jl")
export Zamo, zamo_words, zamo_lhs, zamo_rhs, zamo_relation_endpoints_agree

include("morphism/CircularLightLeaves.jl")
export circular_light_leaf, circular_double_leaf
export circular_light_leaves, circular_double_leaves
export circular_load_double_leaves, circular_reduce_double_leaf

include("morphism/DLBasis.jl")
export LLEntry, DLEntry, LLTable, DLTable
export ll, dl, dl_degrees, dl_ranks, dl_entries
export expected_ranks, ranks_match, is_dl_basis, dl_warnings!, clear_dl_cache!
export CircularDLEntry, CircularDLTable, CircularDLDecomposition
export circular_dl, circular_dl_matrix, circular_dl_triangular, is_circular_dl_basis

# circular/CircularExplosion.jl — when a computation explodes, save the diagram
# to disk. Pure diagnostics; needs `save_circulardecorated` from above.
include("circular/CircularExplosion.jl")
export CIRCULAR_EXPLOSION_DIR, CIRCULAR_EXPLOSION_ENABLED
export circular_capture_explosion, circular_save_explosion
export circular_reset_explosion, circular_explosions

# circular/rules/ — the rule files. The include order below is semantic; several
# files use helpers defined in the file directly above them.
include("circular/rules/CircularRules.jl")
include("circular/rules/CircularMergeRules.jl")
include("circular/rules/CircularBraidRules.jl")
include("circular/rules/CircularBraidChannels.jl")
include("circular/rules/CircularGen12Merge.jl")
export contract_circular_cluster
include("circular/rules/CircularDotOnGen12.jl")
include("circular/rules/CircularGen12Expand.jl")
export expand_gbraid, GBRAID_PREIMAGES
export gbraid_tower_decompositions, gbraid_tower_literal, gbraid_tower_cluster,
       expand_gbraid_tower, expand_gbraid_variants
include("circular/rules/CircularViaPreimage.jl")
export circular_via_preimage
include("circular/rules/CircularBraidOnGen12.jl")
include("circular/rules/CircularGen12OnGen12.jl")
# circular/rules/CircularGen12Null.jl — C22 `gen12_on_gen12_null`: two braid-like
# nodes sharing four or more edges are 0. Needs `_circular_glue3_rewire_at`
# (CircularBraidOnGen12.jl, above), whose surgery it runs to prove each hit.
include("circular/rules/CircularGen12Null.jl")
include("circular/rules/CircularGen12TwoEdges.jl")
include("circular/rules/CircularDotMerge.jl")
export _fr_trivalent_into_gen12, _fr_two_adjacent_dots
export _fr_dot_on_gen12_collapse
export _fr_mixed_into_gen12, _circular_split_mixed_at, _circular_mono_into_gen12_at
export CIRCULAR_MIXED_MERGE
export _fr_pitchfork, _fr_gen12_on_gen12_null
include("circular/rules/CircularWeight.jl")
export circular_weight, circular_arm_weight, circular_sep_weight, circular_chain_weight
export circular_separating_edges, circular_chain_arms
export CIRCULAR_TRIVALENT_MERGE, CIRCULAR_MVALENT_MERGE, CIRCULAR_WEIGHT_MODE
include("circular/rules/CircularDriver.jl")
export CircularRule, CIRCULAR_RULES, reduce_circular, reduce_circular_full, reduce_circular_combo
export CIRCULAR_WEIGHT_ASSERT_EXEMPT
export circular_combo_by_weight
export circular_generic_preimage, CIRCULAR_GENERIC_PREIMAGE_MAXSITES
include("circular/rules/CircularRegionRules.jl")
export CIRCULAR_REGION_RULES, CIRCULAR_RULES_REGION

include("circular/rules/ZamoRules.jl")
export zamo_rule_fixtures, zamo_anchor_roles, zamo_anchors, zamo_prefilter_kind
export find_zamo_matches, apply_zamo_rule

include("circular/rules/CircularDecoratedRules.jl")
include("circular/rules/CircularLeaveDriver.jl")
export CircularDecoratedMorphism, apply_circular_d4, find_circular_d4_match, reduce_to_circular_leave
export CIRCULAR_D4_WORD_TIEBREAK
export CIRCULAR_ZAMO_ENABLED

include("circular/CircularComponents.jl")
export circular_connected_components, circular_extract_floating_components

include("circular/rules/CircularD4Node.jl")
export apply_circular_d4_node, find_circular_d4_node_match, find_circular_d4_any
export with_marks
export circular_extract_scalars
export circular_fusion_step

include("morphism/DLBasisWriter.jl")
export DLCombination, in_circular_dl_basis
export dl_to_circular_matrix, circular_to_dl_matrix
export hflip_dl_matrix, vflip_dl_matrix
export compose_in_basis, structure_constants

include("morphism/CircularPairing.jl")
export CircularPairing, pairing

include("circular/rules/ZamoTermRules.jl")
export ZamoTermRule, zamo_term_rules, zamo_pass, apply_zamo_rule_combo

include("circular/rules/CircularZamoRegion.jl")
include("circular/rules/CircularZamoRegionStep.jl")
include("circular/rules/CircularZamoRegionInverse.jl")
export ZamoRegionPattern, zamo_region_pattern, zamo_region_signature
export circular_zamo_region_matches, circular_zamo_region_step
export zamo_triangles, zamo_regions, zamo_region_components
export zamo_inverse_terms
export circular_zamo_direction_step, zamo_nonzamo_rule_fires
export ZamoInverseLogEntry, zamo_inverse_log, zamo_inverse_log_clear!

include("circular/rules/CircularParallelMerge.jl")
export circular_parallel_merge_step, find_circular_parallel_merge

include("circular/rules/Circular2Parallel.jl")
include("circular/rules/Circular2ParallelApply.jl")
export circular_2parallel_step, find_circular_2parallel, circular_2parallel_apply

include("circular/rules/CircularRegionWord.jl")
export circular_unreduced_transition, circular_regionword_dihedral_step
export circular_region_words, circular_path_edges
export circular_region_mvalent_words, CIRCULAR_REGIONWORD_MVALENT
export circular_region_words_2k, CIRCULAR_REGIONWORD_2K
export CIRCULAR_REGIONWORD_DOT_POLICY


include("circular/rules/CircularLeafDegree.jl")
export circular_term_degree, circular_degree_violations, test1_graded, circular_leaf_degrees

include("circular/rules/CircularInverseRules.jl")
export expand_circular_merge
export expand_circular_braid_relation, expand_circular_braid_relation_variants
export expand_circular_d4


include("circular/rules/CircularRuleFixtures.jl")
export circular_braid_star, circular_glued_braids, circular_two_edge_fixture
export circular_bead_graph, circular_counter_joined, circular_mixed_at_braid

include("circular/rules/CircularDotSlide.jl")
export circular_dot_slide_step, CIRCULAR_DOT_SLIDE_ENABLED

include("morphism/CLBasis.jl")
export cl_cache_dir, cl_cache_path, cl_reduce, cl_basis,
       cl_basis_filter, cl_coordinates

include("circular/rules/CircularRexRelation.jl")
export circular_rex_relation, circular_rex_cache_dir, circular_rex_cache_path

include("circular/rules/CircularSplice.jl")
export circular_splice, circular_splice_identity, circular_splice_sides

include("circular/rules/CircularRexFusion.jl")
export circular_rex_fusion_step, CIRCULAR_REX_FUSION_ENABLED, CIRCULAR_MONO_BREAK

include("circular/rules/CircularDotFusion.jl")
export circular_dot_fusion_step, CIRCULAR_DOT_FUSION_ENABLED

# render/ — SVG renderers. Tutte is the default `display`; the circular
# renderers draw CircularGraph/CircularDecorated; CircularSteps must come last
# (CIRCULAR_STEPS collects all step functions defined above).
include("render/Morphism.jl")
export RedStep, reduction_with_moves, radial_svg, save_radial, linear_svg, save_linear

include("render/Graph.jl")
export diagram_svg, diagram_png, save_diagram
export notebook_display_sink!

include("render/Rect.jl")
export morphism_rect_svg, morphism_rect_png, rect, Rect

include("render/tutte/Avoid.jl")
include("render/tutte/Markers.jl")
include("render/tutte/Tutte.jl")
include("render/tutte/TutteSVG.jl")
export tutte_svg
include("render/tutte/Wiring.jl")
export debug_labels
export wiring_table, show_wiring
include("render/tutte/Display.jl")
export display_rect, display_tutte

include("render/CircularTutte.jl")
include("render/CircularTutteSVG.jl")
export circular_tutte_svg, circular_tutte_positions, display_circular_tutte, debug_circular, show_circular_step

include("render/Cells.jl")
include("render/CellLabels.jl")
include("render/CellDistance.jl")
export cell_number_svg, display_cell_number
export cell_adjacency, cell_distances, cell_distance_svg, display_cell_distance

include("render/CircularCells.jl")
export circular_decorated_svg, display_circular_decorated, show_circular_combo

include("render/CircularRegionTree.jl")

include("render/Decorated.jl")
export decorated_svg, display_decorated, display_combo

include("render/CircularSteps.jl")
export CIRCULAR_STEPS, circular_steps, circular_step_once, circular_trace, show_circular_step_row, show_combo

# words/ tools for the ε-reducibility exploration (independent of diagrams).
include("words/SequenceGraph.jl")
export SequenceGraph, build_reduction_graph, reducible_words, non_increasing_subgraph
export nodes, edges_from, edge_label, all_labels_nonincreasing
export BuildState, run_build!, frontier_remaining
export cached_graph, clear_graph_cache!, graph_bytes
export down_graph, reaches_empty_down, down_reachability_report, down_counts_by_length
export reduces_to_empty_down, down_path
export expands_then_reduces, find_expand_then_reduce
export all_circular_words, expandable_all_around, counts_by_length
export up_component, UP_SEEDS

include("words/Checkpoint.jl")
export save_build, load_build, build_or_resume
export save_graph, load_graph
export save_wordset, load_wordset, up_component_cached

# ---- precompilation ---------------------------------------------------------
#
# One representative run of the two hot paths, so a fresh session does not pay for
# compiling the rule driver and the renderer on its first real call. Kept cheap on
# purpose: deriving the Zamolodchikov relation (`cl_reduce(...; zamo = true)`) costs
# tens of seconds and is left to the first caller who asks for it.
@setup_workload begin
    @compile_workload begin
        g = circular(unit_graph())
        reduce_circular(g)
        fm = circular_morphism(morphism_graph(unit_graph(), 0, 0))
        fdm = CircularDecoratedMorphism(fm, fill(one(SoergelPoly),
                                                region_count(fm.graph)))
        reduce_to_circular_leave(fdm)
        diagram_svg(unit_graph())
        circular_tutte_svg(g)
        for d in double_leaves([1, 2, 1], [1, 2, 1])
            d.degree <= 0 || continue
            cl_reduce(d.morphism; zamo = false)
        end
    end
end

end # module
