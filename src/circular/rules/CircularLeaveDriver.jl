# circular/rules/CircularLeaveDriver.jl — the DRIVER of the circular-leaf reduction.
#
# The file is included directly after CircularDecoratedRules.jl, so every name it uses (`apply_circular_d4`,
# `find_circular_d4_any`, `circular_extract_scalars`, …) is already defined.
#
# What lives here: `reduce_to_circular_leave` (the round-based driver), its growth
# guard `_circular_leave_growth_guard`, the position switch `_CIRCULAR_2PARALLEL_SLOT`, the
# finishing step `_circular_leave_finish` and `circular_fusion_step`.

"""
    _circular_leave_growth_guard(trail, limit, fdm)

The GROWTH GUARD of the circular-leave driver: a safety valve, not a termination
proof, for the case that the "connect → reduce → search again" loop
(2parallel/P1) fails to terminate. The message names the boundary word and arm
sequences so the case can be reproduced in a notebook.

`trail` is the node count **along one recursion path** (each branch gets its
own copy, siblings don't mix). `limit` counts **increases**, not entries:
`limit + 1` entries are needed, and only once ALL of them **strictly
increase** does the driver abort.

An abort does NOT mean "non-terminating", only "growing suspiciously
monotonically for a while" — raising `limit` is a legitimate response.
"""
function _circular_leave_growth_guard(trail::Vector{Int}, limit::Int,
                                 fdm::CircularDecoratedMorphism)
    limit <= 0 && return nothing
    length(trail) >= limit + 1 || return nothing
    tail = trail[(end - limit):end]
    all(k -> tail[k] < tail[k + 1], 1:limit) || return nothing
    g = fdm.m.graph
    trail = join(tail, " -> ")
    arm_lists = [arms(nd) for nd in g.nodes]
    # Save the diagram itself, not only its boundary-word summary
    # (circular/CircularExplosion.jl). If the write fails, that's irrelevant
    # to the diagnosis — the error below must still be raised.
    path = try
        CIRCULAR_EXPLOSION_ENABLED[] ?
            circular_save_explosion(fdm, "growth guard: $(trail)";
                                    depth = length(trail) - 1, name = "growthguard") : ""
    catch
        ""
    end
    error("reduce_to_circular_leave: node count has grown monotonically for " *
          "$(limit) steps ($(trail)) — the loop presumably does not " *
          "terminate. Boundary word $(g.word), arms $(arm_lists). To keep " *
          "computing, raise `node_growth_limit` or set it to 0 (off)." *
          (isempty(path) ? "" : " Diagram saved to $(path)."))
end

"""
    reduce_to_circular_leave(fdm::CircularDecoratedMorphism; maxrounds = 100, maxsteps = 1000)
        -> CircularComboR

Brings a decorated circular morphism into circular-leaf form. Alternates (a)
reducing with [`reduce_circular`](@ref) (`maxsteps` as an emergency brake,
termination backed by the active `circular_weight` assert in CircularDriver.jl), then
(b) searching EVERY resulting term for a D4 match ([`find_circular_d4_any`](@ref) —
over BOTH variants, the edge and the node variant, in the joint candidate
search) and firing the matching rule there — until NO term has a D4 match left
or `maxrounds` is reached. Before every `reduce_circular` round,
[`circular_extract_scalars`](@ref) runs (analogous to `reduce_decorated`'s internal
`extract_scalars`), so labels on distance-0 cells move into the coefficient
instead of staying as a cell label (`reduce_circular`/`CIRCULAR_RULES` know no marking,
see the `circular_extract_scalars` docstring).

TERMINATION (not proven, as in the original — `maxrounds`/`maxsteps` are
emergency brakes): a plausible measure is `Σ_dots dist(cell)` — every D4 step
consumes exactly one dot in a cell of distance `d` against a neighbour cell of
strictly smaller distance `d' < d` (exactly `find_circular_d4_match`'s criterion);
the newly created dot terms (2/3) sit at the same or a cell closer to the
marking (splitting the BC cell can only change the distance non-increasingly,
never larger, since the new connecting edge runs inside the same old BC cell).
`find_circular_d4_match`'s ordering (largest cell distance first) keeps this
measure falling in a controlled way. There is no active `@assert` for it
(unlike `circular_weight` in CircularDriver.jl) — unlike the nine base rules, where the
weight formula is proven for ALL rules in advance, "Σ dist strictly
decreases" here is a heuristic without a complete proof over all cell
topologies; an `@assert` would therefore only generate noise, not a real
safety net.

The NODE variant fits the same measure: it too consumes exactly one dot whose
region has distance `d`, against a region of strictly smaller distance
(`find_circular_d4_node_match` requires `dist[R_between] < dist[R]`), and the newly
created dots of terms 2/3 sit at the far end of an arm bordering that closer
region.
"""
# `circular_extract_scalars` runs AFTER `reduce_circular` TOO, not only before every
# `reduce_circular` round: labels that CIRCULAR_RULES themselves produce (barbell ↦ α_c
# into a region, constants from merges) would otherwise stay as a cell label when no
# further step followed — two otherwise-equal circular leaves could then differ only in
# whether a scalar sits in the coefficient or in the label (a different
# `circular_canonical_key`). So symmetrically: after every
# `reduce_circular` call, each term is run once through `circular_extract_scalars`. The
# cuts come from the starting morphism (`reduce_circular` does not change the boundary).
function _circular_extract_after(combo, m::CircularMorphismGraph)
    out = CircularComboR()
    for (fd, coeff) in pairs_of(combo)
        m2 = CircularMorphismGraph(fd.graph, m.cut1, m.cut2)
        f, ddm = circular_extract_scalars(
            CircularDecoratedMorphism(m2, fd.region_labels, fd.outer_label))
        out = out + (coeff * f) * CircularCombo{SoergelPoly}(_circular_decorated(ddm))
    end
    return out
end

"""
    _CIRCULAR_2PARALLEL_SLOT

WHERE 2parallel SITS IN THE DRIVER. Three positions, in driver order:

* `:before_zamo` — directly after P1, BEFORE Zamo and D4;
* `:after_zamo` — after BOTH Zamo stages (whole-term and region variant), but
  BEFORE the D4 check;
* **`:after_d4` (current default)** — in the `match === nothing` branch of the
  D4 check, i.e. the last attempt before a term counts as finished.

**Why `:after_d4`**: 2parallel GENERATES terms (a two-term sum with `α_s/2`)
and grows the diagram, so it belongs after every step that simplifies without
growing the term count (`CIRCULAR_RULES` incl. F2-3, P1, Zamo). The other two
positions stay reachable through this constant, for checking whether a
repositioning shifts a decomposition.
"""
const _CIRCULAR_2PARALLEL_SLOT = Ref(:after_d4)

function reduce_to_circular_leave(fdm::CircularDecoratedMorphism; maxrounds::Int = 100,
                             maxsteps::Int = 1000,
                             trace = nothing, depth::Int = 0,
                             node_trail::Vector{Int} = Int[],
                             node_growth_limit::Int = 15,
                             after_inverse::Bool = false,
                             rex_alpha_blocked::Bool = false)
    # `after_inverse` (direction choice): this call processes the
    # result of a BACKWARD Zamo step (Z1⁻¹/Z2⁻¹). In its FIRST round no Zamo
    # step then fires (neither whole-term nor region variant, neither forward
    # nor backward) — otherwise forward Z1/Z2 would immediately undo the step
    # (infinite loop). The backward result is accepted only because a
    # non-Zamo rule fires on every one of its terms (`circular_zamo_direction_step`,
    # CircularZamoRegion.jl) — guaranteed to come up first in round 1; from round 2
    # on, the normal direction choice applies again.
    # GROWTH GUARD (see `_circular_leave_growth_guard`). `node_trail` is the node
    # count along THIS recursion path; every recursive call gets its own copy
    # (`_trail`) so sibling branches don't mix.
    _trail = vcat(node_trail, length(fdm.m.graph.nodes))
    # EXPLOSION RECORDING (circular/CircularExplosion.jl):
    # remember the deepest diagram reached, so `circular_capture_explosion`
    # can write it to disk after a StackOverflow/assert. Costs one integer
    # comparison per call.
    _circular_note_depth(fdm, depth)
    _circular_leave_growth_guard(_trail, node_growth_limit, fdm)
    # `trace`: if set, `trace(depth, phase, round, combo)` is called at every
    # intermediate state (`nothing`, the default, costs only this one
    # comparison per call).
    _tr(phase, rnd, c) = trace === nothing ? nothing : trace(depth, phase, rnd, c)
    # FLOATING COMPONENTS FIRST: a connected component
    # hanging neither off a leaf nor off another component (can arise from
    # `compose`) is evaluated ONCE as a whole in isolation and its polynomial
    # entered directly into the geometrically correct region — before any of
    # the incremental `CIRCULAR_RULES` sees it node by node and would have to carry
    # the region bookkeeping across many individual steps
    # (circular/CircularComponents.jl, `circular_extract_floating_components`). Sits BEFORE
    # `circular_extract_scalars`, because a floating component can itself carry
    # non-trivial region labels that only arise once the component is
    # removed.
    fdm_fc, factor_fc = circular_extract_floating_components(fdm)
    # `circular_extract_scalars` BEFORE every `reduce_circular` round: `reduce_circular`/
    # `CIRCULAR_RULES` are not marking-aware (see the `circular_extract_scalars`
    # docstring) — without this step e.g. the alpha_i label of the D4 trivalent
    # term would stay on its own (distance-0) cell instead of moving into the
    # coefficient.
    factor0, fdm0 = circular_extract_scalars(fdm_fc)
    factor0 = factor_fc * factor0

    # FUSION FIRST: `fusion_step` sits BEFORE the plain fallback. As long as a
    # non-constant label sits on a region of distance >= 1, it is pushed along the distance
    # gradient towards the marking; at distance 0 `circular_extract_scalars`
    # then pulls it out. Without this step, an `alpha_i` set by D4 would stay
    # put forever.
    fus = circular_fusion_step(fdm0)
    if fus !== nothing
        out = CircularComboR()
        for (dd, c2) in pairs_of(fus)
            m2 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
            out = out + c2 * reduce_to_circular_leave(
                CircularDecoratedMorphism(m2, dd.region_labels, dd.outer_label);
                maxrounds = maxrounds, maxsteps = maxsteps - 1,
                trace = trace, depth = depth + 1,
                node_trail = _trail, node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
        end
        return factor0 * out
    end

    combo, _ = reduce_circular(_circular_decorated(fdm0); maxsteps = maxsteps)
    combo = factor0 * _circular_extract_after(combo, fdm.m)
    _tr(:after_reduce_circular, 0, combo)
    for rnd in 1:maxrounds
        progressed = false
        out = CircularComboR()
        for (fd, coeff) in pairs_of(combo)
            m2 = CircularMorphismGraph(fd.graph, fdm.m.cut1, fdm.m.cut2)
            # ---- PHASE A AT THE HEAD OF EVERY ROUND -------------------------
            # "as soon as a diagram contains polynomials, the polynomials
            # should first be moved towards the marking" — not inside
            # `reduce_circular` (not marking-aware, has no cuts), but one level up,
            # HERE. Without this check, terms falling out of `reduce_circular`
            # would go UNCHECKED into P1 → 2parallel → Zamo → D4. Example: the
            # double leaf `x=121321, e=100000, y=321323, f=101100` carries `α₂`
            # on a distance-1 region, `reduce_circular` fires `[:merge, :merge]`,
            # and one resulting term still carries the `α₂` — the throw class
            # "the dot's own region carries a non-trivial label". With the check
            # the invariant holds without gaps: NO graph step ever sees a
            # non-trivial label on an inner region.
            #
            # `reduce_circular` produces no labels itself
            # (`_circular_transfer_labels` only moves them, `_fr_barbell` puts
            # its α straight into the coefficient) and `_circular_extract_after`
            # already pulls the scalars, so on ordinary material the check is a
            # no-op; a rule that started producing a label would otherwise fall
            # through silently to the D4 matcher. The `any(!isone, …)` guard
            # keeps it free of cost: all three steps require a non-trivial
            # label, and without the guard every term and round would run extra
            # `circular_region_distances` BFS passes for nothing.
            factor = one(SoergelPoly)
            fdm2 = CircularDecoratedMorphism(m2, fd.region_labels, fd.outer_label)
            fus = nothing
            if any(!isone, fd.region_labels)
                factor, fdm2 = circular_extract_scalars(fdm2)
                fus = circular_fusion_step(fdm2)
            end
            if fus !== nothing
                progressed = true
                for (dd, c2) in pairs_of(fus)
                    m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                    out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                        CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                        maxrounds = maxrounds, maxsteps = maxsteps - 1,
                        trace = trace, depth = depth + 1,
                        node_trail = _trail, node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
                end
                continue
            end
            # From here on use ONLY `fdm2`/`fd2`, never `fd` again — otherwise
            # a scalar just pulled into `factor` would come along twice (once
            # as a factor, once as a label) or not at all.
            fd2 = _circular_decorated(fdm2)
            # P1 BEFORE THE D4 CHECK. Two parallel 1/3 edges are the same
            # morphism as an `[1,1,3,3]` node, but a different
            # `circular_canonical_key` — without this step two circular leaves
            # are counted where there should be one. Sits here,
            # not in CIRCULAR_RULES, because the rule INCREASES `circular_weight`; its
            # termination measure is the falling edge count. See the P1
            # section below.
            par = circular_parallel_merge_step(fdm2)
            if par !== nothing
                progressed = true
                for (dd, c2) in pairs_of(par)
                    m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                    out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                        CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                        maxrounds = maxrounds, maxsteps = maxsteps - 1,
                        trace = trace, depth = depth + 1,
                        node_trail = _trail, node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
                end
                continue
            end
            # ---- 2parallel, EARLIEST POSITION (measurement tool only) -------
            # `:after_d4` (the current default, see `_CIRCULAR_2PARALLEL_SLOT`) is
            # used instead; this position stays reachable via
            # `_CIRCULAR_2PARALLEL_SLOT[] = :before_zamo` purely for measurement.
            if _CIRCULAR_2PARALLEL_SLOT[] === :before_zamo
                same = isone(fdm2.outer_label) ? circular_2parallel_step(fdm2) : nothing
                if same !== nothing
                    progressed = true
                    for (dd, c2) in pairs_of(same)
                        m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                        out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                            CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                            maxrounds = maxrounds, maxsteps = maxsteps - 1,
                            trace = trace, depth = depth + 1,
                            node_trail = _trail,
                            node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
                    end
                    continue
                end
            end
            # ---- ZAMO -----------------------------------------------------------
            # The whole-term rules Z1/Z2 (rules/ZamoTermRules.jl): after the
            # parallel-edge steps, BEFORE D4. Not in CIRCULAR_RULES (the
            # `circular_weight` assert forbids 7 ↦ 7 nodes). Termination: the exact
            # morphism-key matching never hits its own output (checked when
            # the rule was built) — finitely many firings per branch, then the
            # weight-decreasing steps take over again. The boundary-word
            # prefilter in `circular_zamo_step` keeps the step free for ordinary
            # terms (the expensive Z2 derivation only runs when a term
            # carries one of the two fixed LHS boundary words). Resolved at
            # RUNTIME — ZamoTermRules.jl is included after this file (same
            # forward reference as C18/CIRCULAR_RULES).
            # OUTER-LABEL GUARD. Both Zamo variants
            # rebuild their right-hand sides from GRAPHS
            # (`circular_zamo_step` via `circular_decorated(_zamo_denorm(...))`,
            # `circular_zamo_direction_step` gets only `fd2.graph`) — an outer
            # label could only silently drop out there. So, like the existing
            # `all(isone, region_labels)` guard: if the term carries a
            # polynomial outside, Zamo does not fire. An outer label only arises
            # from `_fr_barbell`/floating components, never from Zamo clusters.
            zamo_erlaubt = !(after_inverse && rnd == 1) && isone(fd2.outer_label) &&
                           CIRCULAR_ZAMO_ENABLED[]
            zam = zamo_erlaubt ? circular_zamo_step(fd2, fdm.m.cut1) : nothing
            if zam !== nothing
                progressed = true
                _tr(:zamo_ganzterm, rnd, zam)
                for (dd, c2) in pairs_of(zam)
                    m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                    out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                        CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                        maxrounds = maxrounds, maxsteps = maxsteps - 1,
                        trace = trace, depth = depth + 1,
                        node_trail = _trail, node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
                end
                continue
            end
            # ---- ZAMO, REGION VARIANT (stage 2) -------------------------------
            # `circular_zamo_step` matches only the EXACT whole-term boundary word
            # of the two 7-node fixtures, so it misses terms that structurally
            # contain a Zamo cluster somewhere (6 inner triangular regions with
            # boundary word {1,2,3}, see `circular_zamo_region_matches`). This
            # variant searches LOCALLY over the region structure, independent of
            # term size, so it also fires in large/merged terms.
            # `circular_zamo_region_step(fd2.graph)` only applies to the
            # restricted case (each of the 7 cluster nodes has exactly the
            # fixture degree, see the CircularZamoRegion.jl file header) —
            # otherwise `nothing`, no error.
            #
            # COST PREFILTER: `circular_zamo_region_matches` calls the full
            # region tracer even when nothing ends up matching. A Zamo cluster
            # needs at least 7 nodes in the host;
            # `length(fd2.graph.nodes) < 7` excludes most terms WITHOUT a
            # tracer call (analogous to the boundary-word prefilter of
            # `circular_zamo_step`).
            #
            # DIRECTION CHOICE: forward Z1/Z2 first; if nothing
            # fires and the term has ≥ 7 `zamo_regions`, then Z1⁻¹/Z2⁻¹ —
            # accepted only if a non-Zamo rule then fires, otherwise a fixed
            # point with a log (`zamo_inverse_log`). All inside
            # `circular_zamo_direction_step` (CircularZamoRegion.jl, see the box there);
            # the recursive call on a BACKWARD result gets
            # `after_inverse = true` (guard 2, see the head of this function).
            zamr = nothing
            if zamo_erlaubt && all(isone, fd2.region_labels) && length(fd2.graph.nodes) >= 7
                zamr = circular_zamo_direction_step(fd2.graph, fdm.m.cut1, fdm.m.cut2)
            end
            if zamr !== nothing
                zcombo, zinv = zamr
                progressed = true
                _tr(zinv ? :zamo_backward : :zamo_forward, rnd, zcombo)
                for (dd, c2) in pairs_of(zcombo)
                    m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                    out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                        CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                        maxrounds = maxrounds, maxsteps = maxsteps - 1,
                        trace = trace, depth = depth + 1,
                        node_trail = _trail, node_growth_limit = node_growth_limit,
                        rex_alpha_blocked = rex_alpha_blocked,
                        after_inverse = zinv)
                end
                continue
            end
            # ---- 2parallel: AFTER ZAMO, BEFORE D4 ---------------------------
            # By this point it is certain: no CIRCULAR_RULES rule fires (incl.
            # F2-3), not P1, neither Zamo stage. Only THEN may the step fire
            # that GENERATES terms and initially grows the diagram
            # (two-term sum with `α_s/2`).
            # `isone(fdm2.outer_label)`: `circular_2parallel_step` rebuilds its
            # terms from graphs and doesn't know the outer label — the same
            # loud refusal as with Zamo.
            # Position switchable via `_CIRCULAR_2PARALLEL_SLOT` (see its docstring).
            if _CIRCULAR_2PARALLEL_SLOT[] === :after_zamo
                nz = isone(fdm2.outer_label) ? circular_2parallel_step(fdm2) : nothing
                if nz !== nothing
                    progressed = true
                    for (dd, c2) in pairs_of(nz)
                        m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                        out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                            CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                            maxrounds = maxrounds, maxsteps = maxsteps - 1,
                            trace = trace, depth = depth + 1,
                            node_trail = _trail,
                            node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
                    end
                    continue
                end
            end
            # The pattern of two braids of the same colour pairing with
            # exactly two connecting edges is NOT computed away here: it also
            # occurs in genuine double leaves, and a rule that removed it would
            # alter the DL inventory itself. The double connection stays.
            match = find_circular_d4_any(m2)
            if match === nothing
                # ---- 2parallel ONLY AFTER D4 ----------------------------------
                # Checked only at the very end, as the last thing tried
                # before a term counts as done (current default position, see
                # `_CIRCULAR_2PARALLEL_SLOT`). Two parallel edges of the SAME
                # colour at an INNER region with distance pattern
                # (i−1, i, i+1) become a 4-armed node — a two-term sum with
                # α_s/2. "Simplify the polynomials first" is handled by the
                # recursive call itself: it starts with `circular_extract_scalars`
                # and `circular_fusion_step`. Like P1 it does not sit in
                # CIRCULAR_RULES, because it INCREASES `circular_weight`; its own
                # termination measure is in the docstring of Circular2Parallel.jl
                # and is independent of this ordering (it only counts `s`
                # edges, not steps of other rules).
                # `isone(fdm2.outer_label)`: `circular_2parallel_step`
                # (circular/rules/Circular2Parallel.jl) rebuilds its terms from graphs
                # and doesn't know the outer label — the same loud refusal as
                # with Zamo.
                same = _CIRCULAR_2PARALLEL_SLOT[] === :after_d4 && isone(fdm2.outer_label) ?
                       circular_2parallel_step(fdm2) : nothing
                if same !== nothing
                    progressed = true
                    for (dd, c2) in pairs_of(same)
                        m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                        out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                            CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                            maxrounds = maxrounds, maxsteps = maxsteps - 1,
                            trace = trace, depth = depth + 1,
                            node_trail = _trail,
                            node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
                    end
                    continue
                end
                # ---- REX-FUSION, LAST IN THE GROUP ---------------------------
                # All other rules + Zamo run first, then the parallel rule and
                # the new dot rule — so it sits in the same `:after_d4` group,
                # AFTER 2parallel, because it is the expensive one (derived
                # relation, disk cache) and should only run where the cheap
                # rules find nothing. The trigger is the REGION WORD: it fires
                # exactly where some region word is not reduced, which the other
                # rules do not look at.
                # `isone(fdm2.outer_label)`: the step is label-free and checks
                # this itself, like Zamo/2parallel.
                # Runtime forward reference (CircularRexFusion.jl is included
                # after this file), like `circular_dot_slide_step` above.
                # ⚠ THE α-GUARD (`rex_alpha_blocked`). At `id_1212`, of the
                # four terms the step delivers, three are harmless
                # (two carry no unreduced region word at all, one drops
                # from `1212` to `12`) — but the term with coefficient `−α₁`
                # carries the region word `1212` AGAIN, so the max region word
                # does NOT fall and the rule would call itself forever. What DOES
                # fall there is the degree of the DIAGRAM: the term is
                # `−α₁ ⊗ (diagram of degree −2)`. So the measure is
                # lexicographic — (diagram degree, max unreduced region word) —
                # and a term that picked up an α is handed to the other rules
                # instead: below it the rex fusion stays off for the whole
                # subtree. ⚠ The α arrives in TWO shapes and both must be seen:
                # as a COEFFICIENT (route `:relation`, the `−α₁` term) and as a
                # REGION LABEL `α_s/2` (route `:fusion`, where the surgery writes
                # it into the new region). Watching only the coefficient
                # overflows the stack at `id_11`. Same shape as the
                # `after_inverse` guard against Z1/Z1⁻¹ ping-pong (head of this
                # function).
                rx = CIRCULAR_REX_FUSION_ENABLED[] && !rex_alpha_blocked &&
                     isone(fdm2.outer_label) ? circular_rex_fusion_step(fdm2) : nothing
                if rx !== nothing
                    progressed = true
                    for (dd, c2) in pairs_of(rx.combo)
                        m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                        out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                            CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                            maxrounds = maxrounds, maxsteps = maxsteps - 1,
                            trace = trace, depth = depth + 1,
                            node_trail = _trail,
                            node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked || degree(c2) > 0 ||
                                                any(!isone, dd.region_labels) ||
                                                !isone(dd.outer_label))
                    end
                    continue
                end
                # ---- DOT FUSION, after the rex fusion -----------------------
                # Where the rex fusion declines, the dots get their turn: the
                # shortest `:dot` transition — a dot whose reduced region word
                # ends on the dot's colour, under the dot policy
                # (`circular_dots_reducible`) — is rewritten by the derived
                # relation of `circular_dot_fusion_step`. One dot per step, the
                # terms run on recursively; the same α guard as for the rex
                # fusion applies. `:direct` dots are the classic D4's and the
                # step declines them itself. Runtime forward reference:
                # CircularDotFusion.jl is included after this file.
                df = CIRCULAR_DOT_FUSION_ENABLED[] && !rex_alpha_blocked &&
                     isone(fdm2.outer_label) ? circular_dot_fusion_step(fdm2) : nothing
                if df !== nothing
                    progressed = true
                    for (dd, c2) in pairs_of(df.combo)
                        m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                        out = out + (coeff * factor * c2) * reduce_to_circular_leave(
                            CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label);
                            maxrounds = maxrounds, maxsteps = maxsteps - 1,
                            trace = trace, depth = depth + 1,
                            node_trail = _trail,
                            node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked || degree(c2) > 0 ||
                                                any(!isone, dd.region_labels) ||
                                                !isone(dd.outer_label))
                    end
                    continue
                end
                out = out + (coeff * factor) * CircularCombo{SoergelPoly}(fd2)
                continue
            end
            progressed = true
            hit = _apply_circular_d4_hit(fd2, match)
            hit === nothing && error(
                "reduce_to_circular_leave: find_circular_d4_any found a match of " *
                "variant :$(match.kind) that the rule does not confirm — " *
                "matcher/rule are inconsistent")
            for (dd, c2) in pairs_of(hit)
                m3 = CircularMorphismGraph(dd.graph, fdm.m.cut1, fdm.m.cut2)
                # ⚠ A NAME OF ITS OWN. `factor` above is the scalar pulled out of the
                # term BEFORE the D4 round; this one is pulled out of the D4 result.
                # Both belong in the product, and reusing one name would multiply the
                # inner one twice and drop the outer one — the coefficient comes out
                # squared.
                inner_factor, ddm = circular_extract_scalars(
                    CircularDecoratedMorphism(m3, dd.region_labels, dd.outer_label))
                # FUSION HERE TOO. At the very top of this function all
                # labels are still 1, so fusion never fires there. `alpha_i`
                # only arises IN the D4 rounds — right here it must be pushed
                # outward along the distance gradient, otherwise it stays put
                # and the next D4 step throws over it.
                fus2 = circular_fusion_step(ddm)
                if fus2 !== nothing
                    for (ee, c3) in pairs_of(fus2)
                        m4 = CircularMorphismGraph(ee.graph, fdm.m.cut1, fdm.m.cut2)
                        out = out + (coeff * factor * c2 * inner_factor * c3) * reduce_to_circular_leave(
                            CircularDecoratedMorphism(m4, ee.region_labels, ee.outer_label);
                            maxrounds = maxrounds, maxsteps = maxsteps - 1,
                            trace = trace, depth = depth + 1,
                            node_trail = _trail,
                            node_growth_limit = node_growth_limit,
                            rex_alpha_blocked = rex_alpha_blocked)
                    end
                else
                    dd2, _ = reduce_circular(_circular_decorated(ddm); maxsteps = maxsteps)
                    out = out + (coeff * factor * c2 * inner_factor) *
                                _circular_extract_after(dd2, fdm.m)
                end
            end
        end
        combo = out
        _tr(:d4_round, rnd, combo)
        # After a backward Zamo, no Zamo was allowed in round 1; if nothing
        # else fired there either (only possible if `reduce_circular` at the head
        # already changed the term), the terms get the normal direction
        # choice in round 2 instead of coming back as "done".
        progressed || (after_inverse && rnd == 1) || return _circular_leave_finish(combo, maxsteps, depth)
    end
    return _circular_leave_finish(combo, maxsteps, depth)
end

"""
    _circular_leave_finish(combo::CircularComboR, maxsteps::Int, depth::Int) -> CircularComboR

Final structural pass over the terms `reduce_to_circular_leave` is about to return.

Without this pass the driver can return a term that `CIRCULAR_RULES` still
matches — e.g. two same-coloured nodes joined by TWO edges, on which C7
`mono_double` fires and the term is **0**. Such a term is not a normal form,
and counting it as a circular leaf inflates every leaf count.

The structural pass (`reduce_circular`) runs ONCE before the round loop, and
every step inside the loop that changes a graph recurses and so re-enters that
pass. But the loop can also fall out — `progressed == false`, or `maxrounds`
used up — and then holds terms that were never re-examined. `reduce_circular`
itself stops silently when its own `maxsteps` budget is gone (`for _ in
1:maxsteps`, no error), so an exhausted budget deep in the recursion is
invisible too.

This pass closes both: before returning, every term is offered to
`CIRCULAR_RULES` once more, and the result replaces it. Terms that are already
normal forms pass through unchanged (one failed match attempt each), and a term
that still has a rule gets it, including one that reduces to zero.
"""
function _circular_leave_finish(combo::CircularComboR, maxsteps::Int, depth::Int)
    # TOP LEVEL ONLY.  Inner levels hand their terms upward, so whatever they
    # miss is still in the combo the outermost call returns — checking there is
    # enough, and checking at every level would re-reduce the same terms once
    # per recursion level.
    depth == 0 || return combo
    out = CircularComboR()
    for (dd, coeff) in pairs_of(combo)
        red, hist = reduce_circular(dd; maxsteps = max(maxsteps, 1))
        out = out + coeff * (isempty(hist) ? CircularComboR(dd) : red)
    end
    return out
end

"""
    circular_fusion_step(fdm::CircularDecoratedMorphism) -> Union{Nothing, CircularComboR}

**Exact port of `fusion_step`** onto REGIONS. A
guided fusion step along the marking distances: find a region `R` with a
non-constant label `f` (`degree(f) >= 1`) and distance `>= 1`, plus a
neighbour region `R2` with `dist[R2] == dist[R] - 1` (via
`circular_region_adjacency`); fuse the separating edge with f-side = `R`. Fires
ONLY with a marking set — without a marking all distances are `-1` and the
result is `nothing`.

**Why this is needed.** `CIRCULAR_RULES` contains no fusion (see the
CircularDriver.jl header). Without this step, an `alpha_i` set by D4 could
never migrate outward — it would stay on its own inner region, and the next
D4 step would throw "the dot's own region carries a non-trivial label".

The mechanism: `f` migrates through an edge to `act(i,f)`
(term 1), and the broken-open strand with two fresh dots carries
`demazure(i,f)` (term 2) — both already implemented in
[`_circular_fuse_edge`](@ref).
"""
function circular_fusion_step(fdm::CircularDecoratedMorphism)
    g = fdm.m.graph
    # EDGE DISTANCE, not the node-weighted one. There is ONE mechanism here —
    # push polynomials stepwise through edges — and it needs a distance that
    # drops by exactly 1 per step, which is what
    # `circular_region_distances_edges_only` provides: it counts only real edges,
    # i.e. exactly the steps this rule can take. The node-weighted variant
    # additionally counts node-adjacencies at `:gen2`/`:gen13` (weight =
    # crossed colours) and compresses distances by doing so — at the star
    # node `[1,1,3,3,1,1,3,3]` it leaves no neighbour with `dist − 1` at all,
    # even though the edge is there, so fusion cannot fire.
    dist = circular_region_distances_edges_only(fdm.m)
    all(==(-1), dist) && return nothing              # no marking
    # The boundary circle is a WALL: the fusion never pushes a
    # label out of the disc. It only walks real edges towards the marking; a
    # label that already sits outside is handled by `circular_extract_scalars`, which
    # sees the outer region at distance 0 and turns it into a coefficient.
    adj, _ = circular_region_adjacency(g)
    fd = _circular_decorated(fdm)
    for R in 1:region_count(g)
        dist[R] >= 1 || continue
        f = fdm.region_labels[R]
        degree(f) >= 1 || continue                   # only non-constant labels
        for (R2, ei) in sort(adj[R])
            dist[R2] == dist[R] - 1 || continue
            # The regions come DIRECTLY from the adjacency — no detour via
            # faces. `circular_region_adjacency` already knows the pair (R | R2)
            # and the separating edge; the face->region translation in the
            # `fside` variant of `_circular_fuse_edge` would be ambiguous here
            # anyway (in the throw case ALL 7 regions lie in the same face).
            # f-side is `R`, the region carrying the polynomial.
            return _circular_fuse_edge_regions(fd, ei, R, R2)
        end
    end
    return nothing
end
