# circular/CircularExplosion.jl — WHEN A COMPUTATION EXPLODES, KEEP THE DIAGRAM.
#
# Whenever a diagram blows up during some computation, that is caught and the
# diagram is saved to disk.
#
# THE PROBLEM. A run that blows up (`StackOverflowError` from an unbounded
# recursion, the growth guard, the `circular_weight` assert) tells us the message
# and nothing else — the diagram that caused it is gone with the stack, and the
# case has to be reconstructed by hand from the boundary word. That is the single
# most expensive step in every termination question.
#
# WHAT THIS DOES. `reduce_to_circular_leave` reports every diagram it works on to
# `_circular_note_depth` (two assignments per call). Wrap a computation in
# `circular_capture_explosion`, and any error writes the LAST diagram — and, if
# it is a different one, the DEEPEST by recursion depth — to disk as a `.cd`
# file plus a `.txt` sidecar with the reason, the depth, the marking, the
# boundary word and the switch position: enough to reload it in a notebook with
# `load_circulardecorated` and look at it.
#
# WHY THE LAST ONE. On a `StackOverflowError` the stack unwinds before any catch
# runs, so „the frame that overflowed" cannot be read out; the last diagram
# recorded before the unwinding is the closest honest approximation. The
# `depth` argument is NOT a reliable stand-in for it: something inside the
# reduction (the rex relation derivation) starts a FRESH top-level
# `reduce_to_circular_leave`, so `depth` restarts at 0 there — an overflow can
# arrive at a shallow recorded depth despite thousands of calls having been
# made. Both diagrams are written; the sidecar names the call count.
#
# ⚠ IT IS A DIAGNOSTIC, NOT A GUARD. Nothing here changes what is computed, and
# nothing here catches anything on its own — the error is re-thrown after the
# file is written.

"""
    CIRCULAR_EXPLOSION_DIR

Where [`circular_capture_explosion`](@ref) writes.

Nothing is written inside the package tree. The default is `DH_EXPLOSION_DIR`
if that environment variable is set, otherwise a `diagrammatichecke-explosions`
folder under `tempdir()`.

To collect the dumps somewhere permanent, set either the environment variable
or the `Ref` itself before running:

    CIRCULAR_EXPLOSION_DIR[] = "/path/to/your/sandbox/explosions"
"""
const CIRCULAR_EXPLOSION_DIR =
    Ref{String}(get(ENV, "DH_EXPLOSION_DIR",
                    joinpath(tempdir(), "diagrammatichecke-explosions")))

"""
    CIRCULAR_EXPLOSION_ENABLED

Is the depth bookkeeping on? Default `true`; the cost is one integer comparison
per `reduce_to_circular_leave` call. Set to `false` for a timing run.
"""
const CIRCULAR_EXPLOSION_ENABLED = Ref(true)

# The deepest diagram of the current computation: `(depth, fdm)` or `nothing`.
const _CIRCULAR_DEEPEST = Ref{Any}(nothing)

# The LAST diagram seen, and how many calls there were. ⚠ `depth` is the
# argument of `reduce_to_circular_leave`, so it restarts at 0 whenever something
# (the rex relation derivation, for one) starts a FRESH top-level reduction —
# an overflow can happen at a shallow recorded depth despite many calls
# having been made. The last diagram is therefore the more reliable of
# the two, and both get written.
const _CIRCULAR_LAST  = Ref{Any}(nothing)
const _CIRCULAR_CALLS = Ref(0)

"""
    _circular_note_depth(fdm, depth)

Report a diagram at recursion depth `depth`; the deepest one is kept for
[`circular_capture_explosion`](@ref). Called by
[`reduce_to_circular_leave`](@ref) — not meant to be called by hand.
"""
function _circular_note_depth(fdm, depth::Int)
    CIRCULAR_EXPLOSION_ENABLED[] || return nothing
    _CIRCULAR_CALLS[] += 1
    _CIRCULAR_LAST[]   = (depth, fdm)
    d = _CIRCULAR_DEEPEST[]
    (d === nothing || depth > d[1]) && (_CIRCULAR_DEEPEST[] = (depth, fdm))
    return nothing
end

"""
    circular_reset_explosion()

Forget the deepest diagram seen so far. [`circular_capture_explosion`](@ref)
does this itself before it runs.
"""
function circular_reset_explosion()
    _CIRCULAR_DEEPEST[] = nothing
    _CIRCULAR_LAST[]    = nothing
    _CIRCULAR_CALLS[]   = 0
    return nothing
end

"""
    circular_save_explosion(fdm, reason; depth = -1, name = "explosion", extra = "")
        -> String

Write one decorated morphism to `CIRCULAR_EXPLOSION_DIR[]` and return the path
of the `.cd` file. Next to it goes a `.txt` with `reason`, `depth`, the boundary
word, the arm sequences, the switch position and `extra`.

The marking (`cut1`/`cut2`) is NOT part of the `.cd` format — it is written into
the sidecar, and it matters (GOTCHA: a guessed marking gives the opposite
result). Reload with

```julia
d = load_circulardecorated(pfad)
m = CircularMorphismGraph(d.graph, cut1, cut2)     # cut1/cut2 from the .txt
```
"""
function circular_save_explosion(fdm, reason::AbstractString;
                                 depth::Int = -1, name::AbstractString = "explosion",
                                 extra::AbstractString = "")
    dir = CIRCULAR_EXPLOSION_DIR[]
    mkpath(dir)
    stamp = Libc.strftime("%Y%m%d-%H%M%S", time())
    g     = fdm.m.graph
    base  = joinpath(dir, string(stamp, "-", name))
    # never overwrite a previous capture of the same second
    k = 0
    while isfile(base * ".cd")
        k += 1
        base = joinpath(dir, string(stamp, "-", name, "-", k))
    end
    save_circulardecorated(base * ".cd",
                           CircularDecorated(g, fdm.region_labels, fdm.outer_label))
    open(base * ".txt", "w") do io
        println(io, "reason      ", reason)
        println(io, "depth       ", depth)
        println(io, "word        ", g.word)
        println(io, "cut1        ", fdm.m.cut1)
        println(io, "cut2        ", fdm.m.cut2)
        println(io, "nodes       ", length(g.nodes), "   edges ", length(g.edges),
                    "   regions ", region_count(g), "   degree ", circular_degree(g))
        println(io, "arms        ", [collect(arms(nd)) for nd in g.nodes])
        println(io, "weight      ", circular_weight(g))
        println(io, "version     ", pkgversion(@__MODULE__))
        println(io, "switches    zamo=", CIRCULAR_ZAMO_ENABLED[],
                    " rex_fusion=", CIRCULAR_REX_FUSION_ENABLED[],
                    " mvalent=", CIRCULAR_REGIONWORD_MVALENT[],
                    " dot_policy=", CIRCULAR_DOT_POLICY[])
        isempty(extra) || println(io, "extra       ", extra)
    end
    return base * ".cd"
end

"""
    circular_capture_explosion(f; name = "explosion", verbose = true) -> f()

Run `f()`. If it throws — `StackOverflowError`, the growth guard, the
`circular_weight` assert, anything — write **two** diagrams to disk
([`circular_save_explosion`](@ref)) and then **re-throw**: the LAST one the run
reached (`…-last.cd`) and the DEEPEST one by recursion depth (`…-deepest.cd`,
only if it is a different diagram). The return value and the error behaviour are
unchanged; only files appear.

```julia
c = try
        circular_capture_explosion(name = "id_1213213") do
            reduce_to_circular_leave(fdm)
        end
    catch err
        nothing            # the diagram now lives in CIRCULAR_EXPLOSION_DIR[]
    end
```

⚠ Nothing is written if the run never got as far as one `reduce_to_circular_leave`
call (nothing was recorded), or if `CIRCULAR_EXPLOSION_ENABLED[]` is off.
"""
function circular_capture_explosion(f; name::AbstractString = "explosion",
                                    verbose::Bool = true)
    circular_reset_explosion()
    try
        return f()
    catch err
        reason = first(split(sprint(showerror, err), '\n'))
        last_ = _CIRCULAR_LAST[]
        deep  = _CIRCULAR_DEEPEST[]
        paths = String[]
        if last_ !== nothing
            push!(paths, circular_save_explosion(last_[2], reason; depth = last_[1],
                                                 name = name * "-last",
                                                 extra = "calls=$(_CIRCULAR_CALLS[])"))
            if deep !== nothing &&
               circular_canonical_key(deep[2].m.graph) != circular_canonical_key(last_[2].m.graph)
                push!(paths, circular_save_explosion(deep[2], reason; depth = deep[1],
                                                     name = name * "-deepest",
                                                     extra = "calls=$(_CIRCULAR_CALLS[])"))
            end
            verbose && @info("circular_capture_explosion: diagram(s) saved",
                             paths = paths, calls = _CIRCULAR_CALLS[],
                             maxdepth = deep === nothing ? -1 : deep[1])
        elseif verbose
            @info "circular_capture_explosion: nothing recorded (no reduce_to_circular_leave call)"
        end
        rethrow()
    end
end

"""
    circular_explosions(; dir = CIRCULAR_EXPLOSION_DIR[]) -> Vector{String}

The `.cd` files written so far, newest last.
"""
function circular_explosions(; dir::AbstractString = CIRCULAR_EXPLOSION_DIR[])
    isdir(dir) || return String[]
    return sort([joinpath(dir, f) for f in readdir(dir) if endswith(f, ".cd")])
end
