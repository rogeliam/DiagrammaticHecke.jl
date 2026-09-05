# render/CellLabels.jl — turns cell/region data (number, polynomial, distance) into
# label text and SVG text elements, using the anchors from Cells.jl. Included
# after Cells.jl.


"""
    _region_label_texts(regs, region_labels, distances; debug, region_based=true)
        -> Vector{String}

Like [`_cell_label_texts`](@ref), but per REGION: number, polynomial and distance
are all three properties of THE SAME region, all indexed by `R.id`.

**The labels hang on regions** (`CircularDecorated.region_labels`), so the polynomial
stands where it is stored, with no selection rule.

`region_based = false` only means that `distances` is face-indexed (the plain call
path via the meanwhile removed `circular_cell_distances`; see `circular_region_distances` in
circular/CircularRegion.jl for the current way); the LABELS are always region-indexed.
"""
function _region_label_texts(regs::Vector{Region}, region_labels,
                             distances::Union{Nothing,Vector{Int}};
                             debug::Bool = distances !== nothing,
                             region_based::Bool = true)
    poly(i) = isone(region_labels[i]) ? "" : soergel_str(region_labels[i])
    debug || return [poly(i) for i in eachindex(regs)]
    return map(enumerate(regs)) do (i, R)
        parts = [string(R.id)]
        p = poly(i)
        isempty(p) || push!(parts, p)
        if distances !== nothing
            d = region_based ? distances[i] : distances[R.cell]
            push!(parts, d < 0 ? "·" : string(d))
        end
        length(parts) == 1 ? parts[1] : string("(", join(parts, ", "), ")")
    end
end

function _polynomial_label_svg(txt::String)
    # Render caret exponents (e.g. "α₂^2") as SVG <tspan> superscripts.
    # Use character-array iteration to avoid all UTF-8 byte-index pitfalls.
    out = IOBuffer()
    chars = collect(txt)
    i = 1
    while i <= length(chars)
        c = chars[i]
        if c != '^'
            print(out, c)
            i += 1
            continue
        end
        i += 1
        start = i
        while i <= length(chars) && isdigit(chars[i])
            i += 1
        end
        exp = String(chars[start:i-1])
        print(out, "<tspan baseline-shift=\"super\" font-size=\"0.7em\">", exp, "</tspan>")
    end
    return String(take!(out))
end

"""
    _cell_label_texts(face_labels, distances) -> Vector{String}

The text written into every interior cell.

NORMAL (`distances === nothing`, no debug mode): the polynomial, with the trivial
label `1` hidden as an empty string, as `decorated_svg`/`circular_decorated_svg` do.

DEBUG (`debug_labels()`/`debug_circular()`): **`(number, polynomial)`**, something
like `(1, α₁)`. If the polynomial is `1`, just the number. That way both are visible
and usable.

So `(3, α₂)` for a cell with a nontrivial label and plain `3` when the label is
trivial. The TRIVIAL label here is `1` (the unit), not `0` — a zero polynomial
would never occur as a cell label, it would make the whole diagram 0.

If a DISTANCE is also known (`distances !== nothing`, the path on which
`find_d4_match` decides), it is appended as a third field: `(3, α₂, 1)` resp.
`(3, ·)` for an unreachable cell (`-1`, as in `cell_distance_svg`). Without
distances it stays a pair.

The CELL NUMBER is the index in `inner_faces`/`face_labels` — the same number
`cell_number_svg` draws and the one `apply_circular_d4`/`find_circular_d4_match` name in
their error messages.
"""
function _cell_label_texts(face_labels, distances::Union{Nothing,Vector{Int}};
                           debug::Bool = distances !== nothing)
    debug || return [isone(f) ? "" : soergel_str(f) for f in face_labels]
    return map(enumerate(face_labels)) do (c, f)
        parts = [string(c)]
        isone(f) || push!(parts, soergel_str(f))
        if distances !== nothing
            d = distances[c]
            push!(parts, d < 0 ? "·" : string(d))
        end
        length(parts) == 1 ? parts[1] : string("(", join(parts, ", "), ")")
    end
end

"""
    _region_anchors(g::WordGraph; mark_between = nothing) -> Dict{Int,Tuple}

Anchors per region for a `WordGraph` — collect the pieces for
[`_region_anchors`](@ref).
"""
function _region_anchors(g::WordGraph; mark_between::MarkSpec = nothing)
    walks, cell_of_face, _, _ = _face_walks(g)
    leafxy, nodexy = circular_tutte_positions(circular(g); mark_between = mark_between)
    return _region_anchors(regions(g), walks, cell_of_face, leafxy, nodexy,
                           length(g.word))
end

function _cell_labels_svg(g::WordGraph, labels; size::Int, r::Real,
                            mark_between::MarkSpec = nothing)
    anchors = _region_anchors(g; mark_between = mark_between)
    cx = cy = size / 2
    io = IOBuffer()
    for c in eachindex(labels)
        haskey(anchors, c) || continue          # free circles etc. have no anchor
        txt = string(labels[c])
        isempty(txt) && continue                # hide the label 1
        (x, y) = anchors[c]
        X = cx + r * x; Y = cy + r * y
        svg_txt = _polynomial_label_svg(txt)
        print(io, """<text x="$X" y="$Y" text-anchor="middle" dy="0.35em" """ *
                  """font-size="17" font-style="italic" fill="#555" stroke="white" """ *
                  """stroke-width="4" paint-order="stroke">$svg_txt</text>""")
    end
    return String(take!(io))
end

