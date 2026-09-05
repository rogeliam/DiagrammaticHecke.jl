using DiagrammaticHecke
using Test

# Partial runs: select groups via the environment variable DH_TESTS
# (comma-separated), e.g.
#
#   DH_TESTS=words julia --project=. test/runtests.jl
#   DH_TESTS=circular,zamo julia --project=. test/runtests.jl
#
# Without DH_TESTS everything runs. `--list` prints the group names.
const TEST_GROUPS = Dict(
    "words"    => ["words.jl"],
    "algebra"  => ["algebra.jl"],
    "diagram"  => ["diagram.jl"],
    "faces"    => ["faces.jl"],
    "morphism" => ["morphism.jl"],
    "circular" => ["fixtures_circular.jl", "circulargraph.jl", "circularrules.jl",
                   "circularbraidrules.jl", "circularmisc.jl", "circularweight.jl",
                   "circulardriver.jl", "circularlightleaves.jl",
                   "circulardecoratedrules.jl", "circular2parallel.jl",
                   "circular2parallelapply.jl",
                   "circularmvalentwords.jl", "circulardotmerge.jl",
                   "circularexplosion.jl",                    "circularbraidchannels.jl",                    "circulargen12merge.jl", "circulardotongen12.jl",
                   "circulargen12expand.jl", "circularbraidongen12.jl",
                   "circulargen12ongen12.jl", "circulargen12twoedges.jl",
                   "circulargenericpreimage.jl", "circularbraidedges.jl",
                   "circularregionrules.jl", "circulargraphio.jl",
                   "circulardotslide.jl"],
    "decorated" => ["decorated.jl"],
    "pairing"   => ["pairing.jl"],
    "dlbasis"   => ["dlbasis.jl"],
    "circularpairing" => ["circularpairing.jl"],
    "clbasis"   => ["clbasis.jl"],
    "circularcomponents" => ["circularcomponents.jl"],
    # ⚠️ the first zamo run derives Z2 (~40 s, once per process).
    "zamo"      => ["zamotermrules.jl", "circularzamoregion.jl"],
)

if "--list" in ARGS
    foreach(println, sort(collect(keys(TEST_GROUPS))))
    exit(0)
end

const SELECTED = let raw = get(ENV, "DH_TESTS", "")
    isempty(strip(raw)) ? nothing : Set(strip.(split(raw, ",")))
end

if SELECTED !== nothing
    unknown = setdiff(SELECTED, keys(TEST_GROUPS))
    isempty(unknown) || error("DH_TESTS: unknown group(s) " *
                              join(sort(collect(unknown)), ", ") * ". Known: " *
                              join(sort(collect(keys(TEST_GROUPS))), ", "))
    println("DH_TESTS active — groups only: ", join(sort(collect(SELECTED)), ", "))
end

# `true` if the file `name` is due in the current run (its group is selected).
const _FILE_GROUP = Dict(f => g for (g, fs) in TEST_GROUPS for f in fs)
runs(name) = SELECTED === nothing || get(_FILE_GROUP, name, name) in SELECTED ||
             name in SELECTED

@testset "DiagrammaticHecke" begin
    runs("words.jl") && include("words.jl")
    runs("algebra.jl") && include("algebra.jl")
    runs("diagram.jl") && include("diagram.jl")
    runs("faces.jl") && include("faces.jl")
    runs("morphism.jl") && include("morphism.jl")

    # shared planar fixtures: load before the circular files.
    runs("fixtures_circular.jl") && include("fixtures_circular.jl")
    for f in TEST_GROUPS["circular"][2:end]
        runs(f) && include(f)
    end

    runs("decorated.jl") && include("decorated.jl")
    runs("pairing.jl") && include("pairing.jl")
    runs("dlbasis.jl") && include("dlbasis.jl")
    runs("circularpairing.jl") && include("circularpairing.jl")
    runs("clbasis.jl") && include("clbasis.jl")
    runs("circularcomponents.jl") && include("circularcomponents.jl")
    runs("zamotermrules.jl") && include("zamotermrules.jl")
    runs("circularzamoregion.jl") && include("circularzamoregion.jl")
end
