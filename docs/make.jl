# docs/make.jl — build the manual with Documenter.
#
#     julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#     julia --project=docs docs/make.jl
#
# The API pages are `@autodocs` blocks filtered by source file, so the manual
# follows the same layer order as `src/DiagrammaticHecke.jl` and no exported
# name has to be listed twice.
using Documenter, DiagrammaticHecke

# The repository has no `origin` until it is pushed, and Documenter refuses to
# guess one. `DOCS_REPO` overrides the placeholder without editing this file.
const REPO = get(ENV, "DOCS_REPO", "github.com/<user>/DiagrammaticHecke.jl.git")

makedocs(
    sitename = "DiagrammaticHecke.jl",
    modules  = [DiagrammaticHecke],
    authors  = "Liam Rogel",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://<user>.github.io/DiagrammaticHecke.jl",
        edit_link  = "main",
    ),
    repo = REPO,
    pages = [
        "Home"                => "index.md",
        "Guide" => [
            "Words and rules"      => "words.md",
            "The Soergel ring"     => "ring.md",
            "Diagrams"             => "diagrams.md",
            "Morphisms and leaves" => "morphisms.md",
            "Circular diagrams"    => "circular.md",
            "The rewrite rules"    => "rules.md",
            "Every rule, drawn"    => "rules-gallery.md",
            "The two bases"        => "bases.md",
            "Zamolodchikov"        => "zamolodchikov.md",
            "Rendering"            => "rendering.md",
            "File formats"         => "formats.md",
        ],
        "Pictures" => [
            "A gallery" => "gallery.md",
        ],
        "API reference" => [
            "Words"     => "api/words.md",
            "Algebra"   => "api/algebra.md",
            "Diagrams"  => "api/diagram.md",
            "Morphisms" => "api/morphism.md",
            "Circular"  => "api/circular.md",
            "Rules: core"        => "api/rules-core.md",
            "Rules: braid nodes" => "api/rules-gen12.md",
            "Rules: relations"   => "api/rules-relations.md",
            "Rendering" => "api/render.md",
        ],
    ],
    # Not every docstring is cross-referenced yet, and the guide links to files
    # in the repository rather than to pages. Failing the build on that would
    # stop the manual over a formality.
    warnonly = [:cross_references, :missing_docs, :docs_block],
)

deploydocs(repo = REPO, devbranch = "main")
