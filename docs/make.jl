using Documenter, GLOQ

const ROOT = joinpath(@__DIR__, "..")
makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    sitename = "GLOQ.jl",
    modules = [GLOQ],
    authors = "Daniel Appelo <appeloda@msu.edu> and contributors.",
    pages = Any[
        "Home" => "index.md",
        "Workflow" => "workflow.md",
        "Types" => "types.md",
        "Methods" => "methods.md",
        "Index" => "function-index.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
 deploydocs(
    repo = "github.com/appelo/GLOQ.jl.git",
    devbranch = "main"
)
