using Documenter, GLOQ

const ROOT = joinpath(@__DIR__, "..")
makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    sitename = "GLOQ.jl",
    modules = [GLOQ],
    authors = "Daniel Appelo <appeloda@msu.edu> and contributors.",
    pages = [
        "Home" => "index.md",
        "Model" => Any["Lindblad equation"=>"model.md",
                       "Experiments"=>"experiment.md" 
                      ],
        "Workflow and examples" => 
            Any["Workflow"=>"workflow.md",
                "Example 1" => "example1.md",
                "Example 2" => "example2.md",
                "Example 3" => "example3.md"],
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
