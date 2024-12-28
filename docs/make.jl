push!(LOAD_PATH, "../src/")
using Documenter, GroupFunctions
makedocs(sitename="GroupFunctions documentation",
    pages = [
        "Getting started" => "index.md",
        "Tutorials" => ["Basis states" => "states.md",
                        "Sum rules"=>"sum_rules.md"],
        "Documentation" => "documentation.md",
        "Index" => "docstrings.md"],
    format = Documenter.HTML(
                             assets = ["assets/favicon.ico"]
    ),
        )
deploydocs(
    repo = "github.com/davidamaro/GroupFunctions.jl.git"
)
