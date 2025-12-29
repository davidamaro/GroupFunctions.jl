push!(LOAD_PATH, "../src/")
using Documenter, GroupFunctions
makedocs(sitename="GroupFunctions documentation",
    pages = [
        "Getting started" => "index.md",
        "Tutorials" => ["Basis states" => "states.md",
			"Basic quantum optics" => "quantum_optics.md",
                        "Sum rules"=>"sum_rules.md",
                        "Immanants"=>"immanants.md",
                       ],
        "Documentation" => "documentation.md",
        "Index" => "docstrings.md"],
    format = Documenter.HTML(
                             assets = ["assets/favicon.ico"],
                             sidebar_sitename=false
    ),
        )
deploydocs(
    repo = "github.com/davidamaro/GroupFunctions.jl.git"
)
