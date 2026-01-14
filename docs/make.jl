push!(LOAD_PATH, "../src/")
using Documenter, GroupFunctions
makedocs(sitename="GroupFunctions documentation",
    pages = [
        "Getting started" => "index.md",
        "Tutorials" => ["Basis states" => "states.md",
			"Calculation of group functions" => "group_functions.md",
			"Example: HOM effect" => "quantum_optics.md",
			"Example: qubit transmission" => "qubit_transmission.md",
                        "Characters" => "characters.md",
                        "Sum rules"=>"sum_rules.md",
                        "Immanants"=>"immanants.md",
			"Example: flavor mixing" => "flavor.md"
                       ],
        "Documentation" => "documentation.md"],
    format = Documenter.HTML(
                             assets = ["assets/favicon.ico"],
                             sidebar_sitename=false
    ),
        )
deploydocs(
    repo = "github.com/davidamaro/GroupFunctions.jl.git"
)
