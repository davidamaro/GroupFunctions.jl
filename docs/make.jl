push!(LOAD_PATH, "../src/")
using Documenter, GroupFunctions
makedocs(sitename="GroupFunctions documentation",
    pages = [
        "Getting started" => "index.md",
        "Tutorial" => [
            "Basis states" => "tutorial/states.md",
            "Characters" => "tutorial/characters.md",
            "Calculation of group functions" => "tutorial/group_functions.md",
            "Immanants" => "tutorial/immanants.md",
        ],
        "Applications" => [
            "Example: flavor mixing" => "applications/flavor.md",
            "Example: HOM effect" => "applications/quantum_optics.md",
            "Example: qubit transmission" => "applications/qubit_transmission.md",
            "Sum rules" => "applications/sum_rules.md",
        ],
        "Background" => [
            "Basis states" => "background/states.md",
            "Characters" => "background/characters.md",
            "Calculation of group functions" => "background/group_functions.md",
            "Immanants" => "background/immanants.md",
        ],
        "Documentation" => "documentation.md",
    ],
    format = Documenter.HTML(
                             assets = ["assets/favicon.ico"],
                             sidebar_sitename=false
    ),
        )
deploydocs(
    repo = "github.com/davidamaro/GroupFunctions.jl.git"
)
