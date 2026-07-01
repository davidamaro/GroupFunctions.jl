using Documenter
using GroupFunctions

const run_doctests = !("--nodoctests" in ARGS)

makedocs(
    sitename = "GroupFunctions.jl",
    modules = [GroupFunctions],
    clean = true,
    doctest = run_doctests,
    checkdocs = :exports,
    pages = [
        "Getting started" => "index.md",
        "Background" => [
            "Basis states" => "background/states.md",
            "Computation of group functions" => "background/group_functions.md",
            "Characters" => "background/characters.md",
            "Immanants" => "background/immanants.md",
        ],
        "Tutorial" => [
            "Basis states" => "tutorial/states.md",
            "Computation of group functions" => "tutorial/group_functions.md",
            "Characters" => "tutorial/characters.md",
            "Immanants" => "tutorial/immanants.md",
        ],
        "Applications" => [
            "HOM effect" => "applications/quantum_optics.md",
            "Qubit transmission" => "applications/qubit_transmission.md",
            "Sum rules" => "applications/sum_rules.md",
            "Coincidence computations" => "applications/coincidences.md",
        ],
        "API reference" => "documentation.md",
    ],
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        sidebar_sitename = false,
        size_threshold_warn = 204800,
        size_threshold = 409600,
    ),
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/davidamaro/GroupFunctions.jl.git",
        target = "build",
    )
end
