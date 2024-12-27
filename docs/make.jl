push!(LOAD_PATH, "../src/")
using Documenter, GroupFunctions
makedocs(sitename="GroupFunctions documentation")
deploydocs(
    repo = "github.com/davidamaro/GroupFunctions.jl.git",
)
