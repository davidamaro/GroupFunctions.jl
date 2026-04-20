using Combinatorics
using GroupFunctions
using SymEngine
using Test

import LinearAlgebra: det, eigvals, norm, tr
import RandomMatrices: Haar

include("helpers.jl")
include("internals.jl")
include("double_coset_representatives.jl")
include("Aqua.jl")
include("mathematica_irrep_300.jl")
include("mathematica_additional_su3_irreps.jl")
include("mathematica_even_more_su3_irreps.jl")
include("basic_api.jl")
include("representation_properties.jl")
include("symbolic_reference_checks.jl")
include("applications.jl")
