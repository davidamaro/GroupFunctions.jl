@doc raw"""
GroupFunctions.jl computes symbolic and numeric matrix elements of irreducible
representations of `U(d)` and `SU(d)` using Gelfand-Tsetlin patterns and Young
tableaux.

Repository: <https://github.com/davidamaro/GroupFunctions.jl>

Documentation: <https://davidamaro.github.io/GroupFunctions.jl>
"""
module GroupFunctions

using SymEngine

import SparseArrays: SparseMatrixCSC, spzeros
import SymEngine: expand

include("exports.jl")

@doc raw"""
    expand(expr)

Re-export `SymEngine.expand` through `GroupFunctions` so symbolic group
functions can be expanded without importing `SymEngine` explicitly.
"""
expand

################################################################################
#
# Combinatorial data structures
#
################################################################################

include("generic/GTPatterns.jl")
include("generic/Tableaux.jl")
include("generic/PermGroups.jl")
include("generic/YoungTabs.jl")

################################################################################
#
# Group-function algorithms
#
################################################################################

include("generic/DF.jl")

################################################################################
#
# Matrix constructors and conversion utilities
#
################################################################################

include("generic/Misc.jl")

"""
    doctestsetup()

Return the setup expression shared by the package documentation doctests.
"""
function doctestsetup()
    return :(using GroupFunctions)
end

end # module
