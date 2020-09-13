module GroupFunctions

using Markdown
#import AbstractAlgebra: YoungTableau
using AbstractAlgebra
using SymEngine

import SparseArrays: spzeros, SparseMatrixCSC

#greet() = print("Hello World!")
include("generic/GTPatterns.jl")
include("generic/YoungTabs.jl")
include("generic/DF.jl")
include("generic/PermGroups.jl")

#import GTPattern, basis_states
export YoungTableau, axialdistance, encontrar_posicion
export GTPattern, basis_states, siguientepatron
export primero_lexi, StandardYoungTableaux, generar_matriz
export indice_tablon_semistandard, content, Θ
export group_function

end # module
