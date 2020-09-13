module GroupFunctions

using Markdown
#import AbstractAlgebra: YoungTableau
using AbstractAlgebra

import SparseArrays: spzeros, SparseMatrixCSC

#greet() = print("Hello World!")
include("generic/GTPatterns.jl")
include("generic/YoungTabs.jl")

#import GTPattern, basis_states
export YoungTableau, axialdistance, encontrar_posicion
export GTPattern, basis_states, siguientepatron
export primero_lexi, StandardYoungTableaux, generar_matriz

end # module
