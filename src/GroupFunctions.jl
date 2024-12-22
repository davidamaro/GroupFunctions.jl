module GroupFunctions
  using Markdown
  using AbstractAlgebra
  using SymEngine
  import AbstractAlgebra: YoungTableau
  import SparseArrays: spzeros, SparseMatrixCSC
  import SymEngine: expand 
  include("generic/GTPatterns.jl")
  include("generic/YoungTabs.jl")
  include("generic/DF.jl")
  include("generic/PermGroups.jl")
  include("generic/Misc.jl")
  export YoungTableau, axialdistance, encontrar_posicion
  export GTPattern, basis_states, siguientepatron, siguientepatron!
  export primero_lexi, StandardYoungTableaux, generar_matriz
  export indice_tablon_semistandard, content, Θ
  export group_function, @mma_str, mma_to_julia
  export zweight, pweight
  export gtinicial
  export generar_matriz
  export bloquesun, simplefactorization, simple
  export expand
  export julia_to_mma
end # module
