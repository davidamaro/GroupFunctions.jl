# Documentation

## Functions

```@docs
    group_function
    character_weyl
```

### Example: SU(3) character check

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det, eigvals, tr

irrep = [2, 1, 0]
U = rand(Haar(2), 3)  # 3×3 Haar-random unitary matrix
U_su = U / det(U)^(1/3)

rep, _states = group_function(irrep, U_su)
lambdas = eigvals(U_su)
ideal_character = lambdas[1]/lambdas[2] + lambdas[2]/lambdas[1] +
                  lambdas[1]/lambdas[3] + lambdas[3]/lambdas[1] +
                  lambdas[2]/lambdas[3] + lambdas[3]/lambdas[2] + 2
tr(rep) ≈ ideal_character
```

## Gelfand-Tsetlin patterns

GT patterns are used to denote the basis states.

```@docs
    GTPattern
```

## Construction of matrices

SU(2) blocks are used to construct unitary matrices. For SU(d) characters, use
`character_weyl` to evaluate the Weyl/Schur determinant formula.

```@docs
    su2_block
```
