```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Characters
SU(d) characters can be evaluated directly with `character(irrep, U)`.

```jldoctest
julia> U = ComplexF64[1 0 0; 0 -1 0; 0 0 -1];

julia> iszero(character([2, 1, 0], U))
true
```

## Example: SU(3) character check
Here we present an example with a known formula for the character in the irrep
$\lambda = (2, 1, 0)$, and compare it against `character`.

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det, eigvals

irrep = [2, 1, 0]
U = rand(Haar(2), 3)  # 3×3 Haar-random unitary matrix
U_su = U / det(U)^(1 / 3)

lambdas = eigvals(U_su)
ideal_character = lambdas[1] / lambdas[2] + lambdas[2] / lambdas[1] +
                  lambdas[1] / lambdas[3] + lambdas[3] / lambdas[1] +
                  lambdas[2] / lambdas[3] + lambdas[3] / lambdas[2] + 2
character(irrep, U_su) ≈ ideal_character
```
