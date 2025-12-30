# Characters

SU(d) characters can be evaluated with the Weyl/Schur determinant formula.

```@docs
    character_weyl
```

## Example: SU(3) character check

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det, eigvals, tr

irrep = [2, 1, 0]
U = rand(Haar(2), 3)  # 3×3 Haar-random unitary matrix
U_su = U / det(U)^(1 / 3)

rep, _states = group_function(irrep, U_su)
lambdas = eigvals(U_su)
ideal_character = lambdas[1] / lambdas[2] + lambdas[2] / lambdas[1] +
                  lambdas[1] / lambdas[3] + lambdas[3] / lambdas[1] +
                  lambdas[2] / lambdas[3] + lambdas[3] / lambdas[2] + 2
tr(rep) ≈ ideal_character
```

## Example: variance of Weyl character samples

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det
using Statistics: var

nsamples = 1_000
n = 3
irrep = [2, 1, 0]

values = Vector{ComplexF64}(undef, nsamples)
for i in 1:nsamples
    U = rand(Haar(2), n)
    U_su = U / det(U)^(1 / n)
    values[i] = character_weyl(irrep, U_su)
end

variance = var(values)
# Expect variance ≈ 1 for these samples.
variance
```
