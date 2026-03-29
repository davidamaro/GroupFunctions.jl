# Characters
SU(d) characters are traces of irreducible representation matrices. In
`GroupFunctions.jl`, you can compute the matrix with `group_function(irrep, U)`
and then take its trace.

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

## Example: variance of sampled character values

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det, tr
using Statistics: var

nsamples = 1_000
n = 3
irrep = [2, 1, 0]

values = Vector{ComplexF64}(undef, nsamples)
for i in 1:nsamples
    U = rand(Haar(2), n)
    U_su = U / det(U)^(1 / n)
    rep, _states = group_function(irrep, U_su)
    values[i] = tr(rep)
end

variance = var(values)
variance
```
