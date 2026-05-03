```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Characters
SU(d) characters are traces of irreducible representation matrices. In
`GroupFunctions.jl`, you can compute them directly with `character(irrep, U)`.

## Example: evaluate a sampled character

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det

irrep = [2, 1, 0]
U = rand(Haar(2), 3)  # 3×3 Haar-random unitary matrix
U_su = U / det(U)^(1 / 3)

χ = character(irrep, U_su)
χ
```

## Example: variance of sampled character values

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
    values[i] = character(irrep, U_su)
end

variance = var(values)
variance
```

## Example: Schur polynomial

`schur_polynomial` computes the same character after restricting the matrix to
the diagonal torus. The irrep label uses the same partition convention as the
rest of the package: for `SU(3)`, pass a length-three vector such as `[2, 1, 0]`.

```julia
using GroupFunctions

schur_polynomial([2, 1, 0])
schur_polynomial([2, 1, 0], [2, 3, 5])
```
