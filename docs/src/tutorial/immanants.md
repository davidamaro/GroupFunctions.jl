```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Immanants

Zero-weight example for `λ = [2,1,0]`. Theory: [background page](../background/immanants.md).

```@repl immanants_zero
using GroupFunctions

imm21(A) = 2*A[1,1]*A[2,2]*A[3,3] -
           A[1,2]*A[2,3]*A[3,1] -
           A[1,3]*A[2,1]*A[3,2];

λ = [2, 1, 0];
basis = basis_states(λ);
zero_idx = findall(gt -> zweight(gt) == [0.0, 0.0], basis);
zero_idx

ψ1 = basis[zero_idx[1]];
ψ2 = basis[zero_idx[2]];

ψ1
ψ2

M = ComplexF64[
    1 2 3
    4 5 6
    7 8 10
];

g11 = group_function(λ, ψ1, ψ1, M)
g22 = group_function(λ, ψ2, ψ2, M)
g11 + g22
imm21(M)
g11 + g22 ≈ imm21(M)
```
