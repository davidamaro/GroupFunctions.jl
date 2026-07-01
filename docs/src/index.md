```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

## Overview

`GroupFunctions.jl` computes matrix elements of irreducible representations of the unitary group `U(n)`. The code accepts as input both numerical and symbolic matrices. Given a unitary matrix and an irrep label, the package returns the corresponding transformation matrix in a basis of [Gelfand-Tsetlin patterns](tutorial/states.md).

For symbolic computations, the same API accepts symbolic block matrices such as `su2_block_symbolic`:

```jldoctest
julia> U = su2_block_symbolic(2,1);

julia> U
2×2 Matrix{SymEngine.Basic}:
 v_1_1  v_1_2
 v_2_1  v_2_2

julia> group_function([2,0], U)[1]
3×3 Matrix{SymEngine.Basic}:
             v_2_2^2        sqrt(2)*v_2_2*v_2_1              v_2_1^2
 sqrt(2)*v_1_2*v_2_2  v_1_1*v_2_2 + v_1_2*v_2_1  sqrt(2)*v_1_1*v_2_1
             v_1_2^2        sqrt(2)*v_1_1*v_1_2              v_1_1^2
```

The API also allows to compute `SU(3)` irreps from a symbolic `3×3` matrix:

```jldoctest
julia> U = [GroupFunctions.SymEngine.symbols("u_$(i)_$(j)") for i in 1:3, j in 1:3];

julia> U
3×3 Matrix{SymEngine.Basic}:
 u_1_1  u_1_2  u_1_3
 u_2_1  u_2_2  u_2_3
 u_3_1  u_3_2  u_3_3

julia> group_function([2,0,0], U)[1]
6×6 Matrix{SymEngine.Basic}:
             u_3_3^2        sqrt(2)*u_3_2*u_3_3  …              u_3_1^2
 sqrt(2)*u_2_3*u_3_3  u_2_2*u_3_3 + u_2_3*u_3_2     sqrt(2)*u_2_1*u_3_1
 sqrt(2)*u_1_3*u_3_3  u_1_2*u_3_3 + u_1_3*u_3_2     sqrt(2)*u_1_1*u_3_1
             u_2_3^2        sqrt(2)*u_2_2*u_2_3                 u_2_1^2
 sqrt(2)*u_2_3*u_1_3  u_2_2*u_1_3 + u_2_3*u_1_2     sqrt(2)*u_2_1*u_1_1
             u_1_3^2        sqrt(2)*u_1_2*u_1_3  …              u_1_1^2
```

## Installation

If you are new to Julia, see the language [documentation](https://docs.julialang.org/en/v1/) and [tutorials](https://julialang.org/learning/tutorials/). Then install the package from the Julia registry:

```console
julia> ]
pkg> add GroupFunctions
```


Requires Julia ≥ 1.6.

## Contact

Send questions and suggestions to `david.amaroalcala@savba.sk`.

## References
