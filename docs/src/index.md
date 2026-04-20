```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

## Overview

`GroupFunctions.jl` computes matrix elements of irreducible representations of the unitary group `U(n)`, both numerically and symbolically. Given a unitary matrix and an irrep label, it returns the corresponding transformation matrix, with basis states represented by [Gelfand-Tsetlin patterns](tutorial/states.md).

For symbolic calculations, the same API works with symbolic block matrices such as `su2_block_symbolic`:

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

For symmetric irreps, `occupation_number` translates basis patterns into the corresponding Fock occupations; see the quantum-optics example for a full workflow.

## Installation

If this is your first time using Julia, please refer to the language [documentation](https://docs.julialang.org/en/v1/) and [tutorials](https://julialang.org/learning/tutorials/). Install the package from the Julia registry with:

```console
julia> ]
pkg> add GroupFunctions
```


Requires Julia ≥ 1.6.

## Contact
Questions and suggestions: `david.amaroalcala@savba.sk`

## References
