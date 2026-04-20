```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

## Overview

`GroupFunctions.jl` computes matrix elements of irreducible representations of the unitary group `U(n)`, both numerically and symbolically. Given a unitary matrix and an irrep label, it returns the corresponding transformation matrix, with basis states represented by [Gelfand-Tsetlin patterns](tutorial/states.md).

```jldoctest
julia> values, basis = group_function([2, 0], ComplexF64[1 0; 0 1]);

julia> size(values)
(3, 3)

julia> length(basis)
3

julia> iszero(character([2, 1, 0], ComplexF64[1 0 0; 0 -1 0; 0 0 -1]))
true
```

For symbolic calculations, the same API works with symbolic block matrices such as `su2_block_symbolic`. For symmetric irreps, `occupation_number` translates basis patterns into the corresponding Fock occupations; see the quantum-optics example for a full workflow.

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
