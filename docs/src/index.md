[![Build Status](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://davidamaro.github.io/GroupFunctions.jl/dev)
## Overview

`GroupFunctions.jl` computes matrix elements of irreducible representations of the unitary group `U(n)`, both numerically and symbolically. Given a unitary matrix and an irrep label, it returns the corresponding transformation matrix, with basis states represented by [Gelfand-Tsetlin patterns](tutorial/states.md).

```@repl index_symbolic
using GroupFunctions

U = su2_block_symbolic(2,1);
U
group_function([2,0], U)[1]
```

This computes the SU(2) irreps for symbolic matrices.
For symmetric irreps, `occupation_number` translates basis patterns into the corresponding Fock occupations; see the quantum-optics example for a full workflow.

## Installation

If this is your first time using Julia, please refer to the language [documentation](https://docs.julialang.org/en/v1/) and [tutorials](https://julialang.org/learning/tutorials/). Until the package is registered, please use the manual installation using git URL: open up `julia`, and within the interpreter perform the following:

```console
julia> # here you need to press the character `]`; the prompt turns blue
pkg> add https://github.com/davidamaro/GroupFunctions.jl
```


Requires Julia ≥ 1.6.

## Contact
Questions and suggestions: `david.amaroalcala@savba.sk`

## References
