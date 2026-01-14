[![Build Status](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://davidamaro.github.io/GroupFunctions.jl/dev)
## Overview

A Julia package for computing matrix elements of irreducible representations of the unitary group U(n), commonly known as *D-functions* or *group functions*. Supports both numerical evaluation and symbolic computation. Potential applications include quantum information and quantum optics problems (boson sampling, photonic state preparation), as well as general quantum many-body physics. 

The basic functionality of the function is the following: given a (symbolic or numerical) U(n) matrix, calculate the transformation matrix in a given representation (e.g. fixed particle number subspace), mathematically described via basis states involving Gelfand-Tseltsin patterns.

As a quick example, the following code calculates the probability of photonic $\ket{1,1}$ state transforming to $\ket{2,0}$ photons after going through a 50:50 beamsplitter.
```julia
using GroupFunctions

λ = [2, 0]                           # 2 photons, 2 modes
basis = basis_states(λ)              # enumerate GT patterns

# Check which Fock state each pattern represents
for b in basis
    println(pweight(b), " → ", b)    # pweight gives occupation numbers
end
# Output: [2,0], [1,1], [0,2]

# Beamsplitter unitary and transition amplitude |1,1⟩ → |2,0⟩
BS = su2_block(2, 1, (0., π/2, 0.))
initial = basis[findfirst(b -> pweight(b) == [1,1], basis)]
final   = basis[findfirst(b -> pweight(b) == [2,0], basis)]

amp = group_function(λ, final, initial, BS)
println("Probability: ", abs2(amp))  # ≈ 0.5 (HOM effect)
```

## Installation

If this is your first time using Julia, please refer to the language [documentation](https://docs.julialang.org/en/v1/) and [tutorials](https://julialang.org/learning/tutorials/). Until the package is registered, please use the manual installation using git URL: open up `julia`, and within the interpreter perform the following:

```console
julia> # here you need to press the character `]`; the prompt turns blue
pkg> add https://github.com/davidamaro/GroupFunctions.jl
```


Requires Julia ≥ 1.6.

## Contact
Questions and suggestions: `david.amaroalcala@ucalgary.ca`

## References
