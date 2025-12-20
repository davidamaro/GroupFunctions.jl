[![Build Status](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://davidamaro.github.io/GroupFunctions.jl/dev)
# Getting started

Julia package to compute entries of the irreducible representations of the
unitary group (D-functions or group functions).
It supports both numerical and symbolic group functions.

## [Documentation](https://davidamaro.github.io/GroupFunctions.jl/dev/)
## Installation

### Julia is already installed
You can install the package directly from the repository:

```console
user@machine:~$ mkdir new_code 
user@machine:~$ cd new_code
user@machine:~$ julia --project=.
julia> # here you need to press the character `]`; the prompt turns blue
pkg > add https://github.com/davidamaro/GroupFunctions.jl
```

### Getting Julia

- **Mac**: Use `juliaup`. Installing Julia via `brew` is not recommended.
- **Linux**: Use the appropriate package manager (e.g., `sudo pacman -S julia`).
- **Windows**: Run `winget install julia -s msstore` in your terminal and follow the steps.
## Computing a single group function (symbolical)

```julia
irrep = [2, 1, 0]
basis = basis_states(irrep)
output = group_function(irrep, basis[1], basis[3])
julia_to_mma(output)
```

## Contact
For questions and suggestions: `david.amaroalcala@ucalgary.ca`
