[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://davidamaro.github.io/GroupFunctions.jl/dev)
# Getting started

Julia package to compute entries of the irreducible representations of the
unitary group (D-functions or group functions).
It supports both numerical and symbolical group functions.

```@contents
```

## [Documentation](https://davidamaro.github.io/GroupFunctions.jl/dev/)
## Installation

### If Julia is already installed
Alternatively, you can install the package directly from the repository:

```console
user@machine:~$ mkdir new_code && cd new_code
user@machine:~$ julia --project=.
julia> ] add https://github.com/davidamaro/GroupFunctions.jl
```

### If Julia is not installed

- **Mac**: Use `juliaup`. Installing Julia via `brew` is not recommended.
- **Linux**: Use the appropriate package manager (e.g., `sudo pacman -S julia`).
- **Windows**: Run `winget install julia -s msstore` in your terminal and follow the steps.
## Computing a single group function (symbolical)

```julia
my_fav_irrep = [2, 1, 0]
my_fav_basis = basis_states(my_fav_irrep)
output = group_function(my_fav_irrep, my_fav_basis[1], my_fav_basis[3])
julia_to_mma(output)
```

## Contact
For questions and suggestions : `david.amaroalcala@ucalgary.ca`
