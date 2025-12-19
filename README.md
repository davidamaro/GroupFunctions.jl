[![CI](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/davidamaro/GroupFunctions.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Docs – stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://davidamaro.github.io/GroupFunctions.jl/stable)
[![Docs – dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://davidamaro.github.io/GroupFunctions.jl/dev)
# GroupFunctions.jl

A Julia library to compute D-functions, which are entries of the irreducible representations of the unitary group U(d). These entries can be numeric or symbolic.

## Highlights
- Numerical and symbolic D-functions for U(d) irreps
- Gelfand–Tsetlin basis construction via `basis_states`
- Works with Haar-random unitaries (e.g. using RandomMatrices.jl)
- Export symbolic results to Mathematica with `julia_to_mma`

## Installation

Requires Julia 1.6 or newer.

- **From the Julia registry (recommended):**

```julia
pkg> add GroupFunctions
```

- **From this repository:**

```console
user@machine:~$ mkdir new_code && cd new_code
user@machine:~$ julia --project=.
julia> ] add https://github.com/davidamaro/GroupFunctions.jl
```

### Installing Julia

- **Mac**: Use `juliaup`. Installing Julia via `brew` is not recommended.
- **Linux**: Use the appropriate package manager (e.g., `sudo pacman -S julia`).
- **Windows**: Run `winget install julia -s msstore` in your terminal and follow the steps.

## Quick start

```julia
using GroupFunctions
using RandomMatrices  # Optional: for Haar-random unitaries

irrep = [2, 1, 0]
U = rand(Haar(2), 3)
basis = basis_states(irrep)

# Numerical D-function entry
group_function(irrep, basis[1], basis[3], U)

# Symbolic D-function entry
sym = group_function(irrep, basis[1], basis[3])
julia_to_mma(sym)  # Mathematica-friendly expression
```

For more examples and API details, see the documentation: https://davidamaro.github.io/GroupFunctions.jl/dev/

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. A to-do list is included in the `todo.txt` file.

- Run tests: `julia --project -e 'using Pkg; Pkg.test()'`
- Build docs locally: `julia --project=docs docs/make.jl`

## License

GroupFunctions.jl is distributed under the MIT License (see `LICENSE`). While the package still contains code derived from AbstractAlgebra.jl to handle Young tableaux, it follows the same (MIT) license terms.

## References

1. [J Grabmeier and A Kerber, "The evaluation of irreducible polynomial representations of the general linear groups and of the unitary groups over fields of characteristic 0" Acta Appl. Math, 1987](http://dx.doi.org/10.1007/BF00046717)
2. [A Alex et al, "A numerical algorithm for the explicit calculation of SU(N) and SL(N, C) Clebsch–Gordan coefficients" J. Math. Phys. 2011 ](http://dx.doi.org/10.1063/1.3521562)
3. [D Amaro-Alcala et al "Sum rules in multiphoton coincidence rates" Phys. Lett. A 2020](http://dx.doi.org/10.1016/j.physleta.2020.126459)
4. [AbstractAlgebra.jl](https://nemocas.github.io/AbstractAlgebra.jl/stable/)

## Citation

Citation information is pending. In the meantime, please cite the repository URL and version tag (e.g., “GroupFunctions.jl v0.1.5, 2024, https://github.com/davidamaro/GroupFunctions.jl”).
