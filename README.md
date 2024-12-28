[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://davidamaro.github.io/GroupFunctions.jl/dev)
# GroupFunctions.jl

A Julia library to compute D-functions, which are entries of the irreducible representations of the unitary group U(d). These entries can be numeric or symbolic.

## Installation

### Julia is already installed
Alternatively, you can install the package directly from the repository:

```console
user@machine:~$ mkdir new_code && cd new_code
user@machine:~$ julia --project=.
julia> ] add https://github.com/davidamaro/GroupFunctions.jl
```

### Installing Julia

- **Mac**: Use `juliaup`. Installing Julia via `brew` is not recommended.
- **Linux**: Use the appropriate package manager (e.g., `sudo pacman -S julia`).
- **Windows**: Run `winget install julia -s msstore` in your terminal and follow the steps.

## Usage

```julia
julia> using RandomMatrices # You may be asked to install it. Just answer yes.
julia> using GroupFunctions
julia> my_fav_irrep = [2, 1, 0]
julia> my_fav_matrix = rand(Haar(2), 3)
julia> my_fav_basis = basis_states(my_fav_irrep)
julia> # Check out some cool Gelfand-Tsetlin patterns:
julia> my_fav_basis[1]
julia> my_fav_basis[3]
julia> group_function(my_fav_irrep, my_fav_basis[1], my_fav_basis[3], my_fav_matrix)
julia> # For a symbolic matrix, simply omit the matrix argument
julia> output = group_function(my_fav_irrep, my_fav_basis[1], my_fav_basis[3])
julia> # Translate the symbolic D-function for use in Mathematica
julia> julia_to_mma(output)
```

For more examples, see the "Tutorials" section in the [documentation](https://davidamaro.github.io/GroupFunctions.jl/dev/).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. A to-do list is included in the `todo.txt` file.

## License

Until the code from AbstractAlgebra.jl (to deal with Young tableaux) is removed from this package,
the license will align with AbstractAlgebra.jl's.

## References

1. [J Grabmeier and A Kerber, "The evaluation of irreducible polynomial representations of the general linear groups and of the unitary groups over fields of characteristic 0" Acta Appl. Math, 1987](http://dx.doi.org/10.1007/BF00046717)
2. [A Alex et al, "A numerical algorithm for the explicit calculation of SU(N) and SL(N, C) Clebsch–Gordan coefficients" J. Math. Phys. 2011 ](http://dx.doi.org/10.1063/1.3521562)
3. [D Amaro-Alcala et al "Sum rules in multiphoton coincidence rates" Phys. Lett. A 2020](http://dx.doi.org/10.1016/j.physleta.2020.126459)
4. [AbstractAlgebra.jl](https://nemocas.github.io/AbstractAlgebra.jl/stable/)

## Citation

Pending.

