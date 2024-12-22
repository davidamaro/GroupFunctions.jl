# GroupFunctions.jl

A Julia library to compute D-functions, which are entries of the irreducible representations of the unitary group U(d). These entries can be numeric or symbolic.

## Installation

### If Julia is already installed
Alternatively, you can install the package directly from the repository:

```console
user@machine:~$ mkdir new_code && cd new_code
user@machine:~$ julia --project=.
julia> ] add https://github.com/davidamaro/GroupFunctions.jl
julia> ] test
```

```julia
using Pkg
Pkg.add("GroupFunctions")
# Add other required packages similarly
```

### If Julia is not installed

- **Mac**: Use `juliaup`. Installing Julia via `brew` is not recommended.
- **Linux**: Use the appropriate package manager (e.g., `sudo pacman -S julia`).
- **Windows**: Use the recommended installer (specific details not provided).

## Usage

```julia
julia> using RandomMatrices
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

For more examples, see the file `test/runtests.jl`.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. A to-do list is included in the `todo.txt` file.

## License

Until the code from AbstractAlgebra.jl is removed from this package, the license will align with AbstractAlgebra.jl's.

## Citation

Pending.

