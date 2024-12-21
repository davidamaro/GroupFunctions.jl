# GroupFunctions
[![Build Status](https://travis-ci.org/davidamaro/GroupFunctions.jl.svg?branch=master)](https://travis-ci.org/davidamaro/GroupFunctions.jl) 
[![Build status](https://ci.appveyor.com/api/projects/status/1l0wkv6isffameka?svg=true)](https://ci.appveyor.com/project/davidamaro/groupfunctions-jl)

# Group Functions and Matrix Operations Library

A Julia library for working with group functions, matrix operations, and special mathematical functions like immanants. This library provides functionality for working with SU(n) groups, Young tableaux, and related mathematical structures.

## Features

- Group function calculations for SU(n) groups
- Young tableau operations and manipulations
- Matrix operations with Haar random matrices
- Immanant calculations (using AbstractAlgebra.jl's character computation)
- Zero weight state detection
- Support for GT (Gelfand-Tsetlin) patterns
- Special factorizations and transformations

## Dependencies

The library requires the following Julia packages:

- GroupFunctions
- Test
- SymEngine
- AbstractAlgebra
- Immanants
- RandomMatrices (for Haar random matrices)
- LinearAlgebra

## Installation

To install the package, use Julia's package manager:

```julia
using Pkg
Pkg.add("GroupFunctions")
# Add other required packages similarly
```

## Usage

### Basic Operations

```julia
using GroupFunctions, Test, SymEngine, Immanants
import RandomMatrices: Haar
import LinearAlgebra: norm

# Find zero weight states for a partition
part = [2,1,1,0]
zeroweightstates = findzero(part)

# Work with Young Tableaux
t = YoungTableau([2,1])
fill!(t, [1,2,3])
```

### Group Functions

The library provides functionality to work with group functions and related operations:

```julia
# Create SU(n) matrices
α1, β1, γ1 = rand(Float64, 3)
xx = bloquesun(3,1,(α1,β1,γ1))

# Calculate group functions
result = group_function([2,1,0], state1, state2, matrix)
```

### Matrix Operations

```julia
# Generate Haar random matrices
mat = rand(Haar(2), 3)

# Work with special factorizations
matsimple = simple([α1,β1,γ1,α2,β2,α3,β3,γ3], 3)
```

## Testing

The library includes comprehensive tests covering various functionalities:

- SU(3) and SU(4) operations
- Immanant comparisons
- Sum rules verification
- Factorization tests
- Young tableau operations

Run tests using:

```julia
using Test
include("tests.jl")
```

## Features in Detail

### 1. Group Functions
- Support for various SU(n) group operations
- Implementation of specific irreducible representations
- Handling of complex matrix operations

### 2. Young Tableaux
- Creation and manipulation of Young diagrams
- Fill operations with various number sequences
- Pattern matching and validation

### 3. Matrix Operations
- Support for special unitary group matrices
- Implementation of Haar random matrices
- Custom factorization methods

### 4. Special Functions
- Immanant calculations leveraging AbstractAlgebra.jl's character computation capabilities
- Zero weight state detection
- GT pattern support
- Polynomial operations with SymEngine

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is open source and available under [LICENSE].

## Citation

If you use this software in your research, please cite:

[Citation information to be added]

## Contact

For questions and feedback, please [create an issue](https://github.com/yourusername/yourrepository/issues).
