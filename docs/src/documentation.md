
```@index
```


## Functions

```@docs
    group_function
```

## [Gelfand-Tsetlin patterns](tutorial/states.md)

[GT patterns](tutorial/states.md) are used to denote the basis states.

```@docs
    GTPattern
```

## Construction of matrices

SU(2) blocks are used to construct unitary matrices.

```@docs
    su2_block
```

## Young Tableaux Operations

Young tableaux are fundamental combinatorial objects used in representation theory. These functions provide tools for constructing and manipulating Young tableaux.

```@docs
    YoungTableau
    StandardYoungTableaux
    index_of_semistandard_tableau
```

## [GT Pattern Generation](tutorial/states.md)

Functions for generating and iterating through [Gelfand-Tsetlin patterns](tutorial/states.md), which form a basis for irreducible representations.

```@docs
    basis_states
```

## Weight Functions

Weight functions extract different types of weights from [Gelfand-Tsetlin patterns](tutorial/states.md), useful for quantum optics and physical applications.

```@docs
    zweight
    pweight
```

## SU(n) Matrix Construction

Advanced functions for constructing special unitary matrices using various factorization methods and Euler angle parametrizations.

```@docs
    su2_factorization
    sud_from_angles
    su2_block_symbolic
```


## Utility Functions

Conversion utilities for interfacing with other symbolic computation systems.

```@docs
    julia_to_mma
    mma_to_julia
```
