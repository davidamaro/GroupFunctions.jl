```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# API reference

```@docs
GroupFunctions
```

```@index
```

## Group functions

```@docs
group_function
character
find_tableaux_fillings
find_double_coset_representative_matrices
```

## Gelfand-Tsetlin patterns

```@docs
GTPattern
basis_states
determine_next_pattern
determine_next_pattern!
zweight
pweight
occupation_number
```

## Young tableaux

```@docs
YoungTableau
axialdistance
determine_position
first_young_tableau_lexicographic
StandardYoungTableaux
generate_matrix
index_of_semistandard_tableau
content
Θ
```

## Matrix constructors and symbolic utilities

```@docs
su2_block
su2_block_symbolic
bs_block
bs_block_symbolic
swap_block
swap_block_symbolic
su2_factorization
sud_from_angles
julia_to_mma
mma_to_julia
expand
```

## Backwards compatibility

These names remain available for older callers, but the newer names above are
preferred in new code.

```@docs
bloquesun
simple
simplefactorization
indice_tablon_semistandard
genera_funcion
```

## Internal enumeration helpers

```@docs
GroupFunctions.AllSolutionsMatrix.enumerate_matrices
```
