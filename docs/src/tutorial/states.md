# Basis states and GT patterns

GT patterns are the basis states used by `GroupFunctions.jl`.

## Build one directly

```@repl states_build
using GroupFunctions

gt = GTPattern([[2, 1, 0], [2, 1], [2]]);
gt
```

## Generate all patterns in an irrep

```@repl states_basis
using GroupFunctions

λ = [2, 1, 0];
basis = basis_states(λ);

length(basis)
basis[1]
basis[2]
basis[3]
```

## Read the weights

```@repl states_weights
using GroupFunctions

gt = basis_states([2, 1, 0])[1];

pweight(gt)
zweight(gt)
```

## Bosonic occupation numbers

For symmetric irreps, `pweight` is the occupation vector.

```@repl states_bosons
using GroupFunctions

λ = [2, 0, 0];
basis = basis_states(λ);

foreach(gt -> println(pweight(gt)), basis)
```

## Step to the next pattern

```@repl states_next
using GroupFunctions

gt = GTPattern([[2, 1, 0], [2, 1], [2]]);

determine_next_pattern(gt)
```

## Pick states for a group-function calculation

```@repl states_group_function
using GroupFunctions
using LinearAlgebra: I

λ = [2, 0, 0];
basis = basis_states(λ);

initial = basis[findfirst(gt -> pweight(gt) == [1, 1, 0], basis)];
final   = basis[findfirst(gt -> pweight(gt) == [2, 0, 0], basis)];

U = Matrix{ComplexF64}(I, 3, 3);
group_function(λ, final, initial, U)
```

Next: [group functions](group_functions.md).
