```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Group Functions

`group_function` computes matrix elements between [basis states](states.md). Theory: [background page](../background/group_functions.md).

## Initialization

```@repl gf_pick
using GroupFunctions

λ = [2, 0];
basis = basis_states(λ);

foreach(gt -> println(pweight(gt)), basis)

initial = basis[findfirst(gt -> pweight(gt) == [1, 1], basis)];
final   = basis[findfirst(gt -> pweight(gt) == [2, 0], basis)];

initial
final
```

## Evaluate one numeric entry

```@repl gf_identity
using GroupFunctions
using LinearAlgebra: I

λ = [2, 0];
basis = basis_states(λ);
initial = basis[findfirst(gt -> pweight(gt) == [1, 1], basis)];
final   = basis[findfirst(gt -> pweight(gt) == [2, 0], basis)];

U = Matrix{ComplexF64}(I, 2, 2);
group_function(λ, final, initial, U)
group_function(λ, initial, initial, U)
```

## Beam splitter example

```@repl gf_bs
using GroupFunctions

λ = [2, 0];
basis = basis_states(λ);
initial = basis[findfirst(gt -> pweight(gt) == [1, 1], basis)];
final   = basis[findfirst(gt -> pweight(gt) == [2, 0], basis)];

BS = su2_block(2, 1, (0.0, pi/2, 0.0));
amp_11 = group_function(λ, initial, initial, BS)
amp_20 = group_function(λ, final, initial, BS)
round(abs2(amp_11), digits=6)
round(abs2(amp_20), digits=6)
```

## Compute the full (every entry) representation

```@repl gf_matrix
using GroupFunctions

λ = [2, 0];
BS = su2_block(2, 1, (0.0, pi/2, 0.0));
values, patterns = group_function(λ, BS);

size(values)
round.(values, digits=4)
```

## Compute a (single) symbolic entry

```@repl gf_symbolic
using GroupFunctions

λ = [2, 0];
basis = basis_states(λ);
initial = basis[findfirst(gt -> pweight(gt) == [1, 1], basis)];
final   = basis[findfirst(gt -> pweight(gt) == [2, 0], basis)];

group_function(λ, final, initial)
```

Next: [HOM effect](../applications/quantum_optics.md).
