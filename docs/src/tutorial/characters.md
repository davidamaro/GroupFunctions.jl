```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Characters
SU(d) characters can be evaluated directly with `character(irrep, U)`.
They are the traces of the corresponding irreducible representation matrices, so
`character` is the direct way to ask for that quantity without building the
whole matrix yourself.

## Example: character as the trace of a representation

```@repl chars_trace
using GroupFunctions

λ = [2, 0];
U = su2_block(2, 1, (0.0, pi/3, 0.0));

χ = character(λ, U)

values, basis = group_function(λ, U);
χ_from_matrix = sum(values[i, i] for i in axes(values, 1))

χ ≈ χ_from_matrix
```

## Example: symbolic character

```@repl chars_symbolic
using GroupFunctions

character([2, 0])
```
