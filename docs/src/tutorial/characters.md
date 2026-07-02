```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Characters

`character(λ, U)` returns the trace of the irrep, see the
[background page](../background/characters.md) for the theory.

## Character as a trace

The character is the trace of the representation matrix, so it equals the sum
of the diagonal group functions (and is indeed defined in the code as such):

```@raw html
<p style="text-align: center;">
  <img style="width: 63%; height: auto;" src="../../assets/tutorial/tutorial_character_trace.svg" alt="The character as a sum over diagonal input-to-output transitions">
</p>
```

```@repl chars
using GroupFunctions

λ = [2, 0];
U = su2_block(2, 1, (0.0, pi/3, 0.0));

χ = character(λ, U)

values, basis = group_function(λ, U);
χ ≈ sum(values[i, i] for i in axes(values, 1))
```

## Symbolic character

Omit `U` for the character as a symbolic polynomial:

```@repl chars_sym
using GroupFunctions

character([2, 0])
```

## Schur polynomial on eigenvalues

A character depends only on the eigenvalues of `U`. `schur_polynomial`
evaluates it directly on them (without a costly `group_function` computation), matching `character`:

```@repl chars_schur
using GroupFunctions
using LinearAlgebra: eigvals

λ = [2, 1, 0];
U = su2_block(3, 1, (0.0, pi/3, 0.0));

character(λ, U)
schur_polynomial(λ, eigvals(U))
character(λ, U) ≈ schur_polynomial(λ, eigvals(U))
```
