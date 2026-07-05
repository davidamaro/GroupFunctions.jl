```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Basis states and GT patterns

GT patterns are the basis states used by `GroupFunctions.jl`. Each one labels a single
basis vector of an irrep; `basis_states` enumerates them, and the weight functions read
off what a pattern encodes. For *why* GT patterns are the right basis, see the
[background page](../background/states.md). Here we just build and use them.

We carry one example throughout: the irrep `λ = [2, 0, 0]` of U(3) — two bosons in three
modes.

## Build one directly

When you already know the pattern you want, construct it from its rows (top row first):

```@repl states
using GroupFunctions

gt = GTPattern([[2, 0, 0], [1, 0], [1]])
```

## Generate all patterns in an irrep

More often you want the whole basis. `basis_states` returns every pattern of the irrep,
in a fixed order:

```@repl states
λ = [2, 0, 0];
basis = basis_states(λ);

length(basis)
basis[1]
basis[2]
```

## Read the weights

A pattern's **p-weight** is its sequence of row-sum differences; `occupation_number` is
the same data in mode order (it is `reverse ∘ pweight`). For the symmetric irrep `[2,0,0]`
the occupation numbers are literally the Fock occupations:

```@repl states
gt = basis[1];

pweight(gt)
occupation_number(gt)
```

Listing the occupation numbers across the basis recovers every Fock state of two photons
in three modes:

```@repl states
occupation_number.(basis)
```

## Pick the states you need

To compute a transition amplitude you select an initial and a final pattern. `findfirst`
on the occupation numbers is the convenient way:


```@repl states
initial = basis[findfirst(gt -> occupation_number(gt) == [2, 0, 0], basis)];
final   = basis[findfirst(gt -> occupation_number(gt) == [1, 1, 0], basis)];
```

## Hand them to a group function

The states plug straight into `group_function`. Here `U` is a 50:50 beam splitter on
modes 1–2 (built from an SU(2) block); the result is the amplitude
$\langle 1,1,0 \mid U \mid 2,0,0\rangle$:

```@repl states
θ = float(π) / 2;
U = su2_block(3, 1, (0., θ, 0.));

group_function(λ, final, initial, U)
```

That is the whole pipeline: build patterns, pick two, evaluate. What the number *means*,
and the symbolic and multi-block cases, are covered next in
[group functions](group_functions.md).

## Mixed symmetry

Everything above used the symmetric irrep, where each pattern has a distinct p-weight. That
stops being true for mixed-symmetry irreps. Take `λ = [2, 1, 0]`:

```@repl states
basis21 = basis_states([2, 1, 0]);

occupation_number.(basis21)
```

Two different patterns share the occupation vector `[1, 1, 1]`; it is an *inner multiplicity*.
The weight alone no longer identifies the state, and you have to use the full pattern for unique determination. This is exactly
why the basis is GT patterns and not occupation numbers, and it is where permanents and
determinants give way to the general
[Grabmeier-Kerber formula](../background/group_functions.md#The-general-formula).

Next: [group functions](group_functions.md).
