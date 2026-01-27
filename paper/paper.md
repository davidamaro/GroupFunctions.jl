---
title: "GroupFunctions.jl: code base to compute entries of the irreducible
representations of the special unitary group SU(d)"
tags:
  - Julia
  - group theory
  - computational algebra
  - unitary group
authors:
  - name: "Amaro-Alcala, David"
    orcid: "0000-0001-8137-2161"
    affiliation: 1
affiliations:
  - name: "Research Centre for Quantum Information, Institute of Physics, Slovak Academy of Sciences, Dúbravská cesta 9, Bratislava 845 11, Slovakia"
    index: 1
date: "2026-01-26"
bibliography: references.bib
---

# Summary

<!--Provide a concise overview of the software, its purpose, and the main contributions.-->
<!--Mention the domain context and what is novel or useful about this package.-->

The software computes both symbolic and numerical matrix elements of
irreducible representations of the unitary group.
It is based on well-established formulae.
Group functions of the unitary group are widely useful in quantum optics,
quantum information, and quantum computing.
To our knowledge, this is the first public implementation that can compute
individual group-function entries without constructing the full
representation.

# Statement of need

<!--Explain the problem this software addresses, why existing tools are insufficient,-->
<!--and who benefits from this work.-->

Representations of the unitary group are needed in many branches of physics
because the unitary group appears in many contexts.
However, there is no widely available implementation covering the practical
computations used in applications.
This package addresses that gap by providing tools to compute these
representation-theoretic quantities directly.
Current applications of unitary irreps include boson sampling and related
photonic protocols [@Aaronson2011], as well as recent applications in our
work [@AmaroAlcal2020; @amaroalcala2025].

# Software description

## Features

<!--- Feature 1 and why it matters.-->
<!--- Feature 2 with a short example use case.-->
<!--- Feature 3 and its practical impact.-->


Group functions for numerical and symbolic computations.
This is the main contribution: it enables the computation of individual
irreducible-representation entries for the unitary group.
These matrix elements are common in applications such as boson sampling and
passive interferometry.
Here are two SU(2) examples, common in undergraduate quantum mechanics
classes.
```julia
julia> group_function([1,0])[1]
2×2 Matrix{SymEngine.Basic}:
 u_2_2  u_2_1
 u_1_2  u_1_1

julia> group_function([2,0])[1]
3×3 Matrix{SymEngine.Basic}:
             u_2_2^2        sqrt(2)*u_2_2*u_2_1              u_2_1^2
 sqrt(2)*u_2_2*u_1_2  u_2_1*u_1_2 + u_2_2*u_1_1  sqrt(2)*u_2_1*u_1_1
             u_1_2^2        sqrt(2)*u_1_1*u_1_2              u_1_1^2
```

Writing down states.
This allows enumeration of the basis.
It is useful both for computing matrix elements and for a consistent notation
of states in systems with unitary invariance.
Users of the package may also benefit from gaining familiarity with
Gelfand-Tsetlin patterns.
```julia
julia> states = basis_states([2,1,0])
8-element Vector{GTPattern}:
 GTPattern([[2, 1, 0], [1, 0], [0]], [0])
 GTPattern([[2, 1, 0], [1, 0], [1]], [1])
 GTPattern([[2, 1, 0], [1, 1], [1]], [1])
 GTPattern([[2, 1, 0], [2, 0], [0]], [0])
 GTPattern([[2, 1, 0], [2, 0], [1]], [1])
 GTPattern([[2, 1, 0], [2, 0], [2]], [2])
 GTPattern([[2, 1, 0], [2, 1], [1]], [1])
 GTPattern([[2, 1, 0], [2, 1], [2]], [2])

julia> states[3]
│ 2   1   0 ╲
│   1   1    〉
│     1     ╱
```

Generation of unitary matrices.
This is also useful in the context of boson sampling and related tasks, such
as the characterization of gates.

For example, we can construct an SU(2) matrix from Euler angles with
`sud_from_angles` and evaluate a transition amplitude in the symmetric irrep:

```julia
using GroupFunctions

lambda = [2, 0]
basis = basis_states(lambda)
angles = [0.0, pi/2, 0.0]           # length = 2^2 - 1
U = sud_from_angles(angles, 2)

initial = basis[findfirst(b -> pweight(b) == [1, 1], basis)]
final   = basis[findfirst(b -> pweight(b) == [2, 0], basis)]

amp = group_function(lambda, final, initial, U)
abs2(amp)
```

The same constructor works for higher dimensions. For SU(3), the angle vector
has length `3^2 - 1 = 8`:

```julia
using GroupFunctions

angles3 = [0.3, 1.1, -0.4, 0.7, 0.5, -0.8, 0.9, 0.2]
U3 = sud_from_angles(angles3, 3)

lambda3 = [2, 0, 0]
basis3 = basis_states(lambda3)
group_function(lambda3, basis3[1], basis3[2], U3)
```

## Examples

The following short scripts can be run directly as plain Julia scripts (no notebooks).

### Example: compute a single D-function entry and a character

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det

irrep = [2, 1, 0]
U = rand(Haar(2), 3)                # 3x3 Haar-random unitary
basis = basis_states(irrep)

group_function(irrep, basis[1], basis[3], U)
character_weyl(irrep, U / det(U)^(1/3))
```

### Example: Hong-Ou-Mandel interference (50:50 beamsplitter)

```julia
using GroupFunctions

lambda = [2, 0]
basis = basis_states(lambda)

initial = basis[findfirst(b -> pweight(b) == [1, 1], basis)]
final   = basis[findfirst(b -> pweight(b) == [2, 0], basis)]

BS = su2_block(2, 1, (0.0, pi/2, 0.0))
amp = group_function(lambda, final, initial, BS)
abs2(amp)  # ~0.5 for a 50:50 beamsplitter
```

### Example: SU(3) character check

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det, tr

irrep = [2, 1, 0]
U = rand(Haar(2), 3)
U_su = U / det(U)^(1/3)

rep, _ = group_function(irrep, U_su)
isapprox(tr(rep), character_weyl(irrep, U_su))
```

Additional examples in the documentation and tests include:
- Enumerating Gelfand-Tsetlin basis states and mapping them to p-weights.
- SU(3) character variance checks over Haar-random ensembles.
- SU(2)-invariant qubit transmission via symbolic blocks and `group_function_sym`.
- Immanant/permanent relations and sum-rule rate identities for coincidence counts.
- Construction of `GTPattern` objects and symbolic `group_function` entries.


## Implementation details

<!--Describe the core design choices, data structures, and algorithms.-->
<!--If relevant, mention performance considerations and complexity.-->

Examples are provided in `docs/` and `test/`, and are rendered in the
documentation site built by GitHub Actions.

<!--core design choices-->
The code implements the formula in Grabmeier and Kerber's paper
[@Grabmeier1987].
As such, the main computations involve:
1. Computing the double coset representatives.
2. Computing irreps of the symmetric group.
3. Working with semistandard Young tableaux.

The first task requires obtaining every solution of a linear system of
equations with positive integer solutions.
This is done using `HiGHS.jl`, a robust open-source linear solver.
In previous versions of the code, we used Gurobi.
However, it was inconvenient for CI purposes to require a Gurobi license,
so we adopted `HiGHS.jl` as an alternative.

<!--Data structures-->

The most important structure is `YoungTableau`,
which holds both standard and semistandard tableaux.
These are used, in conjunction with `GTPattern`,
to label states.
However, for the core computation, `YoungTableau` is especially important
because it is needed to compute `Content` values and the double coset
representatives. These form the core of the computation.


<!--and algorithms.-->

There is no single custom algorithm beyond the representation-theoretic
construction described above.
In the future, we may replace `HiGHS.jl` for the integer-system step if a
more specialized approach becomes available.


# Quality control

Tests cover both internal correctness and physics-facing outcomes.
Internally, we test the computation of double coset representatives and the
ordering conventions in permutations.
For outcomes, we use tests based on known results, including comparisons with
character formulae, polynomial identities, and results from the literature.

# Availability

<!--State the license, where to find the source, and how to install.-->
<!--Include a short, reproducible install snippet if useful.-->

The software uses the same license as `AbstractAlgebra.jl`.
The source is hosted on GitHub.
Installation can be done within the Julia package manager: `] add GroupFunctions`.

# Dependencies

- `AbstractAlgebra.jl`, for permutation and Young tableau manipulation.
- `JuMP.jl` and `HiGHS.jl`, for finding all solutions of a system of linear
  equations in $\mathbb{Z}^+$.
- `SymEngine.jl`, to handle symbolic output.

# Acknowledgements

<!--Credit collaborators, funding sources, or helpful discussions.-->
Helpful discussions with:
Dr. Hubert de Guise (facilitating bibliography for this and other versions),
Dr. Alonso Botero (suggesting improvements in the presentation), and
Dr. Konrad Szymanski (improving the documentation).
We also thank Diana Cencer Garafova for the logo design.
Funding: Mitacs CALAREO, DeQHOST APVV-22-0570, and QUAS VEGA 2/0164/25.
Part of the work was done while I was a PhD student in the University of
Calgary,
thus, indirectly, I also acknowledge support from the Government of Alberta
and NSERC.

# References
