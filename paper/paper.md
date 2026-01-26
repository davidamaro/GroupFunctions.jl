---
title: "GroupFunctions.jl: Short, Specific Subtitle Here"
tags:
  - Julia
  - group theory
  - computational algebra
  - passive interferometer
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

Provide a concise overview of the software, its purpose, and the main contributions.
Mention the domain context and what is novel or useful about this package.

The software compute both symbolical and numerical values of the irreducible
representations of the unitary group.
The software is based in well-established formulae.
Using group functions of the unitary group is widely useful in quantum
optics, quantum information, and quantum computing.
This is the first code to compute such group functions that does not need the
whole representation.

# Statement of need

Explain the problem this software addresses, why existing tools are insufficient,
and who benefits from this work.

Representations of the unitary group are needed in many branches in physics,
since the unitary group appears in many contexts. 
However, there is no implementation of the many ways in which it can be
computed.
Thus, this code addresses this problem by providing a function that can compute
such representations.
Current application of unitary irreps are  (benoit, mine, original,
bosonsampling).

# Software description

## Features

- Feature 1 and why it matters.
- Feature 2 with a short example use case.
- Feature 3 and its practical impact.


Group function for numerical and symbolical computations.
This is the main contribution, allows the computation of
individual irreducible representation entries for the unitary group.
These terms are common in applications such as the boson sampling problem
and passive interferometry.

Writing down states.
Allow to enumerate the basis.
Useful for both the computation of the matrix functions as well as the
notation of states in systems with unitary invariance.

Generation of unitary matrices.
This is also useful both in the content of bosonsampling and related tasks,
such as the characterisation of gates.



## Implementation details

Describe the core design choices, data structures, and algorithms.
If relevant, mention performance considerations and complexity.

Examples are provided as scripts (not notebooks), e.g., `examples/`.

core design choices
the code is an implementation of the formula in Grabmeier's paper.
And such, the code main functions involve:
1. Computing the double coset representatives.
2. Computing irreps of the symmetric group.
3. Working with semi standard young tableaux.

The first task requires obtaining every solution of a linear system of equation
with positive integer solutions.
This is done by using HighS.jl, because it is the robust open-source linear
solver. In previous versions of the code, we used Gurobi.
However, it was inconvenient for CI purposes to require a Gurobi license. 
So, I used Highs.jl as an alternative.

Data structures

The most important structure is `YoungTableau`,
which holds both standard and semi-standard gates.
This are used, in conjunction with `GTPattern`,
to label states.
However, for the software itself `YoungTableau`
this is more important as this is needed to compute the `Content`
and also the double coset representatives. These form the core of the
computation.


and algorithms.

There's no specific algorithm that we use.
In the future we expect to ditch Highs.jl to solve the numerical system.


# Quality control

Summarize tests, benchmarks, CI, and any validation against known results.
We can divide the tests by importance: code itself and the outcome.
For the code itself, the code to compute double coset representatives
and the order in the permutations is tested.
For the outcome, we use several test based on known results:
comparison with character formulae,
some polynomials, 
results for the Alex's paper.

# Availability

State the license, where to find the source, and how to install.
Include a short, reproducible install snippet if useful.

Same license as `AbstracAlgebra.jl`.
The source is hosted in Github: ``.
The instalation  can be done with `] add GroupFunctions`.

# Acknowledgements

<!--Credit collaborators, funding sources, or helpful discussions.-->
Helpful discussions with :
Dr Hubert de Guise (facilitating bibliography for
this and other versions),
Dr Alonso Botero (he suggested improvements in the presentation),
Dr Konrad Szymanski (improving the documentation).
Funding: Mitacs CALAREO, DeQHOST APVV-22-0570, and QUAS VEGA 2/0164/25.
Part of the work was done while I was a PhD student in the University of
Calgary,
thus, indirectly I also acknolwedge support from the Government of Albert
and NSERC.

# References
