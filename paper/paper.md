---
title: "GroupFunctions.jl: computing individual entries of the irreducible representations of the unitary group U(d)"
tags:
  - Julia
  - group theory
  - representation theory
  - unitary group
  - quantum optics
authors:
  - name: "David Amaro-Alcalá"
    orcid: "0000-0001-8137-2161"
    affiliation: 1
  - name: "Konrad Szymański"
    orcid: "0000-0001-7676-1605"         
    affiliation: 1   
affiliations:
  - name: "Research Centre for Quantum Information, Institute of Physics, Slovak Academy of Sciences, Dúbravská cesta 9, Bratislava, Slovakia"
    index: 1
date: "2026-07-07"
bibliography: references.bib
---


# Summary
`GroupFunctions.jl`[^author-contributions] is a Julia library for computing individual matrix elements of irreducible representations of $\mathrm{U}(d)$. These matrix elements, called group functions, can be evaluated symbolically or numerically. For $\mathrm{SU}(2)$, they reduce to the Wigner $D$-functions. The library computes these matrix elements in a carrier-space basis enumerated by Gelfand-Tsetlin patterns [@GelfandTsetlin1950]. It can also compute entire representation operators, construct input unitaries from parameterisations common in quantum optics, translate Gelfand–Tsetlin patterns into occupation-number kets, and compute the associated Schur functions. Results can be exported in a form compatible with Mathematica.

[^author-contributions]: David Amaro-Alcalá wrote the package code, implemented the algorithms, developed the tests, and prepared the first version of the documentation. Konrad Szymański substantially revised and expanded the documentation, added examples, requested additions to the code API, and translated several Spanish-language function names and related documentation into English.
     
# Statement of need

Representations of the unitary group $\mathrm{U}(d)$ arise in many subfields of physics and mathematics, and computations often reduce to evaluating their matrix elements, called group functions. For an irrep labelled by $\lambda$, the corresponding object is

$$ D^{(\lambda)}_{\mathrm{out},\mathrm{init}}(U) = \langle \mathrm{out} \mid D^{(\lambda)}(U) \mid \mathrm{init} \rangle ,$$

here, $\mathrm{init}$ and $\mathrm{out}$ denote basis states in the representation carrier space.

In mathematics, summing the diagonal group functions gives the trace of the representation matrix of $U$. This trace is the character of the representation and the Schur polynomial of the eigenvalues of $U$, an object central to algebraic combinatorics and symmetric function theory.

In quantum physics, group functions appear in several settings. In quantum optics, a group function gives the transition amplitude of photons through a linear optical network. The same object helps characterise quantum devices [@amaroalcala2025] and describe the symmetry properties of states [@OttoSzymanski2024]. One important subproblem is boson sampling [@Aaronson2011], where the transition amplitude reduces to a permanent whose evaluation is classically computationally hard. Whether noisy real-world quantum devices can perform this computation remains an active question.

Some of these tasks require only a numerical estimate, whereas others require an exact symbolic group function. `GroupFunctions.jl` addresses the latter need. Although related packages exist (see below), to our knowledge, none is designed to compute individual representation-matrix entries symbolically. When applicable, computing an individual entry is more efficient than assembling the entire matrix.

## Similar software

Several existing packages relate to `GroupFunctions.jl` but address different problems. `SUNRepresentations.jl` [@SUNRepresentations] and the algorithm of @Alex2011 (with appendix code) compute $\mathrm{SU}(d)$ Clebsch-Gordan coefficients. `RepLAB` [@RepLAB] supports manipulating irreducible representations of various groups, including $\mathrm{U}(d)$, but provides only indirect numerical access to group functions. `IntegrateUnitary.jl` [@IntegrateUnitary] performs symbolic integration over compact groups rather than evaluating representation matrices. `haarpy` <https://github.com/polyquantique/haarpy> implements Weingarten-calculus methods. Other libraries focus on quantum optics. `BosonSampling.jl` [@Seron2024], `Perceval` [@Heurtel2023], and `QOptCraft` [@QOptCraft] numerically model linear optical devices, while `The Walrus` [@Gupt2019] helps compute amplitudes for Gaussian boson sampling. These packages do not target the symbolic computation of individual representation-matrix entries, which is the primary purpose of `GroupFunctions.jl`.

# Example use and documentation
The library primarily computes matrix elements of a representation between basis states, with auxiliary functions that support calculations in quantum optics. These matrix elements can be computed from a symbolic matrix. The following example evaluates a matrix element between states of the $\mathrm{U}(10)$ symmetric irrep corresponding to $7$ bosons.

```julia
λ = [7,0,0,0,0,0,0,0,0,0]; basis=basis_states(λ); # integer partition and basis  
init = findfirst(gt -> occupation_number(gt)==[0,0,0,1,1,1,1,1,1,1],basis); 
out  = findfirst(gt -> occupation_number(gt)==[1,1,1,1,1,1,1,0,0,0],basis);
group_function(λ, basis[init],basis[out]) # symbolic matrix element
```
`GroupFunctions.jl` is currently single-threaded. On an AMD Ryzen 7 PRO 4750U laptop CPU, the symbolic computation in the example took approximately 4 seconds after compilation warmup, without reusing cached intermediate or final results.

This functionality supports more complex calculations. For example, the library has been used to numerically evaluate $\mathrm{SU}(d)$ group characters in randomised benchmarking [@amaroalcala2025]. Further examples, applications, and mathematical background are available in the documentation: <https://davidamaro.github.io/GroupFunctions.jl/dev/>.


# Design choices
This library provides a unified method for computing representation matrix elements of $\mathrm{U}(d)$ irreps specified by integer partitions of length at most $d$, including computations with symbolic input matrices. Several computational routes are possible in principle. For symmetric irreps, which model fully indistinguishable bosons, one can manipulate states as polynomials of creation operators applied to the vacuum. Another approach constructs and exponentiates the Lie algebra generators in the chosen representation. Both approaches are computationally expensive. More restricted methods based on generating functions also exist [@prakash1996wigner].

The author of the original package chose the Grabmeier-Kerber formula [@Grabmeier1987] as the most general solution. It expresses the matrix element as a sum of monomials in the entries of the input matrix, weighted by the irreducible representation and the states in question. Our implementation optimises the enumeration of the double cosets that index this sum by grouping permutations that contribute the same monomial. This implementation provides the library's main function, `group_function`.

Internally, the algorithm represents basis states as semistandard Young tableaux. The user-facing functions expose the equivalent Gelfand-Tsetlin patterns through the `GTPattern` data structure and provide utility functions for common quantum optics scenarios, such as `occupation_number`.

# Availability
The library is available under the MIT licence and can be installed through the Julia package manager:

```julia
 ] add GroupFunctions
 ```

# AI usage disclosure
OpenAI Codex (GPT-5.3) assisted with optimising the performance of the double-coset enumeration at a late stage. The authors developed the mathematical design and proof of correctness of the optimised algorithm. Anthropic Claude (Opus 4.8) assisted with code review and language review of the documentation and manuscript. The authors reviewed and edited all AI-assisted changes.

# Acknowledgements
We thank Dr. Hubert de Guise for helpful discussions and suggestions on the bibliography, and Dr. Alonso Botero for suggestions that improved the presentation. Mitacs CALAREO, DeQHOST APVV-22-0570, QUAS VEGA 2/0164/25, Postdokgrant APD0161, and the Stefan Schwarz programme supported this work. David Amaro-Alcalá acknowledges the indirect support of the Government of Alberta and NSERC during his PhD studies at the University of Calgary.
  
# References
