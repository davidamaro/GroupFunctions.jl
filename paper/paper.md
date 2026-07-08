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
`GroupFunctions.jl` is a Julia library for computations involving irreducible representations of $\mathrm{U}(d)$. Its core is the evaluation of individual group functions -- the matrix elements of an irrep between basis states, which for $\mathrm{SU}(2)$ reduce to the Wigner $D$-functions -- either symbolically or numerically. The matrix elements are computed between the vectors belonging to a carrier space basis enumerated via Gelfand-Tsetlin patterns [@GelfandTsetlin1950]. The library can also compute entire representation operators, construct input unitaries from parameterizations common in quantum optics, translate Gelfand–Tsetlin patterns into occupation-number kets, and compute the associated Schur functions. Results can be exported to a form compatible with Mathematica.
     
# Statement of need

Representations of the unitary group $\mathrm{U}(d)$ arise in many subfields of physics and mathematics, and computations often reduce to evaluating their matrix elements, called group functions. For an irrep labelled by $\lambda$, the corresponding object is

$$ D^{(\lambda)}_{\mathrm{out},\mathrm{init}}(U) = \langle \mathrm{out} \mid D^{(\lambda)}(U) \mid \mathrm{init} \rangle ,$$

where $\mathrm{init}$ and $\mathrm{out}$ denote basis states of the representation carrier space.

In mathematics, the trace of the representation matrix of $U$, obtained by summing the diagonal group functions, is the character of the representation; this character is the Schur polynomial of the eigenvalues of $U$. This object is central to algebraic combinatorics and symmetric function theory.

In quantum physics, group functions appear in several settings. In quantum optics, a group function is the transition amplitude of photons through a linear optical network; the same object helps the characterization of quantum devices [@amaroalcala2025] and the description of symmetry properties of states [@OttoSzymanski2024]. An important subproblem is boson sampling [@Aaronson2011], where the transition amplitude reduces to a permanent, whose evaluation is classically computationally hard; whether noisy real life quantum devices can perform this computation is an active question. 

For some of these tasks a numerical estimate suffices; for others an exact and symbolic group function is required. This need is answered by the `GroupFunctions.jl` library: related packages exist (see below), but to our knowledge, none of these tools is designed to compute individual group function entries symbolically, which is more efficient than assembling the entire matrix, when applicable.

## Similar software

Several  existing packages are related to `GroupFunctions.jl`, but address different problems. `SUNRepresentations.jl` [@SUNRepresentations] and the algorithm of @Alex2011 (with appendix code) compute $\mathrm{SU}(d)$ Clebsch-Gordan coefficients. `RepLAB` [@RepLAB] allows for manipulation of irreducible representations of various groups including $\mathrm{U}(d)$, but group functions are only accessible as indirect numerics. `IntegrateUnitary.jl` [@IntegrateUnitary] performs symbolic integration over compact groups rather than evaluating representation matrices. There also exist libraries oriented towards quantum optics: the aim of `BosonSampling.jl` [@Seron2024], `Perceval` [@Heurtel2023], and `QOptCraft` [@QOptCraft] is numerical modelling of linear optical devices rather than symbolic computation, and `The Walrus` [@Gupt2019] helps with computation of amplitudes for Gaussian boson sampling. These packages do not target the symbolic computation of individual group function entries, which is the purpose of `GroupFunctions.jl`.

# Example use and documentation
The core functionality of the library is to compute the matrix element of a representation between basis states, with auxiliary functions that ease computations in quantum optics. The calculation can be done with a symbolic matrix, as in the following example of a matrix element between states of the $\mathrm{U}(10)$ symmetric irrep corresponding to $7$ bosons.

```julia
λ = [7,0,0,0,0,0,0,0,0,0]; basis=basis_states(λ); # integer partition and basis  
init = findfirst(gt -> occupation_number(gt)==[0,0,0,1,1,1,1,1,1,1],basis); 
out  = findfirst(gt -> occupation_number(gt)==[1,1,1,1,1,1,1,0,0,0],basis);
group_function(λ, basis[init],basis[out]) # symbolic matrix element
```
`GroupFunctions.jl` is currently single-threaded. On an AMD Ryzen 7 PRO 4750U laptop CPU, the symbolic computation in the example took approximately 4 seconds after compilation warmup, without reusing cached intermediate or final results.

More complex calculations can be built on top of this, as in [@amaroalcala2025], where the library was used for numerical evaluation of $\mathrm{SU}(d)$ group characters in randomized benchmarking. Further examples, applications, and mathematical background  are available in the documentation, <https://davidamaro.github.io/GroupFunctions.jl/dev/>.


# Design choices
The goal of this library is a unified computation of representation matrix elements for irreps of $\mathrm{U}(d)$ specified by integer partitions of length at most $d$, including symbolic input matrices. Several computational routes are possible in principle: for symmetric irreps (modeling fully indistinguishable bosons), states can be manipulated as polynomials of creation operators applied to the vacuum. One could construct the Lie algebra generators in the chosen representation and exponentiate them. Both approaches mentioned are computationally expensive. Restricted approaches using generating functions also exist [@prakash1996wigner].

The Grabmeier-Kerber formula [@Grabmeier1987] was chosen as the most general solution: the matrix element is a sum of monomials in the entries of the input matrix, with weights stemming from the irreducible representation and the states in question. The authors implemented the formula with an optimized enumeration of the double cosets that index the sum, grouping permutations contributing the same monomial; the result is the main functionality of the library, `group_function`.

Internally, the algorithm represents basis states by semistandard Young tableaux. The user-facing functions expose the equivalent Gelfand-Tsetlin patterns as the `GTPattern` data structure, together with utility functions for common quantum optics scenarios, such as `occupation_number`.

# Availability
The library is available under the MIT licence and can be installed through the Julia package manager:

```julia
 ] add GroupFunctions
 ```

# AI usage disclosure
OpenAI Codex (GPT-5.3) was used to assist with the implementation of Gelfand–Tsetlin pattern pretty-printing and late-stage performance optimization of the double-coset enumeration. The mathematical design and proof of correctness of the optimized algorithm were developed by the authors. Anthropic Claude (Opus 4.8) was used for code review and for language review of the documentation and manuscript. All AI-assisted changes were reviewed and edited by the authors.

# Acknowledgements
We thank Dr. Hubert de Guise for helpful discussions and bibliography suggestions, and Dr. Alonso Botero for suggestions that improved the presentation. This work was supported by Mitacs CALAREO, DeQHOST APVV-22-0570, and QUAS VEGA 2/0164/25, Postdokgrant APD0161, and Stefan Schwarz programme.  David Amaro-Alcalá acknowledges his time as a PhD student at the University of Calgary, supported indirectly by the Government of Alberta and NSERC.
  
# References
