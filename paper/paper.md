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
  - name: "Research Center for Quantum Information, Institute of Physics, Slovak Academy of Sciences, Dúbravská cesta 9, Bratislava, Slovakia"
    index: 1
date: 29 July 2026
bibliography: paper.bib
---


# Summary
[`GroupFunctions.jl`](https://github.com/davidamaro/GroupFunctions.jl)[^author-contributions] is a Julia library for computing individual matrix entries of irreducible representations of $\mathrm{U}(d)$. These entries, called group functions, can be evaluated symbolically or numerically. For $\mathrm{SU}(2)$, they reduce to the Wigner $D$-functions. The library computes these matrix entries in a carrier-space basis enumerated by Gelfand-Tsetlin patterns [@GelfandTsetlin1950]. It can also compute entire representation operators, construct input unitaries from parameterisations common in quantum optics, translate Gelfand-Tsetlin patterns into occupation-number kets, and compute the associated Schur functions.

[^author-contributions]: David Amaro-Alcalá developed the package, its algorithms, and its test suite, and prepared the initial documentation. Konrad Szymański substantially revised and expanded the documentation and developed additional examples, and contributed to the standardisation of the API.
     
# Statement of need

Representations of the unitary group $\mathrm{U}(d)$ arise in many subfields of physics and mathematics, and computations often reduce to evaluating their matrix entries, called group functions. For an irreducible representation labelled by $\lambda$, the corresponding object is -- in the mathematicians' and physicists' inner product notation -- the following:

$$ D^{(\lambda)}_{\mathrm{out},\mathrm{init}}(U) = (D^{(\lambda)}(U) \mathrm{init},\mathrm{out}) =\langle \mathrm{out} \mid D^{(\lambda)}(U) \mid \mathrm{init} \rangle ,$$

where $\mathrm{init}$ and $\mathrm{out}$ denote basis vectors in the representation carrier space.

In mathematics, summing the diagonal group functions gives the trace of the representation matrix of $U$. This trace is the character of the representation and the Schur polynomial of the eigenvalues of $U$, an object central to algebraic combinatorics and symmetric function theory.

In quantum physics, group functions appear in several settings. In quantum optics, a group function gives the transition amplitude of photons through a linear optical network [@deGuise2014]. The same object can be used to characterise the performance of quantum devices [@amaroalcala2025] and to describe the symmetry properties of states [@OttoSzymanski2024]. One important subproblem is boson sampling [@Aaronson2011], where the transition amplitude reduces to a permanent whose evaluation is classically computationally hard. Whether noisy real-world quantum devices can perform this computation remains an open question.

Some of these tasks require only a numerical estimate, whereas others require an exact symbolic group function. `GroupFunctions.jl` addresses the latter. Although related software exists (see below), to our knowledge none is designed to compute a single representation-matrix entry as an exact symbolic $D$-function. We therefore compare the cost of computing one such entry with the conventional route that constructs the Lie-algebra representation matrices and exponentiates the resulting full matrix.

The cost of forming the full matrix depends on the dimension $D_\lambda(d)$ of the representation, the number of rows and columns in $D^{(\lambda)}(U)$. Fix a Young diagram $\lambda$ with $N>1$ boxes and at most $d$ rows. The basis vectors are indexed by semistandard Young tableaux of this shape with entries in $\{1,\ldots,d\}$. The hook-content formula counts these tableaux [@Fulton1997, §4.3, p. 55]. Equivalently, the dimension is the Schur polynomial at $d$ arguments equal to 1, because the character at the identity equals the dimension:

$$D_\lambda(d)=s_\lambda(\underbrace{1,\ldots,1}_{d})
=\prod_{(i,j)\in\lambda}\frac{d+j-i}{h_{ij}}.$$

Here $(i,j)$ denotes a box in row $i$ and column $j$, its content is $j-i$, and its hook length $h_{ij}$ counts the box itself, those to its right in the same row, and those below it in the same column. For fixed $\lambda$, the product has $N$ factors linear in $d$ and a constant denominator, so $D_\lambda(d)=\Theta(d^N)$. For example, $\lambda=(2)$ gives $D_{(2)}(d)=d(d+1)/2$.

To compute one entry directly, we keep the two tableau labels fixed as $d$ grows and supply the corresponding basis vectors, avoiding enumeration of the full basis. In this setting, the current symbolic implementation has time complexity $\Theta(d^3)$ and space complexity $\Theta(d^2)$. The cubic term comes from constructing the tables that index the monomial sum: assembling each length-$d$ row recursively copies successively longer arrays, costing $\Theta(d^2)$ operations per row across $d$ rows. With $N$ fixed, the number of tables and permutations remains bounded as $d$ grows. These estimates use a unit-cost model for array access, copying individual entries, and scalar arithmetic, with space measured by the maximum number of scalar entries stored simultaneously. The dependence on $N$ is absorbed into the constants; bit complexity is excluded.

For comparison, suppose the $D_\lambda(d)\times D_\lambda(d)$ Lie-algebra matrix has already been constructed. Numerical scaling and squaring with a Padé approximant uses dense matrix products and a linear solve [@Higham2005]. With the Padé degree and number of squarings fixed, standard dense multiplication and elimination require $\Theta(D_\lambda(d)^3)=\Theta(d^{3N})$ scalar arithmetic operations and have space complexity $\Theta(D_\lambda(d)^2)=\Theta(d^{2N})$. Here each scalar addition, subtraction, multiplication, or division has unit cost; in a numerical implementation, these operations use fixed-precision arithmetic. For $\lambda=(2)$, these costs grow as $d^6$ and $d^4$, respectively, compared with $d^3$ and $d^2$ for the direct calculation. The numerical route produces the full matrix at a chosen $U$; the direct route produces one exact polynomial in its entries. This comparison explains the cost of forming the full matrix when only one entry is needed. It applies with $N$ fixed and does not compare against numerical methods designed to compute individual entries.

# State of the field

Several existing packages relate to `GroupFunctions.jl` but address different problems. `SUNRepresentations.jl` [@SUNRepresentations] and the algorithm presented in the appendix of [@Alex2011] compute $\mathrm{SU}(d)$ Clebsch-Gordan coefficients. `RepLAB` [@RepLAB] supports manipulating irreducible representations of various groups, including $\mathrm{U}(d)$, but provides only indirect numerical access to group functions. `IntegrateUnitary.jl` [@IntegrateUnitary] performs symbolic integration over compact groups rather than evaluating representation matrices. `haarpy` [@cardin2024haarpy] implements Weingarten-calculus methods. Other libraries focus on quantum optics. `BosonSampling.jl` [@Seron2024], `Perceval` [@Heurtel2023], and `QOptCraft` [@QOptCraft] numerically model linear optical devices, while `The Walrus` [@Gupt2019] helps compute amplitudes for Gaussian boson sampling. These packages do not target the symbolic computation of individual representation-matrix entries, which is the primary purpose of `GroupFunctions.jl`.



# Software design
This library provides a unified method for computing representation matrix entries of $\mathrm{U}(d)$ irreducible representations specified by integer partitions of length at most $d$, including computations with symbolic input matrices. Several computational routes are possible in principle. For symmetric irreducible representations, which model fully indistinguishable bosons, one can manipulate states as polynomials of creation operators applied to the vacuum. Another approach constructs and exponentiates the Lie algebra generators in the chosen representation. Both approaches are computationally expensive. More restricted methods based on generating functions also exist [@prakash1996wigner].

The authors chose the Grabmeier-Kerber formula [@Grabmeier1987] as the most general solution. It expresses the entry of a matrix as a sum of monomials in the entries of the input matrix, weighted by the irreducible representation and the basis vectors in question. Our implementation optimises the enumeration of the double cosets that index this sum by grouping permutations that contribute the same monomial. This implementation provides the library's main function, `group_function`.

Internally, the algorithm represents basis vectors as semistandard Young tableaux. The user-facing functions expose the equivalent Gelfand-Tsetlin patterns through the `GTPattern` data structure and provide utility functions for common quantum optics scenarios, such as `occupation_number`.

# Research impact statement
`GroupFunctions.jl` has been used to evaluate matrix entries and group characters in published work on filtered randomized benchmarking [@amaroalcala2025]. The library has been under continuous development since 2020 and is tested in CI; it is registered in the Julia General registry under the MIT licence, installable with the following command:

```julia
 ] add GroupFunctions
 ```

 The following example evaluates, symbolically, an entry of a matrix between basis vectors of the $\mathrm{U}(4)$ mixed-symmetry irreducible representation.
```julia
λ = [2,2,1,0]; basis=basis_states(λ); # integer partition and basis
group_function(λ, basis[1], basis[end]) # symbolic matrix entry
```

The result is a polynomial in the entries of the $\mathrm{U}(4)$ matrix. For symmetric irreducible representations the entry of a matrix reduces to a permanent, which symbolic algebra packages also compute; for mixed-symmetry ones no existing software computes these entries symbolically. Further examples, applications, and mathematical background are available in [the documentation](https://davidamaro.github.io/GroupFunctions.jl/dev/).


# AI usage disclosure
OpenAI Codex (GPT-5.3) assisted with optimising the performance of the double-coset enumeration at a late stage. The authors developed the mathematical design and proof of correctness of the optimised algorithm. Anthropic Claude (Opus 4.8) assisted with code review and language review of the documentation and manuscript. The authors reviewed and edited all AI-assisted changes.

# Acknowledgements
We thank Dr. Hubert de Guise for helpful discussions and suggestions on the bibliography, and Dr. Alonso Botero for suggestions that improved the presentation. Mitacs CALAREO, DeQHOST APVV-22-0570, QUAS VEGA 2/0164/25, Postdokgrant APD0161, and the Štefan Schwarz programme supported this work. David Amaro-Alcalá acknowledges the indirect support of the Government of Alberta and NSERC during his PhD studies at the University of Calgary.
  
# References
