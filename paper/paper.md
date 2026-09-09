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
[`GroupFunctions.jl`](https://github.com/davidamaro/GroupFunctions.jl)[^author-contributions] is a Julia library that computes individual matrix entries of irreducible representations of the unitary group $\mathrm{U}(d)$. A representation assigns a matrix to each group element and preserves multiplication. An irreducible representation has no nonzero proper subspace invariant under all its matrices. In a fixed basis, each matrix entry, viewed as a scalar function of the group element, is called a group function or $D$-function. For each requested group function, the user can choose either a symbolic expression or a numerical value at a specified input matrix. For $\mathrm{SU}(2)$, they reduce to the Wigner $D$-functions.

The basis vectors are labelled by Gelfand-Tsetlin patterns [@GelfandTsetlin1950]. The library also computes entire representation matrices, constructs input unitaries from parameterisations common in quantum optics, translates Gelfand-Tsetlin patterns into occupation-number kets, and computes the associated Schur functions.

[^author-contributions]: David Amaro-Alcalá developed the package, its algorithms, and its test suite, and prepared the initial documentation. Konrad Szymański substantially revised and expanded the documentation, developed additional examples, and contributed to the standardisation of the API.
     
# Statement of need

The library serves researchers in mathematics and quantum physics who need symbolic expressions or numerical values for representation matrix entries. For $U\in\mathrm{U}(d)$, let $D^{(\lambda)}(U)$ denote the matrix representing $U$ in an irreducible representation labelled by $\lambda$. The entry with row label $\mathrm{out}$ and column label $\mathrm{init}$ has equivalent expressions in inner-product and bra-ket notation:

$$ D^{(\lambda)}_{\mathrm{out},\mathrm{init}}(U) = (D^{(\lambda)}(U) \mathrm{init},\mathrm{out}) =\langle \mathrm{out} \mid D^{(\lambda)}(U) \mid \mathrm{init} \rangle ,$$

where $\mathrm{init}$ and $\mathrm{out}$ are vectors in the chosen orthonormal basis, and the inner product $(\cdot,\cdot)$ is linear in its first argument.

In mathematics, summing the diagonal group functions gives the character of the representation, namely the trace of the representation matrix of $U$. This character is the Schur polynomial evaluated at the eigenvalues of $U$, a central object in algebraic combinatorics and symmetric function theory.

In quantum optics, group functions give photon transition amplitudes through linear optical networks [@deGuise2014]. They also characterise the performance of quantum devices [@amaroalcala2025] and describe the symmetry properties of states [@OttoSzymanski2024]. In boson sampling [@Aaronson2011], the transition amplitude reduces to a permanent whose evaluation is classically computationally hard. Whether noisy real-world quantum devices can perform this computation remains an open question.

These tasks can require numerical estimates or exact symbolic group functions. Both are supported by `GroupFunctions.jl`; the following comparison concerns the symbolic case. To our knowledge, none of the related packages discussed below is designed to compute a single representation-matrix entry as an exact symbolic $D$-function. We therefore compare the cost of computing one entry directly with the conventional procedure of constructing the Lie-algebra representation matrices and exponentiating the resulting full matrix.

The cost of forming the full matrix depends on the representation dimension $r_\lambda(d)$, the number of rows and columns in $D^{(\lambda)}(U)$. Fix a Young diagram $\lambda$ with $N>1$ boxes and at most $d$ rows. The basis vectors are indexed by semistandard Young tableaux of this shape with entries in $\{1,\ldots,d\}$. The hook-content formula counts these tableaux [@Fulton1997, §4.3, p. 55]. The dimension $r_\lambda(d)$ of the irreducible representation labelled by $\lambda$ is also the Schur polynomial evaluated at $d$ arguments equal to 1, because a representation's character at the identity equals its dimension:

$$r_\lambda(d)=s_\lambda(\underbrace{1,\ldots,1}_{d})
=\prod_{(i,j)\in\lambda}\frac{d+j-i}{h_{ij}}.$$

The box $(i,j)$ in row $i$ and column $j$ has content $j-i$. Its hook length $h_{ij}$ counts the box itself, those to its right in the same row, and those below it in the same column. For fixed $\lambda$, the product has $N$ factors linear in $d$ and a constant denominator, giving $r_\lambda(d)=\Theta(d^N)$. Here $\Theta(d^k)$ denotes growth bounded above and below by positive constant multiples of $d^k$ as $d\to\infty$. For example, $\lambda=(2)$ gives $r_{(2)}(d)=d(d+1)/2$.

For the direct calculation, we keep the two tableau labels fixed as $d$ grows and supply the corresponding basis vectors, avoiding enumeration of the full basis. Under these assumptions, the current symbolic implementation has time complexity $\Theta(d^3)$ and space complexity $\Theta(d^2)$. The cubic cost arises from constructing the tables that index the monomial sum: recursively assembling each length-$d$ row copies successively longer arrays, requiring $\Theta(d^2)$ operations per row across $d$ rows. With $N$ fixed, the number of tables and permutations remains bounded as $d$ grows. The estimates assign unit cost to array access, copying individual entries, and scalar arithmetic, and measure space by the maximum number of scalar entries stored simultaneously. Constants may depend on $N$; bit complexity is excluded.

For comparison, consider constructing the representation matrix by exponentiating a Lie-algebra matrix [@BarutRaczka1986, p. 280]. Assume that this $r_\lambda(d)\times r_\lambda(d)$ Lie-algebra matrix has already been constructed. Numerical scaling and squaring with a Padé approximant uses dense matrix products and a linear solve [@Higham2005]. With the Padé degree and number of squarings fixed, standard dense multiplication and elimination require $\Theta(r_\lambda(d)^3)=\Theta(d^{3N})$ scalar arithmetic operations and have space complexity $\Theta(r_\lambda(d)^2)=\Theta(d^{2N})$. Each scalar addition, subtraction, multiplication, or division has unit cost; numerical implementations use fixed-precision arithmetic. For $\lambda=(2)$, these costs grow as $d^6$ and $d^4$, respectively, compared with $d^3$ and $d^2$ for the direct calculation. Exponentiation produces the full matrix at a chosen $U$; the direct calculation produces one exact polynomial in its entries. The comparison quantifies the cost of forming the full matrix when only one entry is needed, with $N$ fixed. Numerical methods designed to compute individual entries are outside this comparison.

# State of the field

Related packages support representation theory and quantum optics, with different computational goals. `SUNRepresentations.jl` [@SUNRepresentations] and the algorithm in the appendix of [@Alex2011] compute $\mathrm{SU}(d)$ Clebsch-Gordan coefficients. `RepLAB` [@RepLAB] manipulates irreducible representations of various groups, including $\mathrm{U}(d)$, but provides only indirect numerical access to group functions. `IntegrateUnitary.jl` [@IntegrateUnitary] performs symbolic integration over compact groups, and `haarpy` [@cardin2024haarpy] implements Weingarten-calculus methods. In quantum optics, `BosonSampling.jl` [@Seron2024], `Perceval` [@Heurtel2023], and `QOptCraft` [@QOptCraft] numerically model linear optical devices, while `The Walrus` [@Gupt2019] computes amplitudes for Gaussian boson sampling. None of these packages targets the symbolic computation of individual representation-matrix entries, the primary purpose of `GroupFunctions.jl`.



# Software design
The library computes matrix entries of $\mathrm{U}(d)$ irreducible representations specified by integer partitions of length at most $d$, using a unified method that also supports symbolic input matrices. Several computational routes are possible. For symmetric irreducible representations, which model fully indistinguishable bosons, states can be manipulated as polynomials of creation operators applied to the vacuum. Another approach constructs and exponentiates the Lie algebra generators in the chosen representation. Both approaches are computationally expensive. More restricted methods use generating functions [@prakash1996wigner].

The authors chose the Grabmeier-Kerber formula [@Grabmeier1987] for its generality. The formula expresses a requested matrix entry as a sum of monomials in the input matrix entries, with coefficients determined by the irreducible representation and the chosen basis vectors. Our implementation optimises the enumeration of the double cosets that index this sum by grouping permutations that contribute the same monomial. The library's main function, `group_function`, evaluates the resulting sum.

Internally, the algorithm represents basis vectors as semistandard Young tableaux. The public interface exposes the equivalent Gelfand-Tsetlin patterns through the `GTPattern` data structure and provides utilities such as `occupation_number` for common quantum optics calculations.

# Research impact statement
`GroupFunctions.jl` has been used to evaluate matrix entries and group characters in published work on filtered randomized benchmarking [@amaroalcala2025]. Development has continued since 2020, with tests run through continuous integration. The package is registered in the Julia General registry, distributed under the MIT licence, and installed with the following command:

```julia
 ] add GroupFunctions
 ```

The following example computes a symbolic matrix entry between two basis vectors of a mixed-symmetry irreducible representation of $\mathrm{U}(4)$.
```julia
λ = [2,2,1,0]; basis=basis_states(λ); # integer partition and basis
group_function(λ, basis[1], basis[end]) # symbolic matrix entry
```

The result is a polynomial in the entries of the $\mathrm{U}(4)$ matrix. For symmetric irreducible representations, matrix entries reduce to permanents, which symbolic algebra packages also compute. For mixed-symmetry representations, no other software computes these entries symbolically. Further examples, applications, and mathematical background are available in [the documentation](https://davidamaro.github.io/GroupFunctions.jl/dev/).


# AI usage disclosure
OpenAI Codex (GPT-5.3) assisted with optimising double-coset enumeration at a late stage of development. The authors developed the mathematical design and proof of correctness of the optimised algorithm. Anthropic Claude (Opus 4.8) assisted with code review and language editing of the documentation and manuscript. The authors reviewed and edited all AI-assisted changes.

# Acknowledgements
We thank Dr. Hubert de Guise for discussions and bibliographic suggestions, and Dr. Alonso Botero for suggestions on the presentation. Mitacs CALAREO, DeQHOST APVV-22-0570, QUAS VEGA 2/0164/25, Postdokgrant APD0161, and the Štefan Schwarz programme supported this work. David Amaro-Alcalá also acknowledges indirect support from the Government of Alberta and NSERC during his PhD studies at the University of Calgary.
  
# References
