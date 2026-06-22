```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Immanants

Immanants generalize both permanents and determinants, and they appear naturally in multiphoton interference when exchange symmetry is not purely bosonic or purely fermionic. In `GroupFunctions.jl`, these quantities are recovered from sums of matrix elements computed by `group_function`.

## Definition and special cases

+For an integer partition $\lambda=(\lambda_1,\ldots,\lambda_k)$ with $\sum_i \lambda_i = n$, the immanant of an $n \times n$ matrix $M$ is

```math
\mathrm{Imm}^{\lambda}(M)=\sum_{\pi \in S_n}\chi^\lambda(\pi)\prod_{i=1}^n M_{i,\pi(i)},
```

where $\chi^\lambda$ is the character of the irrep $\lambda$ of the symmetric group $S_n$ -- a function on permutations, distinct from the U(d) character $\chi^\lambda(U)$ of the [characters page](characters.md).

Two important limits are

```math
\mathrm{Per}(M)=\mathrm{Imm}^{(n)}(M), \qquad
\mathrm{Det}(M)=\mathrm{Imm}^{(1,\ldots,1)}(M).
```

For $3 \times 3$ matrices and mixed symmetry $(2,1)$:

```math
\mathrm{Imm}^{(2,1)}(M)=2M_{11}M_{22}M_{33}-M_{12}M_{23}M_{31}-M_{13}M_{21}M_{32}.
```


## Connection with group functions

A theorem of Kostant (see theorem 3 in [de Guise et al., D-functions and immanants of unitary matrices and
submatrices](https://arxiv.org/pdf/1511.01851)) gives the immanant directly as a sum of group functions.
With $U \in U(d)$ the unitary matrix and $\Gamma^{(\lambda)}(U)$ its matrix in the irrep $\lambda$,

```math
\mathrm{Imm}^{(\lambda)}(U) = \sum_t \langle t| \Gamma^{(\lambda)}(U) |t\rangle,
```
where the sum runs over the *zero-weight* basis states $|t\rangle$ (states with zero `zweight`, see [basis states](states.md)) of the irrep
$\lambda$, and $\langle t| \Gamma^{(\lambda)}(U) |t\rangle$ is a diagonal group function. The permanent
($\lambda$ a single row) and determinant (single column) are the extreme cases.
For the mixed $(2,1)$ immanant of a $3\times3$ matrix the $(11)$ irrep has two
zero-weight states, so the sum has two terms -- the relation verified on the
[tutorial page](../tutorial/immanants.md).
