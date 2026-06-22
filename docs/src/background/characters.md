```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Characters

The character of an irrep is the trace of its representation matrix,
$\chi_\lambda(U) = \operatorname{tr}\rho_\lambda(U)$. Because the trace is
unchanged by conjugation ($\chi(U)=\chi(V^\dagger U V)$ for unitary $V$), $\chi_\lambda$ depends on
`U` only through its eigenvalues (it is a *class function*). `character(λ, U)` computes it by summing the
diagonal group functions over the GT basis (see [group
functions](../tutorial/group_functions.md)).

```@repl chars
using GroupFunctions

λ = [2, 1, 0];
U = su2_block(3, 1, (0.0, pi/3, 0.0));

character(λ, U)
```

## Orthogonality

Characters of distinct irreps are orthonormal under the Haar measure on U(d):

```math
\langle \chi_\lambda, \chi_\mu \rangle
= \int_{U(d)} \overline{\chi_\lambda(U)}\,\chi_\mu(U)\,dU
= \delta_{\lambda\mu}.
```

In particular, for a single irrep $\langle |\chi_\lambda|^2 \rangle = 1$. This
is a statement about averages, and we can exemplified it by sampling.

Here we draw Haar-random U(3) matrices (`Haar(2)` is the unitary ensemble, `2` does *not* denote the dimension) and
average $|\chi_\lambda|^2$. The average approaches `1`; the estimator is noisy,
so a small sample scatters by several percent, use $10^4$ samples or more:

```@repl chars_orth
using GroupFunctions
using RandomMatrices: Haar
using Statistics: mean
using Random: seed!

λ = [2, 1, 0];
nsamples = 10000;
seed!(420);

values = [character(λ, rand(Haar(2), 3)) for _ in 1:nsamples];
mean(abs2, values)
```

## The Weyl character formula

Computing $\chi_\lambda(U)$ does not require building $\rho_\lambda(U)$ (a matrix of size $\dim V_\lambda$) and tracing it. The Weyl character formula recovers the same trace from $\lambda$ and the eigenvalues $x_1,\dots,x_d$ of `U` alone, as a ratio of two $d\times d$ determinants. The denominator depends only on $d$:

```math
\det\begin{pmatrix}
x_1^{d-1} & x_1^{d-2} & \cdots & 1 \\
x_2^{d-1} & x_2^{d-2} & \cdots & 1 \\
\vdots & \vdots & & \vdots \\
x_d^{d-1} & x_d^{d-2} & \cdots & 1
\end{pmatrix}
= \prod_{i<j}(x_i - x_j),
```

the Vandermonde determinant. The numerator is the same matrix with the
partition added into the exponents, so the $(i,j)$-entry is
$x_i^{\lambda_j + d - j}$:

```math
\det\begin{pmatrix}
x_1^{\lambda_1+d-1} & x_1^{\lambda_2+d-2} & \cdots & x_1^{\lambda_d} \\
x_2^{\lambda_1+d-1} & x_2^{\lambda_2+d-2} & \cdots & x_2^{\lambda_d} \\
\vdots & \vdots & & \vdots \\
x_d^{\lambda_1+d-1} & x_d^{\lambda_2+d-2} & \cdots & x_d^{\lambda_d}
\end{pmatrix}.
```

Their ratio is the **Schur polynomial** $s_\lambda$:

```math
\chi_\lambda(U) = s_\lambda(x_1,\dots,x_d)
= \frac{\det\!\left(x_i^{\,\lambda_j+d-j}\right)}{\det\!\left(x_i^{\,d-j}\right)}.
```

For $d=2$, $\lambda=(1,0)$, the top exponents are $(2,0)$ and the ratio is
$\frac{x_1^2 - x_2^2}{x_1 - x_2} = x_1 + x_2$. The size of the irrep never
enters — two $d\times d$ determinants regardless of $\dim V_\lambda$.

This is what `schur_polynomial` computes; `character(λ, U)` instead builds and
traces the representation, so prefer `schur_polynomial` when only the character
is needed.

```@repl chars_schur
using GroupFunctions

schur_polynomial([2, 1, 0])
```
