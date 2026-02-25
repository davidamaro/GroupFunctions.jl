# Immanants

Immanants generalize both permanents and determinants, and they appear naturally in multiphoton interference when exchange symmetry is not purely bosonic or purely fermionic. In `GroupFunctions.jl`, these quantities are recovered from sums of matrix elements computed by `group_function`.

## Setup and notation

For a photon in mode $k$ centered at time $\tau_j$, use

```math
\hat{A}^\dagger_k(\tau_j)=
\int d\omega\,e^{i\omega\tau_j}\varphi(\omega)\hat{a}^\dagger_k(\omega),
```

with bosonic commutator

```math
[\hat{a}_i(\omega),\hat{a}^\dagger_j(\omega')] = \delta_{ij}\delta(\omega-\omega').
```

With a Gaussian spectrum,

```math
|\varphi(\omega)|^2=\frac{e^{-(\omega-\omega_0)^2/(2\sigma^2)}}{\sqrt{2\pi}\sigma},
```

the temporal overlap factor is

```math
\zeta_{ij}=e^{-\sigma^2\tau_{ij}^2}, \qquad \tau_{ij}=\tau_i-\tau_j.
```

For a matrix $M$, $M_{v \to u}$ denotes the submatrix with rows $u$ and columns $v$; repeated row or column labels are allowed (for example, an output pattern like `pp4`).

## Definition and special cases

For a partition $\lambda \vdash n$, the immanant of an $n \times n$ matrix $M$ is

```math
\mathrm{Imm}^{\lambda}(M)=\sum_{\pi \in S_n}\chi^\lambda(\pi)\prod_{i=1}^n M_{i,\pi(i)},
```

where $\chi^\lambda$ is the character of the irrep $\lambda$ of `S_n`.

Two important limits are

```math
\mathrm{Per}(M)=\mathrm{Imm}^{(n)}(M), \qquad
\mathrm{Det}(M)=\mathrm{Imm}^{(1,\ldots,1)}(M).
```

For $3 \times 3$ matrices and mixed symmetry $(2,1)$:

```math
\mathrm{Imm}^{(2,1)}(M)=2M_{11}M_{22}M_{33}-M_{12}M_{23}M_{31}-M_{13}M_{21}M_{32}.
```

## Example: permanent from a symmetric group function

This reproduces the relation tested in the package:

```julia
using GroupFunctions

permanent3(A) = A[1,1]*A[2,2]*A[3,3] +
                A[1,1]*A[2,3]*A[3,2] +
                A[1,2]*A[2,1]*A[3,3] +
                A[1,2]*A[2,3]*A[3,1] +
                A[1,3]*A[2,1]*A[3,2] +
                A[1,3]*A[2,2]*A[3,1]

alpha1, beta1, gamma1 = rand(Float64, 3)
block12_a = su2_block(4, 1, (alpha1, beta1, gamma1))
alpha2, beta2 = rand(Float64, 2)
block23_a = su2_block(4, 2, (alpha2, beta2, alpha2))
alpha3, beta3, gamma3 = rand(Float64, 3)
block12_b = su2_block(4, 1, (alpha3, beta3, gamma3))
alpha4, beta4 = rand(Float64, 2)
block34_a = su2_block(4, 3, (alpha4, beta4, alpha4))
alpha5, beta5 = rand(Float64, 2)
block23_b = su2_block(4, 2, (alpha5, beta5, alpha5))
alpha6, beta6, gamma6 = rand(Float64, 3)
block12_c = su2_block(4, 1, (alpha6, beta6, gamma6))

U = block12_a * block23_a * block12_b * block34_a * block23_b * block12_c

basis = basis_states([3,0,0,0])
state_x = filter(s -> pweight(s) == [0,1,1,1], basis)[1]
state_y = filter(s -> pweight(s) == [1,0,2,0], basis)[1]

M = U[[1,2,3], [2,2,4]]
group_function([3,0,0,0], state_x, state_y, U) ≈ permanent3(M) / sqrt(2)
```

## Two-photon coincidence rate with delay

For two photons entering input modes $k,l$ and detected at outputs $m,n$, with relative delay $\tau_{ab}$:

```math
R(kl \to mn;\tau_{ab})
=\frac{1}{2}(1+\zeta_{ab})\left|\mathrm{Per}(U_{kl\to mn})\right|^2
+\frac{1}{2}(1-\zeta_{ab})\left|\mathrm{Det}(U_{kl\to mn})\right|^2.
```

Limits:

```math
\tau_{ab}=0 \Rightarrow R=\left|\mathrm{Per}(U_{kl\to mn})\right|^2,
```

```math
\tau_{ab}\to\infty \Rightarrow
R=\frac{1}{2}\left|\mathrm{Per}(U_{kl\to mn})\right|^2
+\frac{1}{2}\left|\mathrm{Det}(U_{kl\to mn})\right|^2.
```

So the delay continuously interpolates between fully indistinguishable interference (permanent) and the incoherent permanent/determinant mixture.

## Reading the rate notation

The shorthand

```math
R(\text{inputs} \to \text{outputs}; \text{delays})
```

means: coincidence probability for a fixed input occupation pattern, output occupation pattern, and set of delays.

Examples:

- $R(23 \to 13; \tau_{ab})$: one photon enters input 2 and one enters input 3; detect one at output 1 and one at output 3.
- $R(234 \to \alpha\alpha4; \tau_{13})$: one photon enters each of inputs 2,3,4; detect two photons at output $\alpha$ and one at output 4.

Repeated output labels encode multiplicity.

## Three photons with one delayed input

For inputs $2,3,4$ in a four-mode interferometer, with $\tau_1=\tau_2 \neq \tau_3$, and output pattern `pp4`:

```math
R(234\to pp4;\tau_{13})
=\frac{1}{3}(1+2\zeta_{13})\left|\mathrm{Per}(U_{234\to pp4})\right|^2
+\frac{2}{3}(1-\zeta_{13})\left|\mathrm{Imm}^{(2,1)}(U_{234\to pp4})\right|^2.
```

For distinct outputs $\alpha \neq \beta$:

```math
\begin{aligned}
R(234\to \alpha\beta4;\tau_{13})
&=|A|^2+|B|^2+|C|^2 \\
&\quad+\zeta_{13}\left[(A+B)^*C+(B+C)^*A+(C+A)^*B\right],
\end{aligned}
```

with

```math
U_{xyz}^{(\alpha\beta4)} \equiv U_{xyz\to\alpha\beta4},\quad
P=\mathrm{Per}(U_{234}^{(\alpha\beta4)}),\quad
I_{xyz}=\mathrm{Imm}^{(2,1)}(U_{xyz}^{(\alpha\beta4)}),
```

and

```math
\begin{aligned}
A&=\frac{1}{3}\left(P-I_{243}-I_{324}+I_{342}\right),\\
B&=\frac{1}{3}\left(P-I_{234}+I_{243}-I_{324}-I_{342}\right),\\
C&=\frac{1}{3}\left(P+I_{234}+I_{324}\right).
\end{aligned}
```

Here $\tau_{13}$ is the only independent delay because $\tau_1=\tau_2$. At $\tau_{13}=0$, photons are fully indistinguishable and only the permanent term survives in the `pp4` rate.

## References

- [B. Kostant, "On Macdonald's $\eta$-function formula, the Laplacian and generalized exponents"](http://www.jstor.org/stable/2152885)
- [A. de Guise et al., "Coincidence landscapes for three-channel linear optical networks"](https://arxiv.org/pdf/1511.01851)
- [Immanant (Wikipedia)](https://en.wikipedia.org/wiki/Immanant)
