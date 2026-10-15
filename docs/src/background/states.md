```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Basis vectors and Gelfand-Tsetlin patterns

This library computes matrix entries of U(d) irreducible representations. In quantum optics, for example, these entries describe how passive linear networks transform multiphoton states. Gelfand-Tsetlin (GT) patterns label the basis vectors between which these matrix entries are evaluated and sometimes give those vectors a physical interpretation.

Young diagrams commonly identify U(d) irreducible representations. A GT pattern records both the representation and a specific vector within it. More precisely, its rows record the integer partitions obtained by restricting the irreducible representations along the chain $U(d) \supset U(d-1) \supset \cdots \supset U(1)$. GT patterns thus provide canonical labels for basis vectors, making them useful in subsequent computations.

The following standard definitions establish the basic picture.

**Definitions**. Fix a positive integer $d$. The group $U(d)$ consists of the $d\times d$ unitary matrices. A *representation* is a homomorphism $\rho: U(d) \rightarrow GL(V)$ satisfying $\rho(u \cdot u')=\rho(u)\rho(u')$, where $V$ is a finite-dimensional complex vector space and $GL(V)$ is the group of invertible linear maps from $V$ to itself.

A representation is *irreducible* (an *irrep*) if the zero subspace and $V$ itself are its only subspaces invariant under every $\rho(u)$. A *partition* is a nonincreasing integer list $\lambda=[\lambda_1 \ge \cdots \ge \lambda_d]$, or equivalently, a Young diagram with $\lambda_j$ boxes in row $j$. Every finite-dimensional irreducible representation of $U(d)$ corresponds to exactly one partition $\lambda$. We denote the representation by $\rho_\lambda$ and its representation space by $V_\lambda$; the partition becomes the top row of a GT pattern.

**GT-pattern origin** 

To enumerate an orthonormal basis of $V_\lambda$ with particularly useful properties, we restrict $U(d)$ first to $U(d-1)$ and then recursively to $U(d-2)$, continuing until we reach $U(1)$. Because the irreducible representations of $U(1)$ are one-dimensional, their representation spaces define the desired basis vectors.

To restrict $\rho_\lambda$ to $U(d-1)$, evaluate it only on *block-diagonal matrices* of the form

```math
u = \begin{pmatrix} u' & 0 \\ 0& 1\end{pmatrix}=u'\oplus (1), 
```

where $u'$ is a unitary matrix of dimension $d-1$. The restricted representation sends $u'$ to $\rho(u' \oplus (1))$ and generally decomposes into irreducible representations of $U(d-1)$. This decomposition is multiplicity-free:

```math
V_\lambda\big|_{U(d-1)} \cong \bigoplus_{\mu\,\prec\,\lambda} V_\mu.
```

The standard Gelfand--Tsetlin branching rule identifies the labels in this decomposition as exactly those $\mu=(\mu_1,\ldots,\mu_{d-1})$ that satisfy the betweenness condition[^1]:

```math
\lambda_1 \geq \mu_1 \geq \lambda_2 \geq \mu_2
\geq \cdots
\geq \lambda_{d-1} \geq \mu_{d-1} \geq \lambda_d.
```

Each admissible $\mu$ therefore appears exactly once. For example, for the $U(3)$ highest weight $\lambda=(2,1,0)$, the admissible $U(2)$ labels $\mu=(\mu_1,\mu_2)$ satisfy

```math
2 \geq \mu_1 \geq 1,
\qquad
1 \geq \mu_2 \geq 0,
```

and hence are

```math
(2,1),\quad (2,0),\quad (1,1),\quad (1,0).
```

Iterating this interlacing rule along the chain $U(d)\supset U(d-1)\supset\cdots\supset U(1)$ produces triangular Gelfand-Tsetlin patterns. Each pattern records the successive highest weights, and the admissible patterns index the corresponding Gelfand-Tsetlin basis. Because the branching is multiplicity-free and the final $U(1)$ representations are one-dimensional, these patterns label the basis vectors uniquely, up to normalization and phase.

[^1]: I. M. Gelfand and M. L. Tsetlin, “Finite-dimensional representations of the group of unimodular matrices,” *Doklady Akademii Nauk SSSR* **71**, 825--828 (1950).

## Bra-ket notation

In quantum physics, the complex vectors -- here, of the carrier space -- are often interchangeably called **states**, and are introduced in the ket notation as $|\psi\rangle$. In this notation, the symbol $|\cdot \rangle$ denotes a vector and the symbols inside are labels (vector names, or for a collection of eigenvectors of a particular operator, indices). Similarly, a linear form is written as $\langle \cdot |$ and called a **bra**, so that, if by $\langle \psi |$ we denote a form dual to $|\psi\rangle$, $\langle \psi | \phi \rangle$, the **bra-ket**, is the Hermitean inner product $(\phi,\psi)$. Note that in the bra-ket notation, the second argument is linear, unlike the mathematician's notation.

In this context, the adjoint of an operator $O$ is denoted by $O^\dagger$. An operator is
sandwiched between a bra and a ket as $\langle \psi | O | \phi \rangle$, meaning the bra
$\langle \psi |$ applied to the vector $O |\phi\rangle$. This expression can be interpreted in the mathematicians' notation as $(O \phi, \psi)$.

As the Gelfand-Tsetlin  basis constructed above is orthonormal,
$\langle A | \rho_\lambda(u) | B \rangle$ is the entry of the matrix of
$\rho_\lambda(u)$ in row $A$ and column $B$.

## SU(2) case

The simplest example comes from SU(2), a subgroup of U(2) that differs from it only by phase. A U(2) element may have any determinant on the complex unit circle; multiplying the element by a suitable unimodular scalar produces an SU(2) element with unit determinant. The Lie algebra su(2) is spanned by $J_z, J_+, J_-$ subject to the commutation relations

```math
[J_z, J_\pm] = \pm J_\pm, \qquad [J_+, J_-] = 2 J_z.
```

In physical context, the eigenvectors of $J_z$ in the irreducible representation of dimension $2j+1$ are written as $|j, m\rangle$, with $m = -j, -j+1, \ldots, j$.

They satisfy $J_z |j,m\rangle = m |j,m\rangle$, while $J_\pm$ act as ladder operators. These vectors describe spin and, for $j\in\mathbb{N}$, the orbital angular momentum.

The same algebraic structure also describes two bosonic modes with creation operators $a^\dagger_1$ and $a^\dagger_2$. These operators are defined by their commutation relations $[a_i,a_i^\dagger]=1$ for $i=1,2$ and $[a_i,a_j]=[a_i,a_j^\dagger]=0$ for $i\neq j$. Now, define

```math
J_z = \tfrac{1}{2}(a^\dagger_1 a_1 - a^\dagger_2 a_2), \qquad J_+ = a^\dagger_1 a_2, \qquad J_- = a^\dagger_2 a_1.
```

These operators satisfy the commutation relations above. In the physical context, $J_z$ measures the occupation difference between the modes, while $J_\pm$ move particles between them. A state with $n_1$ particles in mode 1 and $n_2$ in mode 2 corresponds to the irreducible representation with

```math
j = \tfrac{1}{2}(n_1 + n_2), \qquad m = \tfrac{1}{2}(n_1 - n_2).
```

With these identifications, the vector $|n_1, n_2\rangle$, called the Fock state, is mathematically equivalent to the spin state $|j, m\rangle$.


For SU(2), a GT pattern is an inverted triangle of the form

```math
\left|\begin{array}{ccc}
2j&&0\\
&m+j&
\end{array}\right\rangle.
```

The top row specifies the total spin $j$, while the single entry in the bottom row determines $m$. The next section develops the general structure. For now, consider a simple example:

```@repl stateexample
using GroupFunctions

# Irrep [2,0] of SU(2): spin-1, dimension 3
λ = [2, 0];
basis = basis_states(λ);

for b in basis
    print(b)
end
```


## SU(N): Multiple $J^{(l)}_z$ operators


For SU(N), there are $N-1$ commuting operators $J^{(l)}_z$, each with corresponding ladder operators $J^{(l)}_\pm$. Each triple $J^{(l)}_z, J^{(l)}_\pm$ satisfies the su(2) commutation relations. See Alex et al.[^2] for explicit definitions. Because these operators and vectors describe many different physical settings, we retain the abstract viewpoint here.

A GT pattern labels a basis vector by recording the integer partitions obtained through successive group restrictions. This vector is also a simultaneous eigenvector of all $J^{(l)}_z$, whose eigenvalues determine the p-weight defined below. In general, the p-weight fixes the vector only when the irreducible representation is totally symmetric or totally antisymmetric. The pattern forms a triangular array:

```math
\left|\begin{array}{ccccccc}
m_{1,N} & & m_{2,N} & & \cdots & & m_{N,N} \\
& m_{1,N-1} & & \cdots & & m_{N-1,N-1} & \\
& & \ddots & & ~ & & \\
& & m_{1,2} & & m_{2,2} & & \\
& & & m_{1,1} & & &
\end{array}\right\rangle
```


The **top row** $(m_{1,N}, m_{2,N}, \ldots, m_{N,N})$ specifies the irreducible representation and corresponds to a Young diagram, or integer partition. The remaining rows identify a vector within that representation by satisfying the **betweenness condition**:

```math
m_{k, l+1} \geq m_{k, l} \geq m_{k+1, l+1},
```


```@repl stateexample

# Irrep [2,1,0] of SU(3): 8-dimensional
λ = [2, 1, 0];
basis = basis_states(λ);
println("Dimension: ", length(basis))

# Show first few patterns
for b in basis[1:3]
    display(b)
end
```

The **p-weight** of a pattern is the sequence of differences between consecutive row sums: if $\sigma_l = \sum_k m_{k,l}$, then $w_l = \sigma_l - \sigma_{l-1}$. For certain representations, the p-weight has a direct physical interpretation.

The *symmetric representation* $\lambda = [N, 0, \ldots, 0]$, whose Young diagram has a single row, describes $N$ indistinguishable bosons in $d$ modes. Reversing the p-weight of this representation gives the occupation numbers $(n_1, n_2, \ldots, n_d)$. Accordingly, the function `occupation_number` used below is defined as `reverse ∘ pweight`.

```@repl stateexample
# 2 photons in 3 modes
λ = [2, 0, 0];
basis = basis_states(λ);

for b in basis
    println(b, " → occupation ", occupation_number(b)) # will output all Fock states of 3 modes with total number of particles = 2
end
```

Each Fock state corresponds to exactly one GT pattern.


The antisymmetric representation $\lambda = [1, 1, \ldots, 1,0,0,\ldots,0]$, whose Young diagram has a single column, describes fermions. The Pauli exclusion principle restricts each occupation number to 0 or 1.

```@repl stateexample
# 2 fermions in 3 modes
λ = [1, 1, 0];
basis = basis_states(λ);
println("Dimension: ", length(basis))  # Should be 3

for b in basis
    println(b, " → occupation ", occupation_number(b)) # will output all fermionic states -- in Fock space, permutations of [1,1,0] 
end
```

Mixed irreducible representations such as $[2, 1]$ are neither fully symmetric nor fully antisymmetric. They arise physically when particles are *partially distinguishable*.

Consider three photons entering an interferometer. Two arrive simultaneously with identical spectral properties and can therefore interfere fully. Because the third arrives later, it cannot interfere with the earlier pair. The resulting quantum state has components in both the symmetric sector $[3,0,0]$ and the mixed-symmetry sector $[2,1,0]$. This decomposition appears in observable coincidence rates; see Section 3.2 of Amaro-Alcalá et al.[^3] for an explicit computation.

The simpler case of two photons with partial overlap likewise involves both $[2]$ (symmetric, contributing via the permanent) and $[1,1]$ (antisymmetric, contributing via the determinant); see Section 2.1 of the above paper.

[^2]: A. Alex, M. Kalus, A. Huckleberry, and J. von Delft, “A numerical algorithm for the explicit calculation of SU(N) and SL(N,C) Clebsch-Gordan coefficients,” *Journal of Mathematical Physics* **52**, 023507 (2011).

[^3]: D. Amaro-Alcalá, D. Spivak, and H. de Guise, “Sum rules in multiphoton coincidence rates,” *Physics Letters A* **384**, 126459 (2020).

```@repl stateexample
# Mixed irreducible representation [2,1,0]: dimension 8
λ = [2, 1, 0];
basis = basis_states(λ);

for b in basis
    println(occupation_number(b))
end
# Notice: p-weight (1,1,1) appears twice — inner multiplicity!
```

For mixed representations, multiple GT patterns can share the same p-weight. This inner multiplicity shows that occupation numbers alone do not uniquely specify a basis vector when particles have mixed exchange symmetry.


Once GT patterns have labeled the basis vectors, we can ask how their transition amplitudes depend on the unitary transformation $U$. For bosons, the answer involves permanents; for fermions, determinants; and for mixed symmetry, a generalization called the immanant.
