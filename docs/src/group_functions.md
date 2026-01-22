# Group Functions

## The problem: mode mixing in quantum mechanics

Consider a system of $n$ modes, each described by a creation operator $a^\dagger_i$. A Fock state is built by applying creation operators to the vacuum:

$$|m_1, m_2, \ldots, m_n\rangle = \frac{(a_1^\dagger)^{m_1} \cdots (a_n^\dagger)^{m_n}}{\sqrt{m_1! \cdots m_n!}} |0\rangle.$$

Now suppose the modes get mixed by a unitary transformation $U \in \mathrm{U}(n)$:
$$a_i^\dagger \mapsto \sum_j U_{ji} a_j^\dagger.$$

What is the output state? Each creation operator in the original Fock state transforms according to the rule above. Expanding the product, we obtain a superposition of Fock states with coefficients that are polynomials in the matrix elements $U_{ij}$.

The function `group_function` computes these coefficients — the transition amplitudes $\langle m' | U | m \rangle$. The presentation here is slightly more mathematical than a typical quantum optics treatment, for consistency with representation theory literature and to handle more general cases beyond bosons.

## Bosons: the permanent

For bosonic systems, consider the transition amplitude between an input Fock state $|m\rangle$ and an output Fock state $|m'\rangle$, both having $N$ total particles. We follow the derivation of [Scheel (2004)](https://arxiv.org/abs/quant-ph/0406127).

 Label the $N$ particles by $\alpha = 1, \ldots, N$, assigning each to its input mode via $\alpha \mapsto i_\alpha$ such that mode $i$ appears $m_i$ times (so, e.g. for $\ket{2,1,0}$, $(i_1, i_2, i_3)=(1,1,2)$). Then:
 
$$U|m\rangle = \frac{1}{\sqrt{\prod_i m_i!}} \prod_{i=1}^{n} \left(\sum_j U_{ij} a_j^\dagger\right)^{m_i} |0\rangle = \frac{1}{\sqrt{\prod_i m_i!}} \sum_{j_1, \ldots, j_N} \left(\prod_{\alpha=1}^{N} U_{i_\alpha, j_\alpha}\right) \prod_{\alpha=1}^N a_{j_\alpha}^\dagger |0\rangle$$

The first step unfolds powers into labeled factors; the second is distributivity. 

So now we have products of operators acting on the vacuum $\ket{0}$; we wish to close it with a bra $\bra{m'}$ to calculate the matrix element. Let us act with the creatio operators on the vacuum and identify terms proportional to $\ket{m'}$ there; then, the matrix element will be a sum of prefactors leading to $\ket{m'}$.
 For terms where the tuple $(j_1, \ldots, j_N)$ contains mode $k$ exactly $m'_k$ times, the creation operators produce:

$$\prod_{\alpha=1}^N a_{j_\alpha}^\dagger |0\rangle = \sqrt{\prod_j m'_j!}\, |m'\rangle$$

The $\sqrt{m'_j!}$ arises because $m'_j$ identical creation operators acting on vacuum give $(a_j^\dagger)^{m'_j}|0\rangle = \sqrt{m'_j!}|m'_j\rangle_j$. Denoting the constraint on $(j_1, \ldots, j_N)$ as $\sim m'$:

$$\langle m' | U | m \rangle = \frac{\sqrt{\prod_j m'_j!}}{\sqrt{\prod_i m_i!}} \sum_{(j_1,\ldots,j_N) \sim m'} \prod_{\alpha=1}^{N} U_{i_\alpha, j_\alpha}$$

The above expression can be shown to be expressible as a permanent of a properly constructed matrix. This is, by the way, the basis for boson sampling problems; computing the permanent is #P-hard ([Aaronson & Arkhipov, 2011](https://arxiv.org/abs/1011.3245)), making bosonic transition amplitudes classically intractable.

The relation is as follows: construct an $N \times N$ matrix $M$ using output labeling $\beta \mapsto j_\beta$ (analogous to input): $M_{\alpha\beta} = U_{i_\alpha, j_\beta}$. The permanent sums over all permutations $\sigma \in S_N$:

$$\mathrm{perm}(M) = \sum_{\sigma \in S_N} \prod_{\alpha=1}^{N} U_{i_\alpha, j_{\sigma(\alpha)}}$$

Each valid $(j_1, \ldots, j_N) \sim m'$ corresponds to $\prod_j m'_j!$ permutations (permuting indices within each output mode). Therefore:
$$\sum_{(j_1,\ldots,j_N) \sim m'} \prod_{\alpha} U_{i_\alpha, j_\alpha} = \frac{\mathrm{perm}(M)}{\prod_j m'_j!}$$

As a result, we have

$$\langle m' | U | m \rangle = \frac{\mathrm{perm}(M)}{\sqrt{\prod_i m_i!}\sqrt{\prod_j m'_j!}}$$

**Example:** 3 modes, input $|m\rangle = |2,0,1\rangle$, output $|m'\rangle = |1,1,1\rangle$.

Input labeling ($N=3$ particles): $(i_1, i_2, i_3) = (1, 1, 3)$ — two particles from mode 1, one from mode 3.

Output labeling: $(j_1, j_2, j_3) = (1, 2, 3)$ — one particle into each mode.

The matrix $M$ has entries $M_{\alpha\beta} = U_{i_\alpha, j_\beta}$:

$$M = \begin{pmatrix} U_{11} & U_{12} & U_{13} \\ U_{11} & U_{12} & U_{13} \\ U_{31} & U_{32} & U_{33} \end{pmatrix}$$

Note the repeated rows from $m_1 = 2$.

In representation-theoretic terms, bosonic Fock states live in the symmetric subspace, the irrep $\lambda = [N, 0, \ldots, 0]$.

## Fermions: the determinant

For fermions, the derivation parallels the bosonic case with three modifications: occupation numbers are restricted to $m_i, m'_j \in \{0,1\}$ (Pauli exclusion), leading to all normalization factors become $\sqrt{0!}=\sqrt{1!} = 1$, and anticommutation introduces signs.

The expansion step has the aame structure as bosons:

$$U|m\rangle = \sum_{j_1, \ldots, j_N} \left(\prod_{\alpha=1}^{N} U_{i_\alpha, j_\alpha}\right) a_{j_1}^\dagger \cdots a_{j_N}^\dagger |0\rangle$$

To project the result to final bra $\bra{m'}$, consider the following. For a term to contribute to $|m'\rangle$, the tuple $(j_1, \ldots, j_N)$ must have all distinct entries (otherwise $a_j^\dagger a_j^\dagger = 0$), forming a permutation of the occupied output modes. Write $j_\alpha = j'_{\sigma(\alpha)}$ where $(j'_1, \ldots, j'_N)$ is the sorted list and $\sigma \in S_N$. Reordering to standard form introduces signs from the anticommutation relation (rectifying every transposition changes it: $a_2^\dagger a_1^\dagger=-a_1^\dagger a_2^\dagger$):

$$a_{j_1}^\dagger \cdots a_{j_N}^\dagger |0\rangle = \mathrm{sgn}(\sigma) \cdot a_{j'_1}^\dagger \cdots a_{j'_N}^\dagger |0\rangle = \mathrm{sgn}(\sigma) |m'\rangle$$

As a result,

$$\langle m' | U | m \rangle = \sum_{\sigma \in S_N} \mathrm{sgn}(\sigma) \prod_{\alpha=1}^{N} U_{i_\alpha, j'_{\sigma(\alpha)}} = \det(M)$$


where $M_{\alpha\beta} = U_{i_\alpha, j'_\beta}$ is the submatrix of $U$ with rows = occupied input modes, columns = occupied output modes.

Unlike permanents, determinants can be computed efficiently in $O(N^3)$ time, which underlies the tractability of free-fermion systems.

## The general formula

Both results share a common structure: a sum over permutations, weighted by representation-dependent coefficients, times a monomial in matrix elements. For bosons:

$$\langle m' | U | m \rangle \propto \sum_{\sigma \in S_N} 1 \cdot \prod_{\alpha} U_{i_\alpha, j_{\sigma(\alpha)}}$$

and for fermions:

$$\langle m' | U | m \rangle \propto \sum_{\sigma \in S_N} \mathrm{sgn}(\sigma) \cdot \prod_{\alpha} U_{i_\alpha, j_{\sigma(\alpha)}}$$

The weights $1$ and $\mathrm{sgn}(\sigma)$ are the matrix elements of the trivial and sign representations of $S_N$, both one-dimensional.

For a general irrep $\lambda$, the representation has dimension $f^\lambda > 1$, and the weight becomes a matrix element $\omega^\lambda_{i,j}(\sigma)$ of the Young orthogonal representation. The general formula, due to [Grabmeier & Kerber (1985)](https://doi.org/10.1007/BF00046717), is:

$$T^\lambda_{U,V} = \frac{1}{\sqrt{\Theta^\lambda_U \Theta^\lambda_V}} \sum_{\gamma} \left( \sum_{\sigma \in S_\alpha \gamma S_\beta} \omega^\lambda_{i,j}(\sigma) \right) X_{f \circ \gamma, g}$$

Here $U, V$ are semistandard tableaux labeling basis states (equivalent to GT patterns), $\Theta$ is a normalization factor generalizing $\sqrt{m!}$, and $X_{f \circ \gamma, g}$ is the monomial $\prod_k U_{i_k, j_k}$. The double coset decomposition $S_\alpha \backslash S_N / S_\beta$ groups permutations contributing the same monomial.

The partition $\lambda$ labels the symmetry type: $[N]$ for bosons, $[1^N]$ for fermions, mixed shapes for particles with mixed exchange symmetry. The function `group_function(λ, ...)` evaluates this formula.
