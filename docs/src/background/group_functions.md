```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Group Functions

The [basis states and GT patterns page](states.md) introduced these labels for basis states, with different irreps corresponding to different exchange symmetries. Here we sketch the mathematical framework behind the library; actual code is shown in the applications pages on [HOM effect](../applications/quantum_optics.md), [qubit transmission with entangled light](../applications/qubit_transmission.md), [sum rules](../applications/sum_rules.md), as well as the related notes on [characters](characters.md) and [immanants](immanants.md).

## The problem: mode mixing in quantum mechanics

Consider a system of $n$ modes, each described by a creation operator $a^\dagger_i$. A Fock state is built by applying creation operators to the vacuum:

$$|m_1, m_2, \ldots, m_n\rangle = \frac{(a_1^\dagger)^{m_1} \cdots (a_n^\dagger)^{m_n}}{\sqrt{m_1! \cdots m_n!}} |0\rangle.$$

Now suppose the modes get mixed by a unitary transformation $U \in \mathrm{U}(n)$:
$$a_i^\dagger \mapsto \sum_j U_{ji} a_j^\dagger.$$

What is the output state? Each creation operator in the original Fock state transforms according to the rule above. Expanding the product, we obtain a superposition of Fock states with coefficients that are polynomials in the matrix elements $U_{ij}$.

The function `group_function` computes these coefficients — the transition amplitudes $\langle m' | D^{(\lambda)}(U) | m \rangle$, where $\lambda$ labels the space the operator is acting on (the irreducible representation). The presentation here is slightly more mathematical than a typical quantum optics treatment, for consistency with representation theory literature and to handle more general cases beyond bosons. The partition $\lambda$ labels the symmetry type: $[N,0,0,\ldots]$ for bosons, $[1,1,\ldots,1]$ for fermions, mixed shapes for particles with mixed exchange symmetry. 
In particular, we keep the label $\lambda$ in the equations to accommodate the more general cases.

## Bosons: symmetric representation

For bosonic systems, consider the transition amplitude between an input Fock state $|m\rangle$ and an output Fock state $|m'\rangle$, both having $N$ total particles. We follow Scheel's derivation.[^1]

 Label the $N$ particles by $\alpha = 1, \ldots, N$, assigning each to its input mode via $\alpha \mapsto i_\alpha$ such that mode $i$ appears $m_i$ times (so, e.g. for $\ket{2,1,0}$, $(i_1, i_2, i_3)=(1,1,2)$). Then:
 
$$D^{(\lambda)}(U)|m\rangle = \frac{1}{\sqrt{\prod_i m_i!}} \prod_{i=1}^{n} \left(\sum_j U_{ij} a_j^\dagger\right)^{m_i} |0\rangle = \frac{1}{\sqrt{\prod_i m_i!}} \sum_{j_1, \ldots, j_N} \left(\prod_{\alpha=1}^{N} U_{i_\alpha, j_\alpha}\right) \prod_{\alpha=1}^N a_{j_\alpha}^\dagger |0\rangle$$

The first step unfolds powers into labeled factors; the second is distributivity. 

So now we have products of operators acting on the vacuum $\ket{0}$; we wish to close it with a bra $\bra{m'}$ to compute the matrix element. Let us act with the creation operators on the vacuum and identify terms proportional to $\ket{m'}$ there; then, the matrix element will be a sum of prefactors leading to $\ket{m'}$.
 For terms where the tuple $(j_1, \ldots, j_N)$ contains mode $k$ exactly $m'_k$ times, the creation operators produce:

$$\prod_{\alpha=1}^N a_{j_\alpha}^\dagger |0\rangle = \sqrt{\prod_j m'_j!}\, |m'\rangle$$

The $\sqrt{m'_j!}$ arises because $m'_j$ identical creation operators acting on vacuum give $(a_j^\dagger)^{m'_j}|0\rangle = \sqrt{m'_j!}|m'_j\rangle_j$. Denoting the constraint on $(j_1, \ldots, j_N)$ as $\sim m'$:

$$\langle m' | D^{(\lambda)}(U) | m \rangle = \frac{\sqrt{\prod_j m'_j!}}{\sqrt{\prod_i m_i!}} \sum_{(j_1,\ldots,j_N) \sim m'} \prod_{\alpha=1}^{N} U_{i_\alpha, j_\alpha}$$

The above expression can be shown to be expressible as a permanent of a properly constructed matrix. This is, by the way, the basis for boson sampling problems; computing the permanent is #P-hard,[^2] making bosonic transition amplitudes classically intractable.

The relation is as follows: construct an $N \times N$ matrix $M$ using output labeling $\beta \mapsto j_\beta$ (analogous to input): $M_{\alpha\beta} = U_{i_\alpha, j_\beta}$. The permanent sums over all permutations $\sigma \in S_N$:

$$\mathrm{perm}(M) = \sum_{\sigma \in S_N} \prod_{\alpha=1}^{N} U_{i_\alpha, j_{\sigma(\alpha)}}$$

Each valid $(j_1, \ldots, j_N) \sim m'$ corresponds to $\prod_j m'_j!$ permutations (permuting indices within each output mode). Therefore:
$$\sum_{(j_1,\ldots,j_N) \sim m'} \prod_{\alpha} U_{i_\alpha, j_\alpha} = \frac{\mathrm{perm}(M)}{\prod_j m'_j!}$$

As a result, we have

$$\langle m' | D^{(\lambda)}(U) | m \rangle = \frac{\mathrm{perm}(M)}{\sqrt{\prod_i m_i!}\sqrt{\prod_j m'_j!}}$$

**Example:** 3 modes, input $|m\rangle = |2,0,1\rangle$, output $|m'\rangle = |1,1,1\rangle$.

Input labeling ($N=3$ particles): $(i_1, i_2, i_3) = (1, 1, 3)$ — two particles from mode 1, one from mode 3.

Output labeling: $(j_1, j_2, j_3) = (1, 2, 3)$ — one particle into each mode.

The matrix $M$ has entries $M_{\alpha\beta} = U_{i_\alpha, j_\beta}$:

$$M = \begin{pmatrix} U_{11} & U_{12} & U_{13} \\ U_{11} & U_{12} & U_{13} \\ U_{31} & U_{32} & U_{33} \end{pmatrix}$$

Note the repeated rows from $m_1 = 2$.

In representation-theoretic terms, bosonic Fock states live in the symmetric subspace, the irrep $\lambda = [N, 0, \ldots, 0]$.

## Fermions: anti-symmetric representation

For fermions, the derivation parallels the bosonic case with three modifications: occupation numbers are restricted to $m_i, m'_j \in \{0,1\}$ (Pauli exclusion), so all normalization factors become $\sqrt{0!}=\sqrt{1!} = 1$, and anticommutation introduces signs.

The expansion step has the same structure as bosons:

$$D^{(\lambda)}(U) |m\rangle = \sum_{j_1, \ldots, j_N} \left(\prod_{\alpha=1}^{N} U_{i_\alpha, j_\alpha}\right) a_{j_1}^\dagger \cdots a_{j_N}^\dagger |0\rangle$$

To project the result to final bra $\bra{m'}$, consider the following. For a term to contribute to $|m'\rangle$, the tuple $(j_1, \ldots, j_N)$ must have all distinct entries (otherwise $a_j^\dagger a_j^\dagger = 0$), forming a permutation of the occupied output modes. Write $j_\alpha = j'_{\sigma(\alpha)}$ where $(j'_1, \ldots, j'_N)$ is the sorted list and $\sigma \in S_N$. Reordering to standard form introduces signs from the anticommutation relation (rectifying every transposition changes it: $a_2^\dagger a_1^\dagger=-a_1^\dagger a_2^\dagger$):

$$a_{j_1}^\dagger \cdots a_{j_N}^\dagger |0\rangle = \mathrm{sgn}(\sigma) \cdot a_{j'_1}^\dagger \cdots a_{j'_N}^\dagger |0\rangle = \mathrm{sgn}(\sigma) |m'\rangle$$

As a result,

$$\langle m' | D^{(\lambda)}(U) | m \rangle = \sum_{\sigma \in S_N} \mathrm{sgn}(\sigma) \prod_{\alpha=1}^{N} U_{i_\alpha, j'_{\sigma(\alpha)}} = \det(M)$$


where $M_{\alpha\beta} = U_{i_\alpha, j'_\beta}$ is the submatrix of $U$ with rows = occupied input modes, columns = occupied output modes.

Unlike permanents, determinants can be computed efficiently in $O(N^3)$ time, which underlies the tractability of free-fermion systems.

## The general formula

Both results share a common structure: a sum over permutations, weighted by representation-dependent coefficients, times a monomial in matrix elements. For bosons:

$$\langle m' | D^{(\lambda)}(U) | m \rangle \propto \sum_{\sigma \in S_N} 1 \cdot \prod_{\alpha} U_{i_\alpha, j_{\sigma(\alpha)}}$$

and for fermions:

$$\langle m' | D^{(\lambda)}(U) | m \rangle \propto \sum_{\sigma \in S_N} \mathrm{sgn}(\sigma) \cdot \prod_{\alpha} U_{i_\alpha, j_{\sigma(\alpha)}}$$

The weights $1$ and $\mathrm{sgn}(\sigma)$ are the matrix elements of the trivial and sign representations of $S_N$, both one-dimensional.

For a general irrep $\lambda$, the picture complicates a bit: the basis vectors can no longer be uniquely indexed by occupation numbers, and semistandard Young tableaux are used instead (equivalent to [GT patterns](states.md)), here denoted by $A$ and $B$. A semistandard tableau is a filling of the Young diagram with numbers nondecreasing along each row and increasing down each column. The end formula takes a similar form:


$$\langle A \vert D^{(\lambda)}(U) \vert B\rangle = \sum (\text{coefficient}) (\text{monomial in matrix elements of }U),$$

but the construction of coefficients and monomials requires explanation of the notation. First, let us consider the basis states: the user-facing functions mostly use Gelfand-Tsetlin patterns; which are translated in the code to semistandard Young tableaux. For an example of such a tableau, let us take $A$ in the irrep $\lambda=[4,1]$ to be

```@raw html
  <div style="text-align:center;">
    <table style="border-collapse:collapse; margin-top:0.4em;">
      <tr>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">1</td>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">2</td>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">3</td>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">4</td>
      </tr>
      <tr>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">2</td>
      </tr>
    </table>
  </div>
```

The monomials are still constructed as products of matrix elements of $U$, but the indexing depends nontrivially on $A$ and $B$, and the expression involves the Young orthogonal representation of the permutation group $S_N$ too, where $N$ is the number of particles ($N=5$ in the example above). *Standard* Young tableaux enumerate the basis vectors of this representation: fillings that increase both along rows and down columns, using each of $1,\ldots,N$ exactly once. Each semistandard $A$ has an associated standard $\bar A$; for the $A$ above, the associated $\bar A$ is



```@raw html
  <div style="text-align:center;">
    <table style="border-collapse:collapse; margin-top:0.4em;">
      <tr>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">1</td>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">2</td>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">4</td>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">5</td>
      </tr>
      <tr>
        <td style="border:1px solid currentColor; width:2em; height:2em; text-align:center;">3</td>
      </tr>
    </table>
  </div>
```

The monomial is built using the fills of these tableaux. A fill is the sequence of entries read as concatenated rows; for the $A$ above,


$$\phi_A(k) \coloneqq k\text{-th element of }(1,2,3,4,2),$$

and similarly the fill $\phi_{\bar A}(k)$ of $\bar A$ is the $k$-th element of $(1,2,4,5,3)$; fills of $B$ and $\bar B$ are denoted in the same way. The fill of a standard tableau is always a permutation of $1,\ldots,N$, and in particular has an inverse. Define the compositions $f \coloneqq \phi_A \circ \phi_{\bar A}^{-1}$ and $g \coloneqq \phi_B \circ \phi_{\bar B}^{-1}$. The monomial then takes the form $\prod_k U_{f \circ \gamma^{-1}(k),\, g(k)}$, summed over selected permutations $\gamma$ with proper coefficients.



We now have the pieces needed to state the general formula, due to Grabmeier and Kerber:[^3]


$$\langle A \vert D^{(\lambda)}(U) \vert B \rangle = \frac{1}{\sqrt{\Theta^\lambda_A \Theta^\lambda_B}} \sum_{\gamma} \left( \sum_{\sigma \in S_A \gamma S_B} \omega^\lambda_{\bar A,\bar B}(\sigma) \right) \prod_{k} U_{f \circ \gamma^{-1}(k),\; g(k)}.$$


Some new symbols appear in the above formula, and they require careful explanation. 

- The $\Theta$ symbols are normalization factors generalizing $\sqrt{m!}$.
- The outer sum runs over *representatives of the double coset decomposition* $S_A \backslash S_N / S_B$. Here, $S_A$ is the stabilizer of $A$; these are the permutations which do not change the fills. For the semistandard Young tableau $A$ above, the only repeated entry is $2$ at positions $2$ and $5$ of the fill and the stabilizer group is generated by the transposition $(2,5)$ alone. The double coset structure then partitions $S_N$: for any $\gamma, \gamma' \in S_N$, the sets $S_A \gamma S_B$ and $S_A \gamma' S_B$ are either equal or disjoint. Hence, each double coset has well-defined representatives; and the sum runs over these; $\gamma$ denotes the representative.
- The inner sum runs over the coset $S_A \gamma S_B$ defined by $\gamma$, with the summand $\omega^\lambda_{\bar A,\bar B}(\sigma)$ being the matrix element of the Young orthogonal representation between the basis vectors labelled by $\bar A$ and $\bar B$.


The main functionality of this library – the function `group_function(λ, ...)` – evaluates this formula. For the tableau $A$ above (as a GT pattern) and a symbolic $U$, the diagonal matrix element is $\langle A\vert D^{(\lambda)}(U) \vert A\rangle$:

```@repl bggroupfun
using GroupFunctions
# The GT pattern equivalent to the semistandard tableau A above
A = GTPattern([[4, 1, 0, 0, 0], [4, 1, 0, 0], [3, 1, 0], [2, 1], [1]]);
group_function([4, 1, 0, 0, 0], A, A)
```

[^1]: S. Scheel, “Permanents in linear optical networks,” *arXiv preprint quant-ph/0406127* (2004).

[^2]: S. Aaronson and A. Arkhipov, “The computational complexity of linear optics,” in *Proceedings of the 43rd Annual ACM Symposium on Theory of Computing*, 333--342 (2011).

[^3]: J. Grabmeier and A. Kerber, “The evaluation of irreducible polynomial representations of the general linear groups and of the unitary groups over fields of characteristic 0,” *Acta Applicandae Mathematicae* **8**, 271--291 (1987).
