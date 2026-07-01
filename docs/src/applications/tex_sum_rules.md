Sum rules

The coincidence rate for indistinguishable photons is the modulus squared of a permanent built from the scattering matrix U, an element of the d-dimensional unitary group. Permanents are expensive to compute. However, nontrivial constraints apply to sums of rates, and in some cases the sum remains unchanged if U is replaced by a simpler matrix. This can speed up rate-sum computations. This is the result of Amaro-Alcala, Spivak, and de Guise (2020). We reproduce the three-photon case here, following Section 3 of the article.

Rate sum invariance

Three photons enter a four-mode interferometer in the occupation-number ket with entries zero, one, one, and one. We sum the rates over a set of output states, comparing the scattering matrices U1 and U2. The second matrix has an extra block, called addon in the code, that mixes modes 1, 2, and 3 while leaving mode 4 unchanged.

Because the additional block only reshuffles photons within modes 1, 2, and 3, the probability of observing n photons in those modes is the same whether the block is present or absent. The block cannot change the photon number in this subspace. This can be verified directly by evaluating the corresponding rate sums.

The summed-over outputs are the six states with one photon in mode 4 and two photons shared among modes 1, 2, and 3. When we sum over this complete set, the block makes no difference, as shown in Section 3 of the paper.

The leftmost block can therefore be removed without changing the summed rate. In an interferometer, this means the optical elements can be deleted, as shown in Figure 3 of the paper.

Faster evaluation

The point of removing the block is that the simpler matrix has zeros, and zeros make the permanent cheap. After dropping the block that acts on modes 1, 2, and 3, the matrix U1 has a staircase of zeros above the diagonal. Every entry more than one step above the diagonal vanishes. A matrix of this shape, with vanishing superdiagonals except for the one directly above the main diagonal, is called upper Hessenberg.

In the three-by-three submatrix used in the example, the top-right entry is numerically zero. In this case, that entry is the entire second superdiagonal, so the submatrix is upper Hessenberg.

For such a matrix, the permanent, normally an expensive sum over all permutations, equals an ordinary determinant, provided we first flip the sign of the entries just above the diagonal. The example builds this sign-flipped matrix and compares its determinant with the permanent of the original matrix.

The two agree. This allows a faster computation of the summed rate. Instead of evaluating permanents after the unimportant block has been removed, we can use determinants of a modified matrix. Determinants are much more computationally efficient: standard algorithms scale like order n cubed, whereas direct permanent computation scales like order two to the n times n squared. In a three-by-three matrix, this does not make much difference, but the method shown here and in the paper generalises to more complex situations.
