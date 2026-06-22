```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Sum rules

The coincidence rate for indistinguishable photons is the modulus squared of a
permanent (see [the background page](../background/group_functions.md)) built from the scattering matrix $U$ (element of U(d) for $d$ modes). Permanents are expensive to compute; however, there exist nontrivial constraints on *sums* of rates, and in some cases the rate sum is unchanged if $U$ is replaced by a "simpler" one. This allows for quicker calculation of sum rates.  This is the result of
[Amaro-Alcalá, Spivak & de Guise (2020)](https://doi.org/10.1016/j.physleta.2020.126459); we
reproduce its three-photon case here (see section 3 of the article).

## Rate sum invariance

Three photons enter a four-mode interferometer, in the state $\ket{0,1,1,1}$ (written as a ket with occupation numbers). We sum the rates over a set of
output states, comparing the scattering matrix `U1` against `U2` with an extra block `addon` tacked on, mixing modes (1,2,3) while leaving
mode 4 alone.

Since the additional block only reshuffles photons within modes (1,2,3), the probability of observing $n$ photons there is the same, whether it is present or not: it can not change the photon number in this subspace. Let us verify this with code.

The summed-over outputs are the six states with one photon in
mode 4 and two photons shared among modes 1,2,3. Summing over a complete set, the block makes no difference (see paper's section 3):

```@repl sumrules
using GroupFunctions

U1 =  su2_block(4,3,(0.3,0.7,0.3)) *
       su2_block(4,2,(0.4,0.8,0.4))*su2_block(4,1,(0.5,0.9,0.2)); 
addon=su2_block(4,1,(0.6,0.5,0.7)) * su2_block(4,2,(.1,.2,.3))*su2_block(4,1,(.21,.37,.91));
 
U2 = addon*U1; 
basis = basis_states([3,0,0,0]);
input = filter(s -> occupation_number(s) == [0,1,1,1], basis)[1];
outs  = filter(s -> occupation_number(s)[4] == 1, basis); 
#↑ occupation numbers like [1,1,0,1], [2,0,0,1]

ratesum(W) = sum(abs2(group_function([3,0,0,0], o, input, W)) for o in outs);

ratesum(U1)
ratesum(U2)
ratesum(U1) ≈ ratesum(U2)
```

The leftmost block can be removed (which corresponds to deletion of optical elements in an interferometer, see paper Fig. 3), without changing the summed rate.

## Faster evaluation

The point of removing the block is that the simpler matrix has zeros, and zeros
make the permanent cheap. Look at `U1`: dropping the mode-(1,2,3) block left it
with a staircase of zeros above the diagonal — every entry more than one step
above the diagonal vanishes. A matrix of this shape (with vanishing superdiagonals except the one directly above main diagonal) is called *upper Hessenberg*:

```@repl sumrules
M = U1[[2,3,4],[2,3,4]];
M
```

The top-right entry `M[1,3]` is (numerically) zero (in this rather simple case, this entry is the entire superdiagonal); therefore, the matrix is upper Hessenberg.

For such a matrix the permanent, normally an expensive sum over all
permutations, equals an ordinary determinant, provided we first flip the sign
of the entries just above the diagonal. Build that sign-flipped matrix `T` and
compare:

```@repl sumrules
using LinearAlgebra: det

permanent3(A) = A[1,1]*A[2,2]*A[3,3] + A[1,1]*A[2,3]*A[3,2] +
                A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] +
                A[1,3]*A[2,1]*A[3,2] + A[1,3]*A[2,2]*A[3,1];

T = copy(M);
T[1,2] = -T[1,2];     # negate the entries one step above the diagonal
T[2,3] = -T[2,3];

permanent3(M)
det(T)
permanent3(M) ≈ det(T)
```

The two agree! This allows for a faster calculation of the sum rate: instead of permanents, after the removal of some unimportant blocks (`addon`), we can use determinants (of a modified matrix) instead. And determinants are *much more computationally efficient to calculate* ($O(n^3)$ vs $O(2^n n^2)$ for permanents) -- in a $3\times3$ matrix it is not much difference, but the method shown here (and in the paper) generalizes to more complex situations.
