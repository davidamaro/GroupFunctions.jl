```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Immanants

Kostant's theorem ([background page](../background/immanants.md)) writes an immanant as a sum of diagonal group functions over the zero-`zweight` states of an irrep. The irrep `λ = [2,1,0]` has two zero-weight states, so the $(2,1)$ immanant of a $3\times3$ matrix is the sum of their two group functions. We pick the states out by their vanishing `zweight`, evaluate the immanant from its definition, and compare:

```@raw html
<p style="text-align: center;">
  <img style="width: 63%; height: auto;" src="../../assets/tutorial/tutorial_immanant_zero_weight.svg" alt="Kostant's zero-weight diagonal sum for the (2,1) immanant">
</p>
```

```@repl immanants_zero
using GroupFunctions

imm21(A) = 2*A[1,1]*A[2,2]*A[3,3] -
           A[1,2]*A[2,3]*A[3,1] -
           A[1,3]*A[2,1]*A[3,2];

λ = [2, 1, 0];
zero = filter(gt -> all(iszero, zweight(gt)), basis_states(λ));

M1 = su2_block(3,1,(0.5,0.9,0.2)) * su2_block(3,2,(0.4,0.2,0.6)) * su2_block(3,1,(0.5,0.7,0.3)); #certified to be random, picked by me

total = sum(group_function(λ, t, t, M1) for t in zero)
total ≈ imm21(M1)
```

## Permanent as special case

The permanent is the immanant of the fully symmetric irrep, $\mathrm{Per}=\mathrm{Imm}^{(n)}$, and there the zero-weight sum collapses to a single term. For three photons in four modes (`λ = [3,0,0,0]`) the group function reproduces the permanent of the relevant submatrix (with a $\sqrt{2}$ normalization factor -- see [the background page](../background/group_functions.md#Bosons:-the-permanent)):

```@raw html
<p style="text-align: center;">
  <img style="width: 63%; height: auto;" src="../../assets/tutorial/tutorial_permanent_submatrix.svg" alt="Repeated input occupations select repeated submatrix indices in the permanent">
</p>
```

```@repl immanants_perm
using GroupFunctions

permanent3(A) = A[1,1]*A[2,2]*A[3,3] +
                A[1,1]*A[2,3]*A[3,2] +
                A[1,2]*A[2,1]*A[3,3] +
                A[1,2]*A[2,3]*A[3,1] +
                A[1,3]*A[2,1]*A[3,2] +
                A[1,3]*A[2,2]*A[3,1];

U = su2_block(4, 1, (0.3, 0.7, 0.2)) * su2_block(4, 2, (0.5, 0.9, 0.5)) *
    su2_block(4, 1, (0.4, 0.8, 0.1)) * su2_block(4, 3, (0.6, 0.2, 0.6)) *
    su2_block(4, 2, (0.1, 0.3, 0.1)) * su2_block(4, 1, (0.9, 0.4, 0.7));

basis = basis_states([3,0,0,0]);
state_x = filter(s -> occupation_number(s) == [1,1,1,0], basis)[1]; #modes 1,2,3
state_y = filter(s -> occupation_number(s) == [0,2,0,1], basis)[1];#modes 2 (twice), 4

M2 = U[[1,2,3], [2,2,4]]; #the modes from above form the indices to extract matrix elements
group_function([3,0,0,0], state_x, state_y, U) ≈ permanent3(M2) / sqrt(2)
```
