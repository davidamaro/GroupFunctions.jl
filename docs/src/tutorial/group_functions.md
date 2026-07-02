```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Group functions

A group function is a matrix element of an irrep: given an initial and a final
basis state and a unitary `U`, it returns the amplitude between them. We reuse the
states from the [previous page](states.md) — the symmetric irrep `[2, 0, 0]`, two
photons in three modes. For the theory (permanents, determinants, the general
formula) see the [background page](../background/group_functions.md).

```@repl gf
using GroupFunctions

λ = [2, 0, 0];
basis = basis_states(λ);

initial = basis[findfirst(gt -> occupation_number(gt) == [2, 0, 0], basis)]
final   = basis[findfirst(gt -> occupation_number(gt) == [1, 1, 0], basis)]
```

## One numeric entry

A single matrix element is one call. Here `U` is a 50:50 beam splitter on modes
1–2 (an SU(2) block embedded in three modes), and we ask for the amplitude
$\langle 1,1,0 \mid U \mid 2,0,0\rangle$:

```@raw html
<p style="text-align: center;">
  <img style="width: 63%; height: auto;" src="../../assets/tutorial/tutorial_group_function.svg" alt="A group function as one input-to-output transition amplitude">
</p>
```

```@repl gf
U = su2_block(3, 1, (0.0, pi/2, 0.0));
amp = group_function(λ, final, initial, U)
```

The result is a complex amplitude; its squared modulus is the transition
probability:

```@repl gf
abs2(amp)
```

## The whole representation at once

Passing only `λ` and `U` — no states — returns *every* matrix element at once,
together with the patterns indexing its rows and columns. To keep the printed
matrix small we drop here to the two-mode irrep `[2, 0]`, whose representation is
$3\times3$:

```@repl gf_matrix
using GroupFunctions

BS = su2_block(2, 1, (0.0, pi/2, 0.0));
values, patterns = group_function([2, 0], BS);

size(values)
round.(values, digits=4)
```

## One-dimensional irreps

The one-dimensional irreps of U(d) correspond to the constant partitions: betweenness
forces every entry of the pattern to equal the (constant) top row, so there is
exactly one basis vector. The zero partition (`[0,0]`, `[0,0,0]`) is the trivial
irrep, where every group element acts as `1`; the all-ones partition
(`[1,1,1]`) is the determinant rep, acting as $\det U$. The all-entries call
returns that single scalar as its first result:

```@repl gf_onedim
using GroupFunctions
using LinearAlgebra: qr, det
U = Matrix(qr(rand(ComplexF64, 3, 3)).Q);

# trivial irrep: scalar 1, and a single basis state
group_function([0, 0, 0],U)[1]
length(group_function([0, 0,0],U)[2])

# determinant irrep: returns det(U)
group_function([1, 1, 1], U)[1]
det(U)
```

## A symbolic entry

Omit `U` altogether and `group_function` returns the matrix element
*symbolically*, as a function of $U$ with matrix elements `u_i_j` rather than a
number. The ability to compute group functions purely symbolically is this library
most distinctive capability. The same machinery that
gives permanents and determinants for the symmetric and antisymmetric irreps
produces, in general, the [immanants](immanants.md) of the relevant submatrix.

```@repl gf
group_function(λ, final, initial)
```
