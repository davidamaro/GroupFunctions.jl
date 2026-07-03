```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```

# Basic quantum optics using GroupFunctions.jl

## Introduction

Linear optical networks, such as beamsplitters and phase shifters, act on $n$ optical modes via SU($n$) transformations. A key quantity is the transition amplitude between Fock states under these transformations.

Consider $N$ photons distributed across $n$ modes. The Fock state $|m_1, m_2, \ldots, m_n\rangle$ with $\sum_i m_i = N$ transforms under passive linear optics as:
$$|m_1, \ldots, m_n\rangle \xrightarrow{U} \sum_{m'} c_{m,m'}(U) |m'_1, \ldots, m'_n\rangle$$
where $U \in \mathrm{SU}(n)$ mixes the mode creation operators: $a^\dagger_i \mapsto \sum_j U_{ij} a^\dagger_j$.

The amplitudes $c_{m,m'}(U)$ are matrix elements of SU($n$) irreducible representations. Specifically, they are matrix elements of the **symmetric irrep** $\lambda = [N, 0, \ldots, 0]$, which corresponds to a single-row Young tableau. This irrep applies because bosonic states are symmetric under particle exchange. For symmetric irreps, `occupation_number(pattern)` converts a [Gelfand-Tsetlin pattern](../tutorial/states.md) to the corresponding Fock occupation list.

The translation between quantum optics and [GT pattern](../tutorial/states.md) language therefore has four steps:
1. Define the photon-number subspace as `λ = [N, 0, ..., 0]`.
2. Enumerate its basis of [GT patterns](../tutorial/states.md) with `basis = basis_states(λ)`.
3. Map each basis element to the Fock state whose occupation numbers are given by `occupation_number(pattern)`.
4. Compute the transition amplitude with `group_function(λ, final, initial, U)`. Here, $U$ is the $SU(n)$ unitary described above (for example, one provided by `su2_block(n, i, (α, β, γ))`), while `final` and `initial` are the final and initial states written as [GT patterns](../tutorial/states.md).




### Basis states and occupation numbers

The Hilbert space of $N$ photons in $n$ modes has dimension $\binom{N+n-1}{n-1}$. Its basis states can be enumerated as [GT patterns](../tutorial/states.md):

```@repl hom
using GroupFunctions

# Two photons in three modes
λ = [2, 0, 0];
basis = basis_states(λ);

# Each GT pattern corresponds to a Fock state
for b in basis
    occ = occupation_number(b)
    println("|", join(occ, ","), "⟩")
end
```


### Building SU($n$) matrices

Any SU($n$) matrix can be decomposed into SU(2) blocks that act on adjacent modes. The function `su2_block(n, i, (α, β, γ))` embeds an SU(2) rotation, parametrized by Euler angles, into modes $i$ and $i+1$:

```julia
# Beamsplitter on modes 1-2 (θ = π/2 for 50:50)
θ = float(π)/2
BS_12 = su2_block(3, 1, (0., θ, 0.))

# Beamsplitter on modes 2-3
BS_23 = su2_block(3, 2, (0., θ, 0.))

# Compose: first BS_12, then BS_23
U = BS_12 * BS_23
```


### Transition amplitudes

The function `group_function` computes the amplitude $\langle m' | U | m \rangle$:

```julia
λ = [2, 0, 0]
basis = basis_states(λ)

# Initial state: |2,0,0⟩
initial = basis[findfirst(b -> occupation_number(b) == [2,0,0], basis)]

# Final state: |1,1,0⟩
final = basis[findfirst(b -> occupation_number(b) == [1,1,0], basis)]

# Amplitude
amp = group_function(λ, final, initial, U)
prob = abs2(amp)
```

## Example: Hong-Ou-Mandel interference

Consider two photons entering a 50:50 beamsplitter through different input ports:
- Initial state: $|1,1\rangle$
- Beamsplitter: $a^\dagger_1 \mapsto \frac{1}{\sqrt{2}}(a^\dagger_1 + a^\dagger_2)$, $a^\dagger_2 \mapsto \frac{1}{\sqrt{2}}(a^\dagger_1 - a^\dagger_2)$

```@repl hom
λ = [2, 0];
basis = basis_states(λ);

initial = basis[findfirst(b -> occupation_number(b) == [1,1], basis)];

θ = float(π)/2;
BS = su2_block(2, 1, (0., θ, 0.));

for final in basis
    amp = group_function(λ, final, initial, BS)
    println("|", join(occupation_number(final), ","), "⟩: ", round(abs2(amp), digits=4))
end
```

In the expected HOM effect, destructive interference causes the $|1,1\rangle \to |1,1\rangle$ amplitude to vanish. The output state is $\frac{1}{\sqrt{2}}(|2,0\rangle + |0,2\rangle)$.
