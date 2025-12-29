# Basic quantum optics using GroupFunctions.jl

## Introduction

Linear optical networks (beamsplitters, phase shifters, etc. )act on $n$ optical modes via SU($n$) transformations. A key object is the transition amplitude between Fock states under such transformations.

Consider $N$ photons distributed across $n$ modes. The Fock state $|m_1, m_2, \ldots, m_n\rangle$ with $\sum_i m_i = N$ transforms under passive linear optics as:
$$|m_1, \ldots, m_n\rangle \xrightarrow{U} \sum_{m'} c_{m,m'}(U) |m'_1, \ldots, m'_n\rangle$$
where $U \in \mathrm{SU}(n)$ mixes the mode creation operators: $a^\dagger_i \mapsto \sum_j U_{ij} a^\dagger_j$.

The amplitudes $c_{m,m'}(U)$ are matrix elements of SU($n$) irreducible representations, specifically, the **symmetric irrep** $\lambda = [N, 0, \ldots, 0]$, corresponding to a single-row Young tableau. This is because bosonic states are symmetric under particle exchange. Importantly, for symmetric irreps, the p-weight of a Gelfand-Tsetlin pattern (used internally by this library) directly gives the occupation numbers, more natural in quantum optics.  

The basic items to translate between quantum optics and GT pattern language are thus the following:
1. Define photon number subspace  `λ = [N, 0, ..., 0]` 
2. Enumerate the basis of the above subspace using GT patterns  `basis = basis_states(λ)` 
3. Each of the basis elements corresponds to the Fock state with occupation number provided by  `pweight(pattern)` 
4. Transition amplitude is read out by `group_function(λ, final, initial, U)`, where $U$ is the $SU(n)$ unitary in previous paragraphs (e.g. provided by `su2_block(n, i, (α, β, γ))`), `final` and `initial` are the final and initial states (written as GT patterns).




### Basis states and occupation numbers

The Hilbert space of $N$ photons in $n$ modes has dimension $\binom{N+n-1}{n-1}$. We enumerate basis states as GT patterns:

```julia
using GroupFunctions

# Two photons in three modes
λ = [2, 0, 0]
basis = basis_states(λ)

# Each GT pattern corresponds to a Fock state via pweight
for b in basis
    occ = pweight(b)
    println("|", join(occ, ","), "⟩")
end
```

Expected output:
```
|2,0,0⟩
|1,1,0⟩
|1,0,1⟩
|0,2,0⟩
|0,1,1⟩
|0,0,2⟩
```

### Building SU($n$) matrices

Any SU($n$) matrix can be decomposed into SU(2) blocks acting on adjacent modes. The function `su2_block(n, i, (α, β, γ))` embeds an SU(2) rotation (parametrized by Euler angles) into modes $i$ and $i+1$:

```julia
# Beamsplitter on modes 1-2 (θ = π/2 for 50:50)
θ = float(π)/2
BS_12 = su2_block(3, 1, (0., θ, 0.))

# Beamsplitter on modes 2-3
BS_23 = su2_block(3, 2, (0., θ, 0.))

# Compose: first BS_12, then BS_23
U = BS_12 * BS_23
```

TODO: Verify sign/phase conventions for `su2_block`. The beamsplitter mixing $a^\dagger_1 \mapsto \frac{1}{\sqrt{2}}(a^\dagger_1 + a^\dagger_2)$ corresponds to which Euler angles?

### Transition amplitudes

The amplitude $\langle m' | U | m \rangle$ is computed via `group_function`:

```julia
λ = [2, 0, 0]
basis = basis_states(λ)

# Initial state: |2,0,0⟩
initial = basis[findfirst(b -> pweight(b) == [2,0,0], basis)]

# Final state: |1,1,0⟩
final = basis[findfirst(b -> pweight(b) == [1,1,0], basis)]

# Amplitude
amp = group_function(λ, final, initial, U)
prob = abs2(amp)
```

## Example: Hong-Ou-Mandel interference

Two photons entering a 50:50 beamsplitter from different input ports:
- Initial state: $|1,1\rangle$
- Beamsplitter: $a^\dagger_1 \mapsto \frac{1}{\sqrt{2}}(a^\dagger_1 + a^\dagger_2)$, $a^\dagger_2 \mapsto \frac{1}{\sqrt{2}}(a^\dagger_1 - a^\dagger_2)$

```julia
λ = [2, 0]
basis = basis_states(λ)

initial = basis[findfirst(b -> pweight(b) == [1,1], basis)]

θ = float(π)/2
BS = su2_block(2, 1, (0., θ, 0.))

for final in basis
    amp = group_function(λ, final, initial, BS)
    println("|", join(pweight(final), ","), "⟩: ", round(abs2(amp), digits=4))
end
```

Expected (HOM effect): The $|1,1\rangle \to |1,1\rangle$ amplitude vanishes due to destructive interference. Output is $\frac{1}{\sqrt{2}}(|2,0\rangle + |0,2\rangle)$.

TODO: Verify output matches expected HOM signature. Check if beamsplitter convention gives the correct phases.


