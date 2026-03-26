# Example: SU(2)-invariant qubit transmission

This example demonstrates encoding a qubit in four optical modes such that arbitrary polarization mixing during fiber transmission does not affect the logical state. The protocol uses two photons distributed across four modes, where modes 1–2 and 3–4 each experience the same unknown unitary $U$.

The protocol can be summarized as:

1. State initialization of 2 photons across four modes: $\alpha|2,0,0,0\rangle + \beta|0,0,0,2\rangle$,
2. Preprocessing: balanced beamsplitter on modes (2,3), then swap modes (1,2),
3. Transmission through a fiber: first modes (1,2), modelling horizontal and vertical polarization, then (3,4). Both pairs undergo identical transformation $U$, so the total is $U\oplus U$.
4. Potprocessing through swap of modes (1,2), then beamsplitter on modes (1,4).
5. Postselection on detecting 0 photons in mode 1, 1 photon in mode 4.

With this, the resulting state is $\alpha|0,1,0,1\rangle + \beta|0,0,1,1\rangle$ — the original qubit in modes 2–3, independent of $U$. The probability is affected by the exact form of $U$, but the normalization is irrelevant under postselection.

## Implementation

First, set up the Hilbert space and identify the relevant Fock states:

```julia
using GroupFunctions
using SymEngine
using LinearAlgebra: I

# Two photons in four modes
λ = [2, 0, 0, 0]
basis = basis_states(λ)

# Helper: find pattern by occupation
state(occ) = basis[findfirst(b -> pweight(b) == occ, basis)]

# Initial and postselected states
initial_1 = state([2,0,0,0])
initial_2 = state([0,0,0,2])
final_1   = state([0,1,0,1])
final_2   = state([0,0,1,1])
```

Next, construct the optical circuit symbolically. The fiber unitary $U$ is left as a symbolic $2\times 2$ matrix, embedded block-diagonally to act identically on both mode pairs:

```julia
# Build symbolic unitaries
M_23 = bs_block_symbolic(4, 2)           # BS on modes 2-3
S_12 = swap_block_symbolic(4, 1)         # swap modes 1-2
M_14 = bs_block_symbolic(4, (1, 4))      # BS on modes 1-4 (non-adjacent)

# Unknown fiber unitary: symbolic 2×2 block, embedded block-diagonally
U_2x2 = su2_block_symbolic(2, 1, prefix="u")
U_fiber = Matrix{Basic}(I, 4, 4)
U_fiber[1:2, 1:2] = U_2x2[1:2, 1:2]
U_fiber[3:4, 3:4] = U_2x2[1:2, 1:2]

# Full transformation
U_total = M_14 * S_12 * U_fiber * S_12 * M_23
```

Now compute the transition amplitudes. Since the initial state is a superposition $\alpha|2,0,0,0\rangle + \beta|0,0,0,2\rangle$, we compute the amplitudes for each computational basis state separately:

```julia
# Compute symbolic transition amplitudes
amp_11 = group_function_sym(λ, final_1, initial_1, U_total)  # |2,0,0,0⟩ → |0,1,0,1⟩
amp_12 = group_function_sym(λ, final_2, initial_1, U_total)  # |2,0,0,0⟩ → |0,0,1,1⟩
amp_21 = group_function_sym(λ, final_1, initial_2, U_total)  # |0,0,0,2⟩ → |0,1,0,1⟩
amp_22 = group_function_sym(λ, final_2, initial_2, U_total)  # |0,0,0,2⟩ → |0,0,1,1⟩

println("α|2,0,0,0⟩ → |0,1,0,1⟩: ", amp_11)
println("α|2,0,0,0⟩ → |0,0,1,1⟩: ", amp_12)
println("β|0,0,0,2⟩ → |0,1,0,1⟩: ", amp_21)
println("β|0,0,0,2⟩ → |0,0,1,1⟩: ", amp_22)
```

## Interpreting the result

By linearity, the full output state (conditioned on postselection) is obtained by superposing the contributions:

$$|\psi_{\mathrm{out}}\rangle \propto \alpha \bigl(\mathrm{amp}_{11}|0,1,0,1\rangle + \mathrm{amp}_{12}|0,0,1,1\rangle\bigr) + \beta \bigl(\mathrm{amp}_{21}|0,1,0,1\rangle + \mathrm{amp}_{22}|0,0,1,1\rangle\bigr).$$

The symbolic calculation reveals that $\mathrm{amp}_{12} = \mathrm{amp}_{21} = 0$ and $\mathrm{amp}_{11} = \mathrm{amp}_{22} \equiv A$, where $A$ is independent of the fiber unitary $U$. The postselected state therefore simplifies to:

$$|\psi_{\mathrm{out}}\rangle = \frac{1}{\mathcal{N}}\bigl(\alpha|0,1,0,1\rangle + \beta|0,0,1,1\rangle\bigr),$$

where $\mathcal{N} = |A|\sqrt{|\alpha|^2 + |\beta|^2}$ is a normalization factor arising from the postselection probability. The qubit encoded in modes 2–3 is recovered with its original coefficients $\alpha$ and $\beta$, regardless of the unknown transformation $U$ applied during transmission.
