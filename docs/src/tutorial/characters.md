# Characters
SU(d) characters can be evaluated by computing the trace of a representation.

## Example: SU(3) character check
Here we present an example with a known formula for the character in the irrep
$\lambda = (2, 1, 0)$, and compare it against the trace of the representation
computed using `group_function`.

```julia
using GroupFunctions
using RandomMatrices
using LinearAlgebra: det, eigvals, tr

irrep = [2, 1, 0]
U = rand(Haar(2), 3)  # 3×3 Haar-random unitary matrix
U_su = U / det(U)^(1 / 3)

rep, _states = group_function(irrep, U_su)
lambdas = eigvals(U_su)
ideal_character = lambdas[1] / lambdas[2] + lambdas[2] / lambdas[1] +
                  lambdas[1] / lambdas[3] + lambdas[3] / lambdas[1] +
                  lambdas[2] / lambdas[3] + lambdas[3] / lambdas[2] + 2
tr(rep) ≈ ideal_character
```
