# Immanants

First discovered by [Kostant](http://www.jstor.org/stable/2152885)
and then extended by [de Guise et al.](https://arxiv.org/pdf/1511.01851),
[immanants](https://en.wikipedia.org/wiki/Immanant) of unitary matrices can be computed by certain sums of group functions.

Briefly, in this section of the documentation, I give some examples of this
relation.
First, compute the unitary matrix:
```julia
α1,β1,γ1 = rand(Float64,3)
block12_a = su2_block(4,1,(α1,β1,γ1))
α2,β2 = rand(Float64,2)
block23_a = su2_block(4,2,(α2,β2,α2))
α3,β3,γ3 = rand(Float64,3)
block12_b = su2_block(4,1,(α3,β3,γ3))
α4,β4 = rand(Float64,3)
block34_a = su2_block(4,3,(α4,β4,α4))
α5,β5 = rand(Float64,2)
block23_b = su2_block(4,2,(α5,β5,α5))
α6,β6,γ6 = rand(Float64,3)
block12_c = su2_block(4,1,(α6,β6,γ6))

mat4 = block12_a * block23_a * block12_b * block34_a * block23_b * block12_c
```
Then obtain the states with a given p-weight (see the Alex et al. paper).
```julia
basis = basis_states([3,0,0,0])

state_x = filter(x -> pweight(x) == [0,1,1,1], basis)[1]
state_y = filter(x -> pweight(x) == [1,0,2,0], basis)[1]
```
Finally, show the equality between the permanent and a group function.

```julia
group_function([3,0,0,0], state_x, state_y, mat4) ≈ permanent(mat4[[1,2,3], [2,2,4]])/sqrt(2)
```
