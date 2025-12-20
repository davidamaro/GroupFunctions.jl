# Sum rules

First, prepare the matrix using the simple factorization.
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

mat4 = block12_a * block23_a * block12_b
mat4c1 = block12_a * block23_a * block12_b * block34_a * block23_b
```
Then obtain the basis states:
```julia
basis = basis_states([2,0,0,0])
```
Finally, compute the rate:

```julia
rate1 = abs(group_function([2,0,0,0], basis[9], basis[9], mat4))^2 +
        abs(group_function([2,0,0,0], basis[9], basis[4], mat4))^2 +
        abs(group_function([2,0,0,0], basis[9], basis[7], mat4))^2
rate2 = abs(group_function([2,0,0,0], basis[9], basis[9], mat4c1))^2 +
        abs(group_function([2,0,0,0], basis[9], basis[4], mat4c1))^2 +
        abs(group_function([2,0,0,0], basis[9], basis[7], mat4c1))^2
rate1 ≈ rate2
```
