# Sum rules

First, prepare the matrix using the simple factorization.
```julia
    α1,β1,γ1 = rand(Float64,3)
    xx=su2_block(4,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=su2_block(4,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=su2_block(4,1,(α3,β3,γ3))
    α4,β4 = rand(Float64,3)
    xx2=su2_block(4,3,(α4,β4,α4))
    α5,β5 = rand(Float64,2)
    yy2=su2_block(4,2,(α5,β5,α5))
    α6,β6,γ6 = rand(Float64,3)
    zz2=su2_block(4,1,(α6,β6,γ6))

    mat4 = xx*yy*zz
    mat4c1 = xx*yy*zz*xx2*yy2
```
Then obtain the basis states:
```julia
    welcome = basis_states([2,0,0,0])
```
Finally, compute the rate

```julia
    rate1 = abs( group_function([2,0,0,0], welcome[9], welcome[9], mat4) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[4], mat4) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[7], mat4) )^2
    rate2 = abs( group_function([2,0,0,0], welcome[9], welcome[9], mat4c1) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[4], mat4c1) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[7], mat4c1) )^2
    rate1 ≈ rate2
```
