# Immanants

First discovered by [Kostant](http://www.jstor.org/stable/2152885) 
and then extended by [de Guise et al](https://arxiv.org/pdf/1511.01851), 
[immanants](https://en.wikipedia.org/wiki/Immanant) of unitary matrices can be computed by certain sums of group functions.

Briefly, in this section of the documentation, I give some examples of this
relation.
First, I compute the unitary matrix:
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

    mat4 = xx*yy*zz*xx2*yy2*zz2
```
Then I obtain the states with certain p-weight (see Alex et al paper).
```julia
    welcome = basis_states([3,0,0,0])

    edox = filter(x -> pweight(x) == [0,1,1,1] , welcome)[1]
    edoy = filter(x -> pweight(x) == [1,0,2,0], welcome)[1]
```
Finally, I show the equality between the permanent and a group function.

```julia
    group_function([3,0,0,0], edox, edoy, mat4) ≈ permanent(mat4[[1,2,3], [2,2,4]])/sqrt(2)
```
