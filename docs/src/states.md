# States
The preferred method to denote the basis states is by using Gelfand–Tsetlin patterns.  
Without going into detail, this can be seen in the paper by Alex et al.  
The top row of the pattern is the partition that labels the irrep.  
The rest of the entries should satisfy the betweenness condition, which I now explain.  

I start with an example. (Due to $\LaTeX$ limitations, I denote the pattern as left-justified.)  
The pattern is:  
```math
\begin{pmatrix}
c_{0,0} & c_{0,1} & c_{0,2} \\
c_{1,0} & c_{1,1} & \\
c_{2,0} & &  
\end{pmatrix}
```


The betweenness condition imposes the requirement: 
```math
c_{i-1,j} \geq c_{i,j} \geq c_{i-1,j+1}.
```


To introduce a `GTPattern`, simply write an array of arrays with the content of
the pattern as argument:
```julia
GTPattern([[2,1,0], [2,1], [2]])
```
