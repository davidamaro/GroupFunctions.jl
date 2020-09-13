@doc Markdown.doc"""
    descomp_total(g::Perm)
> Return the adjacent transposition decomposition of 'g'

```jldoctest; setup = :(using AbstractAlgebra)
julia> G = SymmetricGroup(5); g = Perm([3,4,5,2,1])
(1,3,5)(2,4)

julia> descomp_total(g)
true
```
"""
function descomp_total(original::Perm)
    lista_ciclos = cycles(original)
    #@show lista_ciclos
    completa = []
    for elem in lista_ciclos
        #push!(completa, descomponer_ciclo(elem))
        descomponer_ciclo!(elem, completa)
    end
    completa
end

function descomponer_ciclo!(ciclo, output)
    lista = Tuple{Int64,Int64}[]
    for i in 1:length(ciclo) - 1
        push!(lista, (ciclo[i], ciclo[i+1]))
    end
    #lista
    #output = Tuple{Int64,Int64}[]
    for elem in lista
        individual!(elem, output)
    end
    #output
end

function individual!(transp,lista)
    j,k = transp
    if k < j
      j,k = k,j
    end
    for i in k:-1:j+1
        #@show (i-1, i)
        push!(lista, (i-1,i))
        #prepend!(lista, (i-1,i))
    end
    #println("luego")
    for i in j+1:k-1
        #@show (i,i+1)
        push!(lista, (i, i+1))
        #prepend!(lista, (i, i+1))
    end
end
