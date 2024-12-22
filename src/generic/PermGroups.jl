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
    completa = Tuple{Int64,Int64}[]
    for elem in lista_ciclos
        descomponer_ciclo!(elem, completa)
    end
    completa
end

function descomponer_ciclo!(ciclo::SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true}, output::Array{Tuple{Int64,Int64},1})
    lista = Tuple{Int64,Int64}[]
    for i in 1:length(ciclo) - 1
        push!(lista, (ciclo[i], ciclo[i+1]))
    end

    for elem in lista
        individual!(elem, output)
    end
end

# function individual!(tranposicion::Tuple{Int64,Int64},lista::Array{Tuple{Int64,Int64},1})
    # j,k = tranposicion
    # if k < j
      # j,k = k,j
    # end
    # for i in k:-1:j+1
        # push!(lista, (i-1,i))
    # end
    # for i in j+1:k-1
        # push!(lista, (i, i+1))
    # end
# end
function individual!(tranposicion::Tuple{Int64,Int64}, lista::Array{Tuple{Int64,Int64},1})
    j, k = tranposicion
    if k < j
        j, k = k, j
    end
    
    # Pre-allocate space for all upcoming pushes
    needed_capacity = 2(k - j)
    sizehint!(lista, length(lista) + needed_capacity)
    
    # Combine both loops into a single pass with @inbounds for performance
    @inbounds begin
        # First pass (k down to j+1)
        for i in k:-1:j+1
            push!(lista, (i-1, i))
        end
        
        # Second pass (j+1 up to k-1)
        for i in j+1:k-1
            push!(lista, (i, i+1))
        end
    end
end
