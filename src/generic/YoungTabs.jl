#import AbstractAlgebra: YoungTableau
using AbstractAlgebra

@doc Markdown.doc"""
    axialdistance(Y::YoungTableau, i, j)
> Return the hook-length of an element in `Y` at position `(i,j)`, i.e the
> number of cells in the `i`-th row to the right of `(i,j)`-th box, plus the
> number of cells in the `j`-th column below the `(i,j)`-th box, plus `1`.
>
> Return `0` for `(i,j)` not in the tableau `Y`.

# Examples
```
julia> y = YoungTableau([4,3,1])
┌───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │
├───┼───┼───┼───┘
│ 5 │ 6 │ 7 │
├───┼───┴───┘
│ 8 │
└───┘

julia> axialdistance(y, 1,1)
6

julia> axialdistance(y, 1,3)
3

julia> axialdistance(y, 2,4)
0
```
"""
function axialdistance(Y::AbstractAlgebra.Generic.YoungTableau{Int64}, u::Int64, v::Int64)
  #fila, columna
  i,j = encontrar_posicion(Y, u)
  k,l = encontrar_posicion(Y, v)

  l - k - (j - i)
end

@doc Markdown.doc"""
    matrix_repr(Y::YoungTableau)
> Construct sparse integer matrix representing the tableau.

# Examples:
```
julia> y = YoungTableau([4,3,1]);


julia> matrix_repr(y)
3×4 SparseArrays.SparseMatrixCSC{Int64,Int64} with 8 stored entries:
  [1, 1]  =  1
  [2, 1]  =  5
  [3, 1]  =  8
  [1, 2]  =  2
  [2, 2]  =  6
  [1, 3]  =  3
  [2, 3]  =  7
  [1, 4]  =  4
```
"""
function encontrar_posicion(Y::AbstractAlgebra.Generic.YoungTableau{Int64}, entrada::Int64)
   k=1
   for (idx, p) in enumerate(Y.part)
      for (idy, q) in enumerate(Y.fill[k:k+p-1])
        if q == entrada
          return idx, idy
        end
      end
      k += p
   end
   0,0
end

function determinar_coeficiente_irrep_yamanouchi(Y::AbstractAlgebra.Generic.YoungTableau{Int64}, u::Integer)
  #fila, columna
  v = u - 1
  i,j = encontrar_posicion(Y, u)
  k,l = encontrar_posicion(Y, v)
  
  if i == k
    return 1
  elseif j == l
    return -1 
  else
    return 0
  end
end

@doc Markdown.doc"""
    primero_lexi(YoungTableau) -> YoungTableau
    Computes the first ---in lexicographic order---
    Standard Tableaux.
# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
julia> pat = YoungTableau([2,1])
julia> primero_lexi(pat)
┌───┬───┬
│ 1 │ 3 │
├───┼───┼
│ 2 │
├───┼
```
"""
function primero_lexi(pat::AbstractAlgebra.Generic.YoungTableau{Int64})
    mat = Matrix(matrix_repr(pat))

    orco = zeros(Int, size(mat))
    inx = 1
    for (ind,val) in enumerate(mat)
        if val != 0
            orco[ind] = inx
            inx += 1
        end
    end
    fill!(pat, filter(x -> x > 0,reshape(orco', *(size(orco)...)) ) )
    pat
end
@doc Markdown.doc"""
    StandardYoungTableaux(part::Array{Int64,1}) -> list of YoungTableaux
> Return a list of Standard YoungTableaux.
# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
julia> StandardYoungTableaux([2,1])
[1 3; 2 0]
[1 2; 3 0]
"""
function StandardYoungTableaux(part::Array{Int64,1}) 
    part = filter(x -> x > 0, part)
    len = length(part)
    flat = zeros(Int, sum(part))
    flat[1:len] = part
    primero = primero_lexi(YoungTableau(part))
    lista = [deepcopy(primero)]
    for i in 2:dim(YoungTableau(part))
        push!(lista, deepcopy(encontrar_malo_imp!(primero)))
    end
    lista
end

function generar_matriz(lista_tablones::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, m::Int, irrep::Array{Int64,1})
    mat = spzeros(Int64,length(lista_tablones), length(lista_tablones))

    for i in 1:length(lista_tablones), k in 1:length(lista_tablones)
        if i == k && mat[i,i] == 0
            mat[i,i] = determinar_coeficiente_irrep_yamanouchi(lista_tablones[i], m)
        elseif i < k
            bloque_yamanouchi(mat, lista_tablones, i, k, m, irrep)
        end
    end


    mat
end

function encontrar_malo_imp!(tab::AbstractAlgebra.Generic.YoungTableau{Int64})
    i_max = 0
    jm = 0
    malo = 0

    lista_recorridos = typeof((1,2))[]

    for l in 1:sum(tab.part)
        i,j = encontrar_posicion(tab, l)
        push!(lista_recorridos,(i,j))
        if i >= i_max
            i_max = i
            jm = l
        else
            malo = l
            break
        end

    end

    mat = matrix_repr(tab)
    pos_maximo = encontrar_posicion(tab, jm)
    pos_malo = encontrar_posicion(tab, malo)

    pos_maximo = encontrar_esquina(lista_recorridos)

    mat[pos_maximo...] = malo

    filtrado = filter(x -> x != pos_maximo, sort(lista_recorridos,by = x->(x[2], x[1])) )
    for (ind, val) in enumerate(filtrado)
        if val == pos_maximo
            continue
        else
            mat[val...] = ind
        end
    end

    yup = *(size(mat|>Matrix)...)
    fill!(tab, filter(x -> x > 0, reshape(mat', yup)))

    tab
end
@doc Markdown.doc"""
    encontrar_esquina(lista tuplas)
"""
function encontrar_esquina(lista::Array{Tuple{Int64,Int64},1})
    (fila_malo, col_malo) = lista[end]
    i_max, jm = (1,1)
    if col_malo == 1
        return lista[end-1]
    else
       pre_col = col_malo - 1
        while (i_max <= fila_malo && pre_col > 0)
            (i_max, jm) = filter(x -> x[2] == pre_col, lista[1:end-1])[end]
            pre_col -= 1
        end
    end
    (i_max, jm)
end

function bloque_yamanouchi(mat::SparseMatrixCSC{Int64,Int64}, lista_tablones::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, i::Int64, k::Int64, m::Int64, irrep::Array{Int64,1})
    tab_1 = lista_tablones[i]
    tab_2 = lista_tablones[k]
    if i < k && lista_tablones[i] == intercambio(lista_tablones[k], m,irrep)
        mat[i,i] = -((axialdistance(tab_1,m, m-1))^(-1))
        mat[i,k] = sqrt(1-((axialdistance(tab_1,m, m-1))^(-2)))
        mat[k,i] = sqrt(1-((axialdistance(tab_1,m, m-1))^(-2)))
        mat[k,k] = ((axialdistance(tab_1,m, m-1))^(-1))
    end
    mat
end

function intercambio(Y::AbstractAlgebra.Generic.YoungTableau{Int64}, m::Int64, irrep::Array{Int64,1})
    V = deepcopy(Y)
    relleno = V.fill
    @assert m > 1 && m <= maximum([length(irrep), length(relleno)])
    posiciones = Int[]
    for (indx, elem) in enumerate(relleno)
        if elem == m || elem == m - 1
            push!(posiciones, indx)
        end
        if length(posiciones) == 2
            break
        end
    end

    if length(posiciones) > 1
      tmp = relleno[posiciones[1]]
      relleno[posiciones[1]] = relleno[posiciones[2]]
      relleno[posiciones[2]] = tmp
      fill!(V, relleno)
    else
      relleno[posiciones[1]] = m
      fill!(V, relleno)
    end

    V
end
