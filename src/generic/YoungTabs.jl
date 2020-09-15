#import AbstractAlgebra: YoungTableau
using AbstractAlgebra
import Primes: prime
import LinearAlgebra:  dot, I
import Combinatorics:  permutations

const Content = Vector{T} where T <: Integer
const Irrep = Vector{T} where T <: Integer

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
    mat = Array(matrix_repr(pat))

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

@doc Markdown.doc"""
  generar_matriz(Y::Array{YoungTableau}, p::Perm) -> SparseMatrixCSC
> Return non-zero entries of the orthogonal irrep given by the permutation 'p'
> The information of the irrep is introduced via 'Y' which is a list of
> Standard Young tableaux

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> guilty = StandardYoungTableaux([3,2])
julia> generar_matriz(guilty, Perm([2,1,3,4,5]))
[1, 1]  =  -1.0
[2, 2]  =  1.0
[3, 3]  =  -1.0
[4, 4]  =  1.0
[5, 5]  =  1.0

julia> generar_matriz(guilty, Perm([1,3,2,4,5]))
[1, 1]  =  0.5
[2, 1]  =  0.866025
[1, 2]  =  0.866025
[2, 2]  =  -0.5
[3, 3]  =  0.5
[4, 3]  =  0.866025
[3, 4]  =  0.866025
[4, 4]  =  -0.5
[5, 5]  =  1.0
```
"""
function generar_matriz(patrones::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, p::Perm, irrep::Array{Int64,1})
    descom_en_trans = descomp_total(p)
    len = length(patrones)
    mat = I#spzeros(len,len)
    for (a,b) in descom_en_trans # a + 1 = b
        mat = generar_matriz(patrones, b, irrep)*mat
    end
    mat
end
function generar_matriz(lista_tablones::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, m::Int, irrep::Array{Int64,1})
    mat = spzeros(Float64,length(lista_tablones), length(lista_tablones))

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

    yup = *(size(mat|>Array)...)
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

function bloque_yamanouchi(mat::SparseMatrixCSC{Float64,Int64}, lista_tablones::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, i::Int64, k::Int64, m::Int64, irrep::Array{Int64,1})
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
##############################################################################
#
#   Codigo para representates double coset
#
##############################################################################
@doc Markdown.doc"""
    indice_tablon_semistandard(p::YoungTableau)
> Returns the index of the standard YoungTableau such that the function mapping
> the filling of the semistandard to the standard Tableau is non decreasing

# Examples:
```
julia> patrones_gt_prueba = genera_patrones([2,1,0]);
julia> young_prueba = map(YoungTableau, patrones_gt_prueba);
julia> for elemento in young_prueba
           icons_abstract_thee(elemento)
       end
```
"""
function indice_tablon_semistandard(tablon_semistandard::AbstractAlgebra.Generic.YoungTableau{Int64})
    particion = (tablon_semistandard.part) |> collect
    tablon_semistandard_v = YoungTableau(particion)# |> primero_lexi
    orden = sortperm(tablon_semistandard.fill)
    diccionario = generate_dictionary(orden)#Dict()
    nuevo_fill = [diccionario[x] for x in tablon_semistandard_v.fill]
    tablon_resultado_permutacion = AbstractAlgebra.fill!(tablon_semistandard_v, nuevo_fill)
    tablones_standard = StandardYoungTableaux(particion)
    etiquetas_semi = map(gen_etiqueta, map(x->x.fill, tablones_standard))
    
    return findfirst(x -> x ≈ gen_etiqueta(tablon_resultado_permutacion.fill), etiquetas_semi)
end

function generate_dictionary(lista::Array{Int64,1})
    fvars = Dict()
    for (n, f) in enumerate(lista)
        fvars[f] = n
    end
    fvars
end

function gen_etiqueta(lista::Array{Int64,1})
    len = length(lista)
    dot([Base.sqrt(prime(i)) for i in 1:len], lista)
end

@doc Markdown.doc"""
    content(p::Partition, λ::Irrep)
> Return the size of the vector which represents the partition.
> ADVERTENCIA content sin λ ignora los 0s de la irrep.

# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
julia> ss = YoungTableau(GTPattern([[2,1,0,0],[2,1,0],[2,1],[2]]));
julia> content(ss, [2,1,0,0])
[2,2,0,0]

julia> ss = YoungTableau(GTPattern([[2,1,0,0],[2,1,0],[2,1],[2]]));
julia> content(ss, [2,1,0,0])
[2,2,0]
```
"""
#function content(y::YoungTableau, λ::Irrep)
function content(y::AbstractAlgebra.Generic.YoungTableau{Int64}, λ::Irrep)
    relleno = y.fill
    if length(relleno) <= length(λ)
      len = length(λ)
    else
      len = length(relleno)
    end
    tablon_nuevo = tablon_standar_asociado_a_semiestandar(y).fill
    nuevo_relleno = deepcopy(relleno)
    anterior = relleno[1]
    for i in 1:len
      if i ∉ tablon_nuevo
        push!(nuevo_relleno, anterior)
      end
      if i > length(relleno)
        break
      end
      relleno[i] > anterior ? anterior = relleno[i] : continue
    end

    [count(y -> x == y,nuevo_relleno) for x in 1:len]
end

#function content(p::YoungTableau)
function content(p::AbstractAlgebra.Generic.YoungTableau{Int64})
    relleno = p.fill
    len = length(relleno)
    [count(y -> x == y,relleno) for x in 1:len]
end

@doc Markdown.doc"""
    tablon_standar_asociado_a_semiestandar(tablon_semistandard::YoungTableau)
> Returns a YoungTableau corresponding to the standard YoungTableau such that
> f is non decreasing.
"""
function tablon_standar_asociado_a_semiestandar(tablon_semistandard::AbstractAlgebra.Generic.YoungTableau{Int64})
    particion = (tablon_semistandard.part) |> collect
    tablon_standard = YoungTableau(particion)# |> primero_lexi

    orden = sortperm(tablon_semistandard.fill)
    diccionario = generate_dictionary(orden)#Dict()
    nuevo_fill = [diccionario[x] for x in tablon_standard.fill]

    tablon_resultado_permutacion = AbstractAlgebra.fill!(tablon_standard, nuevo_fill)
    tablones_standard = StandardYoungTableaux(particion)
    etiquetas_semi = map(gen_etiqueta, map(x->x.fill, tablones_standard))
    
    pos = findfirst(x -> x ≈ gen_etiqueta(tablon_resultado_permutacion.fill), etiquetas_semi)
    tablones_standard[pos]
end

@doc Markdown.doc"""
    Θ(patron_semi::YoungTableau)
> Computes coefficient Θ. Returns a Float64
"""
function Θ(patron_semi::AbstractAlgebra.Generic.YoungTableau{Int64}, irrep::Array{Int64,1})
    relleno_semi = patron_semi.fill
    tablones_standard = StandardYoungTableaux((patron_semi.part) |> collect)
    tablon_standard = tablon_standar_asociado_a_semiestandar(patron_semi)
    relleno_standard = tablon_standard.fill
    n = ((patron_semi.fill) |> collect |> length)
    n = irrep |> length
    
    parejas = zip(relleno_standard, relleno_semi) |> collect
    prod = 1.0
    for k in 1:n
        α = map(first,filter((xx -> (last(xx) == k)), parejas))
        if length(α) > 1
            for indx in 1:length(α), indy in indx+1:length(α)
                x,y = α[indx], α[indy]
                if x > y
                  prod *= (1 + (1/axialdistance(tablon_standard, y, x)))
                else
                  prod *= (1 + (1/axialdistance(tablon_standard, x, y)))
                end
            end
        end
    end
    prod
end

@doc Markdown.doc"""
calcula_proto_permutacion(proto::Array{Int64,1})
Notese que el orden en que se ingresa la matriz es igual a la traspuesta.
Esto es debido a la forman en la que Julia recorre matrices.
> calcula_proto_permutacion([2 0 1 1])
1
1
2
2
"""
function calcula_proto_permutacion(proto::Array{Int64,1})
    len = length(proto) |> sqrt |> Int
    new_proto = reshape(proto, len, len)'

    mult = 0
    yard = Array{Int64,1}[]
    for i in 1:len
        densidad = vcat(fill.(1:len, new_proto[mult*len + 1:len*(mult + 1)])...)
        push!(yard, densidad)
        mult += 1
    end
    vcat(yard...)
end

@doc Markdown.doc"""
    genera_funcion(patron_semi::YoungTableau, irrep::Vector{T})
    genera_funcion(patron_semi::YoungTableau)
> Genera un diccionario con la función entre un tablón standar
> y uno semistandar.
# Examples:
```
julia> t_u = YoungTableau([2,2, 1]);
julia> fill!(t_u, [1,2,3,3,4]);
julia> genera_funcion(t_u)
Dict(4=>4, 2=>2, 3=>3, 5=>5, 1=>1)
```
"""
function genera_funcion(patron_semi::AbstractAlgebra.Generic.YoungTableau{Int64}, irrep::Array{Int64,1}) 
    len = length(irrep)
    relleno_semi = patron_semi.fill

    if length(relleno_semi) >= len
      return genera_funcion(patron_semi)
    end

    tablones_standard = StandardYoungTableaux((patron_semi.part) |> collect)
    tablon_standard = tablon_standar_asociado_a_semiestandar(patron_semi)
    relleno_standard = tablon_standard.fill
    
    parejas = zip(relleno_standard, relleno_semi) |> collect
    dd = Dict{Int64, Int64}()

    for (va, viene) in parejas
      dd[va] = viene
    end
    nuevos = setdiff!(collect(1:len), relleno_standard)
    for i in nuevos
      cercano = sort(relleno_standard, by=(x -> abs(x-i))) |> first
      dd[i] = dd[cercano]
    end

    dd
end

function genera_funcion(patron_semi::AbstractAlgebra.Generic.YoungTableau{Int64})
    relleno_semi = patron_semi.fill
    tablones_standard = StandardYoungTableaux((patron_semi.part) |> collect)
    tablon_standard = tablon_standar_asociado_a_semiestandar(patron_semi)
    relleno_standard = tablon_standard.fill
    
    parejas = zip(relleno_standard, relleno_semi) |> collect
    dd = Dict{Int64, Int64}()

    for (va, viene) in parejas
      #dd[viene] = va
      dd[va] = viene
    end
    dd
end

@doc Markdown.doc"""
    calcular_sα(c::Content)
> Return the size of the vector which represents the partition.

# Examples:
```
julia> c = [0,1,2]; calcular_sα(c)
```
"""
function calcular_sα(c::Content)
    inferior = 1
    superior = length(c) 
    list_perm_output = Perm{Int64}[]
    for sub_cjto in c
        if sub_cjto == 0
            continue
        elseif sub_cjto == 1
            inferior += 1
        elseif sub_cjto > 1 && length(list_perm_output) == 0
            conjunto = collect(inferior:(inferior + sub_cjto - 1))
            permutaciones = permutations(conjunto) |> collect

            for p in permutaciones
                lista = sortperm(vcat(collect(1:inferior - 1), p, collect((inferior + sub_cjto):superior)))
                push!(list_perm_output, Perm(lista))
            end
            inferior += sub_cjto
        else
            conjunto = collect(inferior:(inferior + sub_cjto - 1))
            permutaciones = permutations(conjunto)

            tmp = Perm{Int64}[]
            for p in permutaciones
                lista = sortperm(vcat(collect(1:inferior - 1), p , collect((inferior + sub_cjto):superior)))
                for elem in list_perm_output
                  push!(tmp, elem*Perm(lista))
                end
            end
            list_perm_output = tmp
            inferior += sub_cjto
        end
        push!(list_perm_output, Perm(1:superior |> collect))
        unique!(list_perm_output)
    end
    if length(list_perm_output) == 0
      return [Perm(1:superior |> collect)]
    end
    list_perm_output
end

function calcular_sα(tablon::AbstractAlgebra.Generic.YoungTableau{Int64})
  contenido = content(tablon)
  inferior = 1
  superior = length(contenido)
  list_perm_output = Perm[]
  for sub_cjto in contenido
      if sub_cjto == 0
          continue
      elseif sub_cjto == 1
          inferior += 1
      elseif sub_cjto > 1
          conjunto = collect(inferior:(inferior + sub_cjto - 1))
          permutaciones = permutations(conjunto)

          for p in permutaciones
              lista = vcat(collect(1:inferior - 1), p, collect((inferior + sub_cjto):superior))
              push!(list_perm_output, Perm(lista))
          end
          inferior += 1
      end
      unique!(list_perm_output)
  end
  if length(list_perm_output) == 0
    return [Perm(1:superior |> collect)]
  end
  list_perm_output
end

function YoungTableau(tab::GTPattern)
    filas = tab.filas
    len = filas[1] |> length
    conjunto_contenido = [obtener_diferencias_patron(tab, x) for x in 1:len]
    p = Partition(filter(x -> x>0, filas[1]))
    Generic.YoungTableau(p, vcat(conjunto_contenido...))
end
