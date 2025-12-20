using AbstractAlgebra
import LinearAlgebra:  dot, I
import Combinatorics:  permutations
import SymEngine: Basic

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
  i,j = determine_position(Y, u)
  k,l = determine_position(Y, v)

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
function determine_position(tableau::AbstractAlgebra.Generic.YoungTableau{Int64}, entry::Int64)
   row_start = 1
   @inbounds for (row_idx, row_len) in enumerate(tableau.part)
      for col_idx in 1:row_len
        if tableau.fill[row_start + col_idx - 1] == entry
          return row_idx, col_idx
        end
      end
      row_start += row_len
   end
   return 0, 0
end

function determinar_coeficiente_irrep_yamanouchi(Y::AbstractAlgebra.Generic.YoungTableau{Int64}, u::Integer)
  v = u - 1
  i,j = determine_position(Y, u)
  k,l = determine_position(Y, v)
  
  if i == k
    return Basic(1)
  elseif j == l
    return Basic(-1)
  else
    return Basic(0)
  end
end

@doc Markdown.doc"""
    first_young_tableau_lexicographic(YoungTableau) -> YoungTableau
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
function first_young_tableau_lexicographic(pat::AbstractAlgebra.Generic.YoungTableau{Int64})
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
    first_tab = first_young_tableau_lexicographic(YoungTableau(part))
    tableaux = [deepcopy(first_tab)]
    for _ in 2:dim(YoungTableau(part))
        push!(tableaux, deepcopy(reposition_out_of_order!(first_tab)))
    end
    tableaux
end

@doc Markdown.doc"""
    generate_matrix(Y::Array{YoungTableau}, p::Perm, λ::Array{Int64,1}) -> SparseMatrixCSC
> Return non-zero entries of the orthogonal irrep given by the permutation 'p'
> The information of the irrep is introduced via 'Y' which is a list of
> Standard Young tableaux

# Examples
```
julia> guilty = StandardYoungTableaux([3,2])
julia> generate_matrix(guilty, Perm([2,1,3,4,5]), [3,2])
[1, 1]  =  -1.0
[2, 2]  =  1.0
[3, 3]  =  -1.0
[4, 4]  =  1.0
[5, 5]  =  1.0

julia> generate_matrix(guilty, Perm([1,3,2,4,5]), [3,2])
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
const GENERATE_MATRIX_CACHE = Dict{Tuple{Perm{Int64},Tuple}, SparseMatrixCSC{Basic,Int64}}()

function generate_matrix(patrones::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, p::Perm{Int64}, irrep::Array{Int64,1})
    cache_key = (p, Tuple(irrep))
    if haskey(GENERATE_MATRIX_CACHE, cache_key)
        return GENERATE_MATRIX_CACHE[cache_key]
    end

    transposition_factors = descomp_total(p)
    len::Int64 = length(patrones)
    mat::SparseMatrixCSC{Basic,Int64} = spzeros(Basic, len, len)
    @simd for i in 1:len
      @inbounds mat[i,i] = Basic(1)
    end
    for (_, transposition) in transposition_factors # a + 1 = b
        mat = generate_matrix(patrones, transposition, irrep) * mat
    end

    GENERATE_MATRIX_CACHE[cache_key] = mat
    mat
end
function generate_matrix(lista_tablones::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, m::Int, irrep::Array{Int64,1})
    mat = spzeros(Basic,length(lista_tablones), length(lista_tablones))

    for row_idx in 1:length(lista_tablones), col_idx in 1:length(lista_tablones)
        if row_idx == col_idx && mat[row_idx,row_idx] == 0
            mat[row_idx,row_idx] = determinar_coeficiente_irrep_yamanouchi(lista_tablones[row_idx], m)
        elseif row_idx < col_idx
            yamanouchi_block!(mat, lista_tablones, row_idx, col_idx, m, irrep)
        end
    end


    mat
end


@doc Markdown.doc"""
    reposition_out_of_order!(tab::YoungTableau)
> Attempt to fix the first out-of-order entry encountered in a Young tableau
> by bumping it into a corner and renumbering the remaining path.

"""
function reposition_out_of_order!(tab::AbstractAlgebra.Generic.YoungTableau{Int64})
    max_row_seen = 0
    last_label_at_max_row = 0
    violating_label::Int = 0

    traversal_path = Tuple{Int64,Int64}[]

    for label in 1:sum(tab.part)
        row, col = determine_position(tab, label)
        push!(traversal_path,(row,col))
        if row >= max_row_seen
            max_row_seen = row
            last_label_at_max_row = label
        else
            # first label that violates the row-increasing scan; triggers bump
            violating_label = label
            break
        end

    end

    mat = matrix_repr(tab)
    pos_max = determine_position(tab, last_label_at_max_row)
    pos_violation = determine_position(tab, violating_label)

    pos_max = determine_corner(traversal_path)

    mat[pos_max...] = violating_label

    filtered = filter(x -> x != pos_max, sort(traversal_path,by = x->(x[2], x[1])) )
    for (ind, val) in enumerate(filtered)
        if val == pos_max
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
    determine_corner(list of two-integer tuples, e.g. [(1,2),(3,4)])
    was: encontrar_esquina(lista tuplas)
    TODO: see what it actually expects as an input. vector of tuples is a bit weird.
"""
function determine_corner(path::Array{Tuple{Int64,Int64},1})
    (violation_row, violation_col) = path[end]
    corner_row, corner_col = (1, 1)
    if violation_col == 1
        return path[end-1]
    else
       prev_col = violation_col - 1
        while (corner_row <= violation_row && prev_col > 0)
            (corner_row, corner_col) = filter(x -> x[2] == prev_col, path[1:end-1])[end]
            prev_col -= 1
        end
    end
    (corner_row, corner_col)
end

@doc Markdown.doc"""
    yamanouchi_block!(mat, tableaux, row_idx, col_idx, swap_entry, irrep)
> Fill the 2×2 Yamanouchi block for the given pair of tableaux when they differ
> by swapping `swap_entry` and `swap_entry - 1`.

Updates `mat` in place; no effect if the pair is not related by the swap.
"""
function yamanouchi_block!(mat::SparseMatrixCSC{Basic,Int64}, tableaux::Array{AbstractAlgebra.Generic.YoungTableau{Int64},1}, row_idx::Int64, col_idx::Int64, swap_entry::Int64, irrep::Array{Int64,1})
    if row_idx >= col_idx
        return mat
    end

    base_tab = tableaux[row_idx]
    target_tab = tableaux[col_idx]
    if base_tab != swap_adjacent_entries(target_tab, swap_entry, irrep)
        return mat
    end

    axial = Basic(axialdistance(base_tab, swap_entry, swap_entry - 1))
    inv_axial = axial^(-1)
    sqrt_term = sqrt(1 - inv_axial^2)

    mat[row_idx,row_idx] = -inv_axial
    mat[row_idx,col_idx] = sqrt_term
    mat[col_idx,row_idx] = sqrt_term
    mat[col_idx,col_idx] = inv_axial
    mat
end

@doc Markdown.doc"""
    swap_adjacent_entries(tableau::YoungTableau, entry::Int64, irrep::Vector{Int})
> Swap the two occurrences of `entry` and `entry - 1` in a Young tableau, or
> set the lone occurrence to `entry` if only one is present.

Intended for use when building Yamanouchi blocks; assumes `entry > 1` and
`entry` does not exceed the larger of `length(irrep)` and the tableau fill length.

# Examples:
```
julia> t = YoungTableau([2,1]); fill!(t, [1,2,3]);
julia> swap_adjacent_entries(t, 3, [2,1]).fill
3-element Vector{Int64}:
 1
 3
 2

julia> t = YoungTableau([2,1]); fill!(t, [1,1,2]);
julia> swap_adjacent_entries(t, 3, [3]).fill
3-element Vector{Int64}:
 1
 1
 3
```
"""
function swap_adjacent_entries(tableau::AbstractAlgebra.Generic.YoungTableau{Int64}, entry::Int64, irrep::Array{Int64,1})
    fill_entries = copy(tableau.fill)
    @assert entry > 1 && entry <= maximum((length(irrep), length(fill_entries)))

    swap_positions = Int[]
    for (idx, value) in enumerate(fill_entries)
        if value == entry || value == entry - 1
            push!(swap_positions, idx)
        end
        if length(swap_positions) == 2
            break
        end
    end

    if length(swap_positions) > 1
        pos1, pos2 = swap_positions
        fill_entries[pos1], fill_entries[pos2] = fill_entries[pos2], fill_entries[pos1]
    else
        fill_entries[swap_positions[1]] = entry
    end

    new_tableau = YoungTableau(tableau.part)
    fill!(new_tableau, fill_entries)
    new_tableau
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

"""
# TODO add examples
function indice_tablon_semistandard(tablon_semistandard::AbstractAlgebra.Generic.YoungTableau{Int64})
    particion = (tablon_semistandard.part) |> collect
    tablones_standard = StandardYoungTableaux(particion)
    target = standard_tableau_from_semistandard(tablon_semistandard)

    idx = findfirst(x -> x.fill == target.fill, tablones_standard)
    idx === nothing && error("No matching standard tableau found for fill=$(tablon_semistandard.fill) and part=$(tablon_semistandard.part)")
    return idx
end

function generate_dictionary(lista::Array{Int64,1})
    fvars = Dict{Int64,Int64}()
    for (n, f) in enumerate(lista)
        fvars[f] = n
    end
    fvars
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
@doc Markdown.doc"""
    content_length(fill_values::Vector{Int64}, irrep::Irrep) -> Int
> Return the length used for content calculations, i.e. `max(length(fill_values), length(irrep))`.

# Examples:
```
julia> content_length([1,2,3], [2,1])
3

julia> content_length([1,2], [2,1,0,0])
4
```
"""
function content_length(fill_values::Vector{Int64}, irrep::Irrep)
    return max(length(fill_values), length(irrep))
end

@doc Markdown.doc"""
    expanded_fill_values(fill_values::Vector{Int64}, standard_fill::Vector{Int64}, len::Int) -> Vector{Int64}
> Expand `fill_values` by inserting the most recent seen label for missing indices
> in `standard_fill`, up to `len`.

# Examples:
```
julia> expanded_fill_values([1,2,2], [1,3], 4)
4-element Vector{Int64}:
 1
 2
 2
 2
```
"""
function expanded_fill_values(fill_values::Vector{Int64}, standard_fill::Vector{Int64}, len::Int)
    expanded_fill = copy(fill_values)
    last_seen = fill_values[1]

    present = falses(len)
    @inbounds for value in standard_fill
        if value <= len
            present[value] = true
        end
    end

    for i in 1:len
      if !present[i]
        push!(expanded_fill, last_seen)
      end
      if i > length(fill_values)
        break
      end
      if fill_values[i] > last_seen
        last_seen = fill_values[i]
      end
    end

    return expanded_fill
end

@doc Markdown.doc"""
    count_entries(values::Vector{Int64}, len::Int) -> Vector{Int}
> Count occurrences of integers `1:len` in `values`.

# Examples:
```
julia> count_entries([1,1,3,2,3], 3)
3-element Vector{Int64}:
 2
 1
 2
```
"""
function count_entries(values::Vector{Int64}, len::Int)
    counts = zeros(Int, len)
    @inbounds for value in values
        if value <= len
            counts[value] += 1
        end
    end
    return counts
end

@doc Markdown.doc"""
    content(y::YoungTableau, λ::Irrep) -> Vector{Int}
> Return the content vector sized to `max(length(y.fill), length(λ))`.
> In simple terms: it counts how many times each label (1, 2, 3, ...) appears
> in the tableau, padding as needed to match the irrep length.

# Examples:
```
julia> t = YoungTableau([2,1]); fill!(t, [1,2,2]);
julia> content(t, [2,1,0])
3-element Vector{Int64}:
 1
 2
 0

julia> t = YoungTableau([2,1]); fill!(t, [1,1,2]);
julia> content(t, [3])
3-element Vector{Int64}:
 2
 1
 0
```
"""
function content(y::AbstractAlgebra.Generic.YoungTableau{Int64}, λ::Irrep)
    fill_values = y.fill
    len = content_length(fill_values, λ)
    standard_fill = standard_tableau_from_semistandard(y).fill
    expanded = expanded_fill_values(fill_values, standard_fill, len)
    return count_entries(expanded, len)
end

@doc Markdown.doc"""
    content(p::YoungTableau) -> Vector{Int}
> Return the content vector of the tableau fill (counts of each label).
> In simple terms: it counts how many times each label (1, 2, 3, ...) appears
> in the tableau.

# Examples:
```
julia> t = YoungTableau([2,1]); fill!(t, [1,2,2]);
julia> content(t)
2-element Vector{Int64}:
 1
 2
```
"""
function content(p::AbstractAlgebra.Generic.YoungTableau{Int64})
    fill_values = p.fill
    len = length(fill_values)
    return count_entries(fill_values, len)
end

@doc Markdown.doc"""
    standard_tableau_from_semistandard(tablon_semistandard::YoungTableau)
> Returns a YoungTableau corresponding to the standard YoungTableau such that
> f is non decreasing.
"""
function standard_tableau_from_semistandard(tablon_semistandard::AbstractAlgebra.Generic.YoungTableau{Int64})
    particion = (tablon_semistandard.part) |> collect
    tablon_standard = YoungTableau(particion)

    orden = sortperm(tablon_semistandard.fill)  # stable by default
    diccionario = generate_dictionary(orden)
    nuevo_fill = [diccionario[x] for x in tablon_standard.fill]

    tablon_resultado_permutacion = AbstractAlgebra.fill!(tablon_standard, nuevo_fill)
    tablones_standard = StandardYoungTableaux(particion)

    pos = findfirst(x -> x.fill == tablon_resultado_permutacion.fill, tablones_standard)
    pos === nothing && error("No matching standard tableau found for fill=$(tablon_semistandard.fill) and part=$(tablon_semistandard.part)")
    tablones_standard[pos]
end

@doc Markdown.doc"""
    Θ(patron_semi::YoungTableau)
> Computes coefficient Θ. Returns a Float64
"""
function Θ(patron_semi::AbstractAlgebra.Generic.YoungTableau{Int64}, _::Array{Int64,1})
    tablon_standard = standard_tableau_from_semistandard(patron_semi)
    relleno_standard = tablon_standard.fill
    relleno_semi = patron_semi.fill

    # group positions in the standard filling by their semistandard label
    grupos = Dict{Int64,Vector{Int64}}()
    @inbounds for (std, semi) in zip(relleno_standard, relleno_semi)
        push!(get!(grupos, semi, Int[]), std)
    end

    prod = one(Basic)
    for posiciones in values(grupos)
        lenp = length(posiciones)
        lenp <= 1 && continue
        @inbounds for a in 1:lenp-1, b in a+1:lenp
            u = posiciones[a]; v = posiciones[b]
            u, v = u < v ? (u, v) : (v, u)
            dist = axialdistance(tablon_standard, u, v)
            prod *= (one(Basic) + inv(Basic(dist)))
        end
    end
    prod
end

function Θn(patron_semi::AbstractAlgebra.Generic.YoungTableau{Int64}, _::Array{Int64,1})
    tablon_standard = standard_tableau_from_semistandard(patron_semi)
    relleno_standard = tablon_standard.fill
    relleno_semi = patron_semi.fill

    grupos = Dict{Int64,Vector{Int64}}()
    @inbounds for (std, semi) in zip(relleno_standard, relleno_semi)
        push!(get!(grupos, semi, Int[]), std)
    end

    prod = one(Float64)
    for posiciones in values(grupos)
        lenp = length(posiciones)
        lenp <= 1 && continue
        @inbounds for a in 1:lenp-1, b in a+1:lenp
            u = posiciones[a]; v = posiciones[b]
            u, v = u < v ? (u, v) : (v, u)
            dist = axialdistance(tablon_standard, u, v)
            prod *= (1.0 + 1.0 / dist)
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
function calcula_proto_permutacion(proto::AbstractArray{Int64})
    if ndims(proto) == 1
        len = length(proto) |> sqrt |> Int
        new_proto = reshape(proto, len, len)'
    elseif ndims(proto) == 2
        len = size(proto, 1)
        new_proto = proto'
    else
        error("Unsupported input dimension: ", ndims(proto))
    end

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
    standard_to_semistandard_map(semi_tableau::YoungTableau, irrep::Vector{T})
    standard_to_semistandard_map(semi_tableau::YoungTableau)
> Build a dictionary mapping entries of the associated standard tableau
> to the entries of the given semistandard tableau.

# Examples:
```
julia> t_u = YoungTableau([2,2, 1]);
julia> fill!(t_u, [1,2,3,3,4]);
julia> standard_to_semistandard_map(t_u)
Dict(4=>4, 2=>2, 3=>3, 5=>5, 1=>1)
```
"""
function standard_to_semistandard_map(semi_tableau::AbstractAlgebra.Generic.YoungTableau{Int64}, irrep::Array{Int64,1}) 
    len = length(irrep)
    semi_fill = semi_tableau.fill

    length(semi_fill) >= len && return standard_to_semistandard_map(semi_tableau)

    standard_tableau = standard_tableau_from_semistandard(semi_tableau)
    standard_fill = standard_tableau.fill
    
    mapping = Dict{Int64, Int64}()
    @inbounds for (std, semi) in zip(standard_fill, semi_fill)
        mapping[std] = semi
    end

    missing_labels = setdiff!(collect(1:len), standard_fill)
    sorted_standard = sort(standard_fill)
    # Missing-label policy: map unseen labels to the nearest existing standard label by absolute distance.
    # This may send multiple missing labels to the same target; kept for backward compatibility.
    @inbounds for label in missing_labels
        nearest = nearest_standard_label(sorted_standard, label)
        mapping[label] = mapping[nearest]
    end

    mapping
end

@doc Markdown.doc"""
    nearest_standard_label(sorted_standard::Vector{Int64}, label::Int) -> Int
> Return the nearest label in `sorted_standard` to `label` (ties choose the lower neighbor).

# Examples:
```
julia> nearest_standard_label([1,4,7], 6)
7

julia> nearest_standard_label([1,4,7], 5)
4
```
"""
function nearest_standard_label(sorted_standard::Vector{Int64}, label::Int)
    idx = searchsortedlast(sorted_standard, label)
    if idx == 0
        return sorted_standard[1]
    elseif idx == length(sorted_standard)
        return sorted_standard[end]
    else
        left = sorted_standard[idx]
        right = sorted_standard[idx + 1]
        return (label - left) <= (right - label) ? left : right
    end
end

function standard_to_semistandard_map(semi_tableau::AbstractAlgebra.Generic.YoungTableau{Int64})
    semi_fill = semi_tableau.fill
    standard_tableau = standard_tableau_from_semistandard(semi_tableau)
    standard_fill = standard_tableau.fill
    
    mapping = Dict{Int64, Int64}()
    @inbounds for (std, semi) in zip(standard_fill, semi_fill)
      mapping[std] = semi
    end
    mapping
end

@deprecate genera_funcion(semi_tableau::AbstractAlgebra.Generic.YoungTableau{Int64}, irrep::Array{Int64,1}) standard_to_semistandard_map(semi_tableau, irrep)
@deprecate genera_funcion(semi_tableau::AbstractAlgebra.Generic.YoungTableau{Int64}) standard_to_semistandard_map(semi_tableau)

@doc Markdown.doc"""
    stabilizer_permutations(c::Content)
> Return the list of permutations that stabilize each content block (product of symmetric groups).

# Examples:
```
julia> c = [0,1,2]; stabilizer_permutations(c)
```
"""
function stabilizer_permutations(c::Content)
    start_idx = 1
    total_len = length(c)
    perms = Perm{Int64}[]

    for block_size in c
        if block_size == 0
            continue
        elseif block_size == 1
            start_idx += 1
            continue
        end

        block = collect(start_idx:(start_idx + block_size - 1))
        block_perms = collect(permutations(block))

        if isempty(perms)
            for p in block_perms
                push!(perms, Perm(sortperm(vcat(collect(1:start_idx - 1), p, collect((start_idx + block_size):total_len)))))
            end
        else
            new_perms = Perm{Int64}[]
            for p in block_perms
                perm_vec = vcat(collect(1:start_idx - 1), p, collect((start_idx + block_size):total_len))
                perm = Perm(perm_vec)
                for existing in perms
                    push!(new_perms, existing * perm)
                end
            end
            perms = new_perms
        end

        start_idx += block_size
    end

    if isempty(perms)
        return [Perm(collect(1:total_len))]
    end

    unique!(perms)
    perms
end

function stabilizer_permutations(tableau::AbstractAlgebra.Generic.YoungTableau{Int64})
  content_vec = content(tableau)
  stabilizer_permutations(content_vec)
end

function YoungTableau(tab::GTPattern)
    filas = tab.rows
    len = filas[1] |> length
    conjunto_contenido = [calculate_pattern_differences(tab, x) for x in 1:len]
    p = Partition(filter(x -> x>0, filas[1]))
    Generic.YoungTableau(p, vcat(conjunto_contenido...))
end
