#export GTPattern
#export basis_states, obtener_diferencias_patron#, prematuro_pesos#, YoungTableau
import Base.isvalid


using Base.Iterators

##############################################################################
#
#   Partition type, AbstractVector interface
#
##############################################################################


@doc Markdown.doc"""

  GTPattern(Array of arrays, Array)

# Examples:

```julia
julia> GTPattern([[2,1,0],[2,1],[2]],[2])
```
"""
mutable struct GTPattern
    filas::Array{Array{Int64,1},1}
    ultima_fila::Array{Int64,1}
end
const Fila = Array{Int64,1}

function Base.show(io::IO, ::MIME"text/plain", G::GTPattern)
    list_of_rows = G.filas
    pattern = ""
    n = length(list_of_rows)
    mitad = n ÷ 2
    between = isodd(n)
    cont = 1
    i = 1
    length_first_row = length(list_of_rows[1])*2
    for row in list_of_rows
        pattern *= "│ "
        for number in row
            pattern *= string(number)
            pattern *= " "
        end
        pattern *= repeat(" ", length_first_row - 2*length(row))
        if i <= mitad
            pattern *= repeat(" ", cont - 1)
            pattern *= "╲\n"
            cont += 1
        elseif between
            pattern *= repeat(" ", cont - 1)
            pattern *= "〉\n"
            between = false
        else
            if !between && cont > mitad
                cont -= 1
            end
            pattern *= repeat(" ", cont - 1)
            pattern *= "╱\n"
            cont -= 1
        end
        i += 1
    end

    print(io, pattern)
end

function primera!(row::Fila, pattern::GTPattern)
    isempty(pattern.ultima_fila) && (pattern.ultima_fila = row)
    push!(pattern.filas, row)
    return pattern
end

function determinar_siguientes(row::Fila)
    # Pre-allocate the array with known size
    n_ranges = length(row) - 1
    ranges = Vector{UnitRange{Int64}}(undef, n_ranges)
    
    # Fill ranges directly without push!
    @inbounds for i in 1:n_ranges
        ranges[i] = row[i+1]:row[i]
    end
    
    return product(ranges...)
end

function generar_siguiente_fila(pattern::GTPattern)
    # Pre-calculate size for allocation
    next_possibilities = determinar_siguientes(pattern.ultima_fila)
    patterns = Vector{GTPattern}(undef, length(collect(next_possibilities)))
    
    # Fill patterns array directly
    for (idx, next_row) in enumerate(next_possibilities)
        next_row_vec = collect(next_row)
        new_pattern = GTPattern(
            [pattern.filas..., next_row_vec],
            next_row_vec
        )
        patterns[idx] = new_pattern
    end
    
    return patterns
end

function generar_siguiente_fila(patterns::Vector{GTPattern})
    # Calculate total size for pre-allocation
    total_patterns = sum(pattern -> 
        length(collect(determinar_siguientes(pattern.ultima_fila))), 
        patterns
    )
    
    result_patterns = Vector{GTPattern}(undef, total_patterns)
    current_idx = 1
    
    # Process each pattern
    for pattern in patterns
        next_possibilities = determinar_siguientes(pattern.ultima_fila)
        
        # Add new patterns directly to pre-allocated array
        for next_row in next_possibilities
            next_row_vec = collect(next_row)
            new_pattern = GTPattern(
                [pattern.filas..., next_row_vec],
                next_row_vec
            )
            result_patterns[current_idx] = new_pattern
            current_idx += 1
        end
    end
    
    return result_patterns
end

function basis_states(weight::Fila)
    # Initialize with first pattern
    initial_pattern = GTPattern([], [])
    primera!(weight, initial_pattern)
    
    # Generate patterns iteratively
    current_patterns = generar_siguiente_fila(initial_pattern)
    
    # Pre-calculate iterations needed
    iterations = length(weight) - 2
    
    # Generate subsequent patterns
    @inbounds for _ in 1:iterations
        current_patterns = generar_siguiente_fila(current_patterns)
    end
    
    return current_patterns
end

##############################################################################
#
#   Codigo para la traduccion
#
##############################################################################

function obtener_diferencias_patron(tab::GTPattern, numerofila::Int64)
    filas = tab.filas
    n_filas = length(filas)
    
    # Input validation with descriptive error
    numerofila > n_filas && throw(BoundsError("Row number $numerofila exceeds pattern size $n_filas"))
    
    # Calculate the relevant rows we need to process
    relevant_rows = n_filas - numerofila + 1
    
    # Pre-allocate the differences array with known size
    diferencias = zeros(Int64, relevant_rows + 1)
    
    # Process rows in reverse order more efficiently
    max_value = 0
    for (idx, fila) in enumerate(view(filas[1:relevant_rows], reverse(1:relevant_rows)))
        current_value = fila[numerofila]
        diferencia = current_value - max_value
        
        if diferencia > 0
            diferencias[idx + 1] = diferencia
            max_value = current_value
        end
    end
    
    # Pre-calculate final array size to avoid resizing
    total_elements = sum(diferencias)
    contenido = Vector{Int64}(undef, total_elements)
    
    # Fill the content array more efficiently
    pos = 1
    current_value = numerofila
    
    # Process all differences except the first (which is always 0)
    for diff_count in view(diferencias, 2:length(diferencias))
        if diff_count > 0
            # Use range assignment instead of repeated push!
            range_end = pos + diff_count - 1
            contenido[pos:range_end] .= current_value
            pos = range_end + 1
        end
        current_value += 1
    end
    
    return contenido
end

function prematuro_pesos(tab::GTPattern)
    n_filas = length(tab.filas)
    
    # Early return for single-row patterns
    n_filas == 1 && return true
    
    # Calculate row sums efficiently
    totales = Vector{Int}(undef, n_filas + 1)
    
    # Fill totals array with row sums
    @inbounds for i in 1:n_filas
        totales[i] = sum(tab.filas[i])
    end
    totales[end] = 0  # Last element is always 0
    
    # Check if differences between consecutive elements are all 1
    # Using direct array indexing instead of creating a new array
    @inbounds for i in 1:n_filas-1
        totales[i] - totales[i+1] != 1 && return false
    end
    
    return true
end

@doc Markdown.doc"""
> Custom `isvalid` for GTPattern

  isvalid(x::GTPattern)

# Examples:

```julia
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2])
julia> isvalid(t)
>true


julia> t = GTPattern([[2,1,0],[2,2],[2]],[2])
julia> isvalid(t)
>false
```
"""
function isvalid(x::GTPattern)
    rows = x.filas
    n_rows = length(rows)
    
    # Early return for single-row patterns
    n_rows == 1 && return true
    
    # Vectorized comparison instead of nested loops
    for i in 1:n_rows-1
        upper = rows[i]
        lower = rows[i+1]
        
        # Check if any element violates the GT pattern conditions
        # Using broadcasting for element-wise comparison
        if any(!(lower[j] <= upper[j] && lower[j] >= upper[j+1]) 
               for j in eachindex(lower))
            return false
        end
    end
    return true
end

function siguientepatron!(x::GTPattern)
    dim_result = disminuible(x)
    isnothing(dim_result) && return nothing
    
    row, col = dim_result
    rows = x.filas
    
    # Decrease the value at the found position
    rows[row][col] -= 1
    
    # Update subsequent elements in the current row
    copyto!(rows[row], 1, rows[row-1], 1, col-1)
    
    # Update subsequent rows
    for next_row in row+1:length(rows)
        copyto!(rows[next_row], rows[next_row-1])
    end
    
    return x
end

function siguientepatron(x::GTPattern)
    dim_result = disminuible(x)
    isnothing(dim_result) && return nothing
    
    row, col = dim_result
    rows = deepcopy(x.filas)
    
    rows[row][col] -= 1
    
    # Use copyto! for better performance
    copyto!(rows[row], col+1, rows[row-1], col+1, length(rows[row])-col)
    
    # Update subsequent rows efficiently
    for next_row in row+1:length(rows)
        copyto!(rows[next_row], rows[next_row-1])
    end
    
    return GTPattern(rows, rows[end])
end

function disminuible(x::GTPattern)
    rows = deepcopy(x.filas)
    
    # Iterate from bottom to top
    for row in length(rows):-1:2
        row_vals = rows[row]
        
        # Check each position in the current row
        for col in eachindex(row_vals)
            row_vals[col] == 0 && continue
            
            # Try decreasing the current value
            row_vals[col] -= 1
            
            # Check if the resulting pattern is valid
            if isvalid(GTPattern(rows, rows[end]))
                return row, col
            end
            
            # Restore the value if invalid
            row_vals[col] += 1
        end
    end
    
    return nothing
end

function gtinicial(irrep::Fila)
    n = length(irrep)
    # Pre-allocate array with known size
    rows = Vector{Vector{Int64}}(undef, n)
    rows[1] = irrep
    
    # Build subsequent rows more efficiently
    for i in 2:n
        rows[i] = rows[i-1][1:end-1]
    end
    
    return GTPattern(rows, rows[end])
end
