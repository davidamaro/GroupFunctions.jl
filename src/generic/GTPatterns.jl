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

  GTPattern(array_of_arrays)
Stucture to hold Gelfand-Tsetlin patterns. 
Data fields:
- rows
- last_row, copied from rows during initialization (for internal algorithm purposes, mathematically irrelevant)
TODO: either make it immutable, or somehow auto-change last_row upon change in rows (and vice versa).
# Example:

```julia
julia> gt=GTPattern([[2,1,0],[2,1],[2]])
│ 2 1 0 ╲
│ 2 1    〉
│ 2     ╱


julia> gt.rows
3-element Vector{Vector{Int64}}:
 [2, 1, 0]
 [2, 1]
 [2]

julia> gt.last_row
1-element Vector{Int64}:
 2
```
"""
mutable struct GTPattern
    rows::Array{Array{Int64,1},1}   #all rows of the 
    last_row::Array{Int64,1}        #cheap solution
end

function GTPattern(rows::Array{Array{Int64,1},1})
  return GTPattern(rows, rows[end])
end

const Row = Array{Int64,1}

function Base.show(io::IO, ::MIME"text/plain", G::GTPattern)
    list_of_rows = G.rows
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

"""
    create_first_pattern!(row::Row, pattern::GTPattern)

Attach the first row to an empty `GTPattern`, recording it as both `rows[1]`
and `last_row` for subsequent construction.
"""
function create_first_pattern!(row::Row, pattern::GTPattern)
    isempty(pattern.last_row) && (pattern.last_row = row)
    push!(pattern.rows, row)
    return pattern
end

"""
    determine_next_possibilities(row::Row)
    Given one row of a GT pattern, calculates all the candidates for the next row.
    TODO: check if the returned items are valid 
"""
function determine_next_possibilities(row::Row)
    # Pre-allocate the array with known size
    n_ranges = length(row) - 1
    ranges = Vector{UnitRange{Int64}}(undef, n_ranges)
    
    # Fill ranges directly without push!
    @inbounds for i in 1:n_ranges
        ranges[i] = row[i+1]:row[i]
    end
    
    return product(ranges...)
end

"""
    determine_next_row(pattern::GTPattern)
    Given a GTPattern, returns a vector of all GTPatterns extended by a single row.
"""
function determine_next_row(pattern::GTPattern)
    # Pre-calculate size for allocation
    next_possibilities = determine_next_possibilities(pattern.last_row)
    patterns = Vector{GTPattern}(undef, length(collect(next_possibilities)))
    
    # Fill patterns array directly
    for (idx, next_row) in enumerate(next_possibilities)
        next_row_vec = collect(next_row)
        new_pattern = GTPattern(
            [pattern.rows..., next_row_vec],
            next_row_vec
        )
        patterns[idx] = new_pattern
    end
    
    return patterns
end

function determine_next_row(patterns::Vector{GTPattern})
    # Calculate total size for pre-allocation
    total_patterns = sum(pattern -> 
        length(collect(determine_next_possibilities(pattern.last_row))), 
        patterns
    )
    
    result_patterns = Vector{GTPattern}(undef, total_patterns)
    current_idx = 1
    
    # Process each pattern
    for pattern in patterns
        next_possibilities = determine_next_possibilities(pattern.last_row)
        
        # Add new patterns directly to pre-allocated array
        for next_row in next_possibilities
            next_row_vec = collect(next_row)
            new_pattern = GTPattern(
                [pattern.rows..., next_row_vec],
                next_row_vec
            )
            result_patterns[current_idx] = new_pattern
            current_idx += 1
        end
    end
    
    return result_patterns
end

function basis_states(weight::Row)
    # Initialize with first pattern
    initial_pattern = GTPattern([], [])
    create_first_pattern!(weight, initial_pattern)
    
    # Generate patterns iteratively
    current_patterns = determine_next_row(initial_pattern)
    
    # Pre-calculate iterations needed
    iterations = length(weight) - 2
    
    # Generate subsequent patterns
    @inbounds for _ in 1:iterations
        current_patterns = determine_next_row(current_patterns)
    end
    
    return current_patterns
end

##############################################################################
#
#   Codigo para la traduccion 
#   (between Young tableaux and GT patterns - Konrad, December 2025)
#
##############################################################################


"""
    calculate_pattern_differences(tab::GTPattern, row_number::Int64)
    TODO: refactorize
    was: obtener_diferencias_patron(tab::GTPattern, numerofila::Int64)
    Returns filling of a Young tableau, given a GT number.
"""
function calculate_pattern_differences(tab::GTPattern, row_number::Int64)
    rows = tab.rows
    number_of_rows = length(rows)
    
    # Input validation with descriptive error
    row_number > number_of_rows && throw(BoundsError("Row number $row_number exceeds pattern size $number_of_rows"))
    
    # Calculate the relevant rows we need to process
    relevant_rows = number_of_rows - row_number + 1
    
    # Pre-allocate the differences array with known size
    differences = zeros(Int64, relevant_rows + 1)
    
    # Process rows in reverse order more efficiently
    max_value = 0
    for (idx, row) in enumerate(view(rows[1:relevant_rows], reverse(1:relevant_rows)))
        current_value = row[row_number]
        current_difference = current_value - max_value
        
        if current_difference > 0
            differences[idx + 1] = current_difference
            max_value = current_value
        end
    end
    
    # Pre-calculate final array size to avoid resizing
    total_elements = sum(differences)
    content = Vector{Int64}(undef, total_elements)
    
    # Fill the content array more efficiently
    pos = 1
    current_value = row_number
    
    # Process all differences except the first (which is always 0)
    for diff_count in view(differences, 2:length(differences))
        if diff_count > 0
            # Use range assignment instead of repeated push!
            range_end = pos + diff_count - 1
            content[pos:range_end] .= current_value
            pos = range_end + 1
        end
        current_value += 1
    end
    
    return content
end

"""
    are_row_sums_decreasing_by_one(tab::GTPattern)
    was: prematuro_pesos
    Returns true if the sums of rows decrease exactly by 1 with each row. false otherwise.
    TODO: for deletion? Not used anywhere, not exported.
    
"""
function are_row_sums_decreasing_by_one(tab::GTPattern)
    number_of_rows = length(tab.rows)
    # Early return for single-row patterns
    number_of_rows == 1 && return true
    
    # Calculate row sums efficiently
    row_sums = Vector{Int}(undef, number_of_rows + 1)
    
    # Fill totals array with row sums
    @inbounds for i in 1:number_of_rows
        row_sums[i] = sum(tab.rows[i])
    end
    row_sums[end] = 0  # Last element is always 0
    
    # Check if differences between consecutive elements are all 1
    # Using direct array indexing instead of creating a new array
    @inbounds for i in 1:number_of_rows-1
        row_sums[i] - row_sums[i+1] != 1 && return false
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
    rows = x.rows
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

"""
    determine_next_pattern!(x::GTPattern)
was: siguientepatron!
Determine next GT pattern by decreasing one entry (if possible) and propagating the change to rows below.
"""
function determine_next_pattern!(x::GTPattern)
    dim_result = decreasable(x)
    isnothing(dim_result) && return nothing
    
    row, col = dim_result
    rows = x.rows
    
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

"""
    determine_next_pattern(x::GTPattern)
was: siguientepatron
Determine next GT pattern by decreasing one entry (if possible) and propagating the change to rows below.
Note: With incomplete patterns e.g. GTPattern([[2,1,0],[2,0]]) outputs of determine_next_pattern and determine_next_pattern! differ.
"""
function determine_next_pattern(x::GTPattern)
    dim_result = decreasable(x)
    isnothing(dim_result) && return nothing
    
    row, col = dim_result
    rows = deepcopy(x.rows)
    
    rows[row][col] -= 1
    
    # Use copyto! for better performance
    copyto!(rows[row], col+1, rows[row-1], col+1, length(rows[row])-col)
    
    # Update subsequent rows efficiently
    for next_row in row+1:length(rows)
        copyto!(rows[next_row], rows[next_row-1])
    end
    
    return GTPattern(rows, rows[end])
end

"""
    decreasable(x::GTPattern)
    was: disminuible
    Determines if it's possible to decrease any entry of a GT pattern by 1. 
    If not, does not return nothing.
    If yes, returns (row index, column index) of the first found solution.
    Examples:
    ```julia
    julia> gt=GTPattern([[2,1,0],[2,0]])
    │ 2 1 0 ╲
    │ 2 0   ╱

    julia> GroupFunctions.decreasable(gt)
    (2, 1)

    julia> gt=GTPattern([[2,1,0],[1,0]])
    │ 2 1 0 ╲
    │ 1 0   ╱

    julia> GroupFunctions.decreasable(gt)
    ```

"""
function decreasable(x::GTPattern)
    rows = deepcopy(x.rows)
    
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

function initial_gt(irrep::Row)
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
