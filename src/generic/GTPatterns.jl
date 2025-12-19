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
│ 2   1   0 ╲
│   2   1    〉
│     2     ╱


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

"""
    Base.show(io::IO, ::MIME"text/plain", G::GTPattern)

Pretty-print a `GTPattern` as a triangular diagram with slashes indicating
the GT interlacing structure.
"""
function Base.show(io::IO, ::MIME"text/plain", G::GTPattern)
    list_of_rows = G.rows
    n = length(list_of_rows)
    num_width = maximum(length(string(num)) for row in list_of_rows for num in row)
    spacing = num_width + 2                       # keep entries apart and guarantees even padding
    first_row_len = length(list_of_rows[1])
    first_row_width = first_row_len * num_width + (first_row_len - 1) * spacing
    conn_upper = '╲'
    conn_middle = '〉'
    conn_lower = '╱'

    mitad = n ÷ 2
    between = isodd(n)
    cont = 1
    i = 1

    buffer = IOBuffer()
    for row in list_of_rows
        row_width = length(row) * num_width + (length(row) - 1) * spacing
        padding = (first_row_width - row_width) ÷ 2

        print(buffer, "│ ", repeat(" ", padding))
        for (j, number) in enumerate(row)
            num_str = string(number)
            print(buffer, rpad(num_str, num_width))
            j < length(row) && print(buffer, repeat(" ", spacing))
        end

        print(buffer, repeat(" ", padding + 1))
        if i <= mitad
            print(buffer, repeat(" ", cont - 1), conn_upper, "\n")
            cont += 1
        elseif between
            print(buffer, repeat(" ", cont - 1), conn_middle, "\n")
            between = false
        else
            cont > mitad && (cont -= 1)
            print(buffer, repeat(" ", cont - 1), conn_lower, "\n")
            cont -= 1
        end
        i += 1
    end

    print(io, String(take!(buffer)))
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

struct NextRowIterator
    ranges::Vector{UnitRange{Int}}
    len::Int
end

Base.IteratorSize(::Type{NextRowIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{NextRowIterator}) = Base.HasEltype()
Base.eltype(::Type{NextRowIterator}) = Row
Base.length(it::NextRowIterator) = it.len

@inline function Base.iterate(it::NextRowIterator)
    it.len == 0 && return nothing
    n = length(it.ranges)
    if n == 0
        return (Int[], nothing)
    end
    current = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        current[i] = first(it.ranges[i])
    end
    return (copy(current), current)
end

@inline function Base.iterate(it::NextRowIterator, state::Union{Nothing, Vector{Int}})
    current = state
    current === nothing && return nothing

    n = length(it.ranges)
    @inbounds for i in n:-1:1
        rng = it.ranges[i]
        if current[i] < last(rng)
            current[i] += 1
            # reset lower positions
            for j in i+1:n
                current[j] = first(it.ranges[j])
            end
            return (copy(current), current)
        end
    end

    return nothing
end

"""
    determine_next_possibilities(row::Row)
    Given one row of a GT pattern, calculates all the candidates for the next row.
"""
function determine_next_possibilities(row::Row)
    n_ranges = length(row) - 1
    ranges = Vector{UnitRange{Int}}(undef, n_ranges)

    @inbounds for i in 1:n_ranges
        ranges[i] = row[i+1]:row[i]
    end

    len = n_ranges == 0 ? 1 : begin
        prod = 1
        @inbounds for r in ranges
            rlen = length(r)
            prod *= rlen
            rlen == 0 && break
        end
        prod
    end

    return NextRowIterator(ranges, len)
end

"""
    determine_next_row(pattern::GTPattern)
    Given a GTPattern, returns a vector of all GTPatterns extended by a single row.
"""
function determine_next_row(pattern::GTPattern)
    next_possibilities = determine_next_possibilities(pattern.last_row)
    count = length(next_possibilities)
    patterns = Vector{GTPattern}(undef, count)

    rows_prefix = pattern.rows
    prefix_len = length(rows_prefix)

    idx = 1
    for next_row_vec in next_possibilities
        new_rows = Vector{Row}(undef, prefix_len + 1)
        copyto!(new_rows, 1, rows_prefix, 1, prefix_len)
        new_rows[end] = next_row_vec
        patterns[idx] = GTPattern(new_rows, next_row_vec)
        idx += 1
    end

    return patterns
end

"""
    determine_next_row(patterns::Vector{GTPattern})

Given a collection of GT patterns, extend each by one row in all valid ways.
Returns a flat vector containing every resulting `GTPattern`.
"""
function determine_next_row(patterns::Vector{GTPattern})
    total_patterns = sum(pattern -> length(determine_next_possibilities(pattern.last_row)), patterns)
    result_patterns = Vector{GTPattern}(undef, total_patterns)
    current_idx = 1
    
    # Process each pattern
    for pattern in patterns
        next_possibilities = determine_next_possibilities(pattern.last_row)
        rows_prefix = pattern.rows
        prefix_len = length(rows_prefix)
        
        # Add new patterns directly to pre-allocated array
        for next_row_vec in next_possibilities
            new_rows = Vector{Row}(undef, prefix_len + 1)
            copyto!(new_rows, 1, rows_prefix, 1, prefix_len)
            new_rows[end] = next_row_vec
            result_patterns[current_idx] = GTPattern(new_rows, next_row_vec)
            current_idx += 1
        end
    end
    
    return result_patterns
end

function basis_states(weight::Row)
    isempty(weight) && return GTPattern[]

    initial_pattern = GTPattern(Row[weight], weight)
    current_patterns = determine_next_row(initial_pattern)

    iterations = max(length(weight) - 2, 0)
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
    row_num, n_rows = _validate_row_number(row_number, rows)
    return _calculate_pattern_differences(rows, row_num, n_rows)
end

@inline function _validate_row_number(row_number::Integer, rows::Vector{Row})
    row_number < 1 && throw(BoundsError("row_number must be at least 1, got $row_number"))
    n_rows = length(rows)
    row_number > n_rows && throw(BoundsError("Row number $row_number exceeds pattern size $n_rows"))
    return Int(row_number), n_rows
end

@inline function _calculate_pattern_differences(rows::Vector{Row}, row_number::Int, n_rows::Int)
    relevant_rows = n_rows - row_number + 1
    differences = Vector{Int}(undef, relevant_rows)

    # Walk the relevant rows from bottom to top, accumulating positive jumps
    max_value = 0
    @inbounds for i in 1:relevant_rows
        row = rows[relevant_rows - i + 1]
        val = row[row_number]
        diff = val - max_value
        if diff > 0
            differences[i] = diff
            max_value = val
        else
            differences[i] = 0
        end
    end

    total_elements = sum(differences)
    content = Vector{Int}(undef, total_elements)

    pos = 1
    current_value = row_number
    @inbounds for diff_count in differences
        if diff_count > 0
            for j in 0:diff_count-1
                content[pos + j] = current_value
            end
            pos += diff_count
        end
        current_value += 1
    end

    return content
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
    rows = x.rows
    n_rows = length(rows)

    # Scan from bottom to top, checking GT inequalities in-place
    @inbounds for row_idx in n_rows:-1:2
        upper = rows[row_idx - 1]
        current = rows[row_idx]
        for col_idx in eachindex(current)
            val = current[col_idx]
            val == 0 && continue

            val -= 1
            # Check GT constraints locally without allocations
            left_ok = col_idx == 1 || val <= upper[col_idx - 1]
            right_ok = val >= upper[col_idx]
            if left_ok && right_ok
                return row_idx, col_idx
            end
        end
    end

    return nothing
end
