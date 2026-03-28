import Base: fill!, show, ==

"""
    YoungTableau(part::AbstractVector{<:Integer})
    YoungTableau(part::AbstractVector{<:Integer}, fill::AbstractVector{<:Integer})

Minimal Young tableau type used by GroupFunctions. The `part` argument is the
partition describing the tableau shape by row lengths.
"""
mutable struct YoungTableau
    part::Vector{Int64}
    fill::Vector{Int64}
end

function normalize_partition(part::AbstractVector{<:Integer})
    normalized = filter(>(0), Int64.(part))
    issorted(normalized; rev = true) || throw(ArgumentError("Partition must be non-increasing"))
    normalized
end

function YoungTableau(part::AbstractVector{<:Integer})
    normalized = normalize_partition(part)
    YoungTableau(normalized, collect(Int64, 1:sum(normalized)))
end

function YoungTableau(part::AbstractVector{<:Integer}, fill_values::AbstractVector{<:Integer})
    normalized = normalize_partition(part)
    fill_int = Int64.(fill_values)
    expected = sum(normalized)
    length(fill_int) == expected || throw(ArgumentError("Expected $expected fill values, got $(length(fill_int))"))
    YoungTableau(normalized, fill_int)
end

function fill!(tableau::YoungTableau, values::AbstractVector{<:Integer})
    expected = sum(tableau.part)
    length(values) == expected || throw(ArgumentError("Expected $expected fill values, got $(length(values))"))
    tableau.fill = Int64.(values)
    tableau
end

==(lhs::YoungTableau, rhs::YoungTableau) = lhs.part == rhs.part && lhs.fill == rhs.fill

function tableau_rows(tableau::YoungTableau)
    max_cols = isempty(tableau.part) ? 0 : maximum(tableau.part)
    rows = Vector{Vector{Int64}}()
    idx = 1
    for row_len in tableau.part
        row = zeros(Int64, max_cols)
        for col_idx in 1:row_len
            row[col_idx] = tableau.fill[idx]
            idx += 1
        end
        push!(rows, row)
    end
    rows
end

function show(io::IO, tableau::YoungTableau)
    rows = tableau_rows(tableau)
    print(io, "[")
    for (row_idx, row) in enumerate(rows)
        row_idx > 1 && print(io, "; ")
        print(io, join(row, " "))
    end
    print(io, "]")
end

show(io::IO, ::MIME"text/plain", tableau::YoungTableau) = show(io, tableau)

function matrix_repr(tableau::YoungTableau)
    nrows = length(tableau.part)
    ncols = isempty(tableau.part) ? 0 : maximum(tableau.part)
    matrix = spzeros(Int64, nrows, ncols)
    idx = 1
    for (row_idx, row_len) in enumerate(tableau.part)
        for col_idx in 1:row_len
            matrix[row_idx, col_idx] = tableau.fill[idx]
            idx += 1
        end
    end
    matrix
end

function hook_length(part::Vector{Int64}, row_idx::Int, col_idx::Int)
    right = part[row_idx] - col_idx
    below = 0
    for next_row in (row_idx + 1):length(part)
        below += part[next_row] >= col_idx
    end
    right + below + 1
end

function dim(part::AbstractVector{<:Integer})
    normalized = normalize_partition(part)
    n = sum(normalized)
    n == 0 && return 1

    hooks = big(1)
    for (row_idx, row_len) in enumerate(normalized), col_idx in 1:row_len
        hooks *= hook_length(normalized, row_idx, col_idx)
    end

    Int(factorial(big(n)) ÷ hooks)
end

dim(tableau::YoungTableau) = dim(tableau.part)
