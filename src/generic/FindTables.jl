module FindTables

export find_tablaeux_fillings

using ..AllSolutionsMatrix

function rows_to_matrix(rows::Vector{Vector{T}}) where {T<:Integer}
    nrows = length(rows)
    nrows == 0 && return Matrix{T}(undef, 0, 0)
    ncols = length(rows[1])
    matrix = Matrix{T}(undef, nrows, ncols)
    @inbounds for i in 1:nrows
        row = rows[i]
        @inbounds for j in 1:ncols
            matrix[i, j] = row[j]
        end
    end
    return matrix
end

"""
    find_tablaeux_fillings(A::Vector{Int}, B::Vector{Int})
    was: encontrar_prototablones
    see the following discussion for context:
        https://discourse.julialang.org/t/right-solver-for-jump-to-find-every-solution-of-a-linear-system-of-equations-with-integer-solutions/44709/6
"""
function find_tablaeux_fillings(A::Vector{Int}, B::Vector{Int})
    if length(A) != length(B) || sum(A) != sum(B)
        println("No solution exists: row sums must equal column sums")
        return Matrix{Int}[]
    end

    return map(rows_to_matrix, enumerate_matrices(A, B))
end

end  # module
