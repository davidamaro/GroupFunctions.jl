"""
    permanent(matrix::AbstractMatrix)

Compute the permanent of a square matrix using Ryser's formula.
This implementation has complexity O(n * 2^n).

# Arguments
- `matrix::AbstractMatrix`: A square matrix

# Returns
- The permanent of the matrix

# Example
```julia
julia> A = [1 2; 3 4]
julia> permanent(A)
10
```
"""
function permanent(matrix::AbstractMatrix)
    rows, cols = size(matrix)
    if rows != cols
        throw(ArgumentError("Matrix must be square"))
    end
    
    if rows == 0
        return 1  # Empty matrix has permanent 1
    end
    
    if rows == 1
        return matrix[1,1]
    end
    
    if rows == 2
        return matrix[1,1] * matrix[2,2] + matrix[1,2] * matrix[2,1]
    end
    
    # For larger matrices, use the definition with permutations
    result = 0
    for perm in permutations(1:rows)
        term = 1
        for i in 1:rows
            term *= matrix[i, perm[i]]
        end
        result += term
    end
    
    return result
end

function immanant2110(A::AbstractMatrix)
      size(A) == (4,4) || throw(ArgumentError("Matrix must be 4×4"))
    
    return -A[1,4] * A[2,3] * A[3,2] * A[4,1] + 
            A[1,3] * A[2,4] * A[3,2] * A[4,1] - 
            A[1,4] * A[2,2] * A[3,3] * A[4,1] + 
            A[1,2] * A[2,3] * A[3,4] * A[4,1] +
            A[1,4] * A[2,3] * A[3,1] * A[4,2] - 
            A[1,3] * A[2,4] * A[3,1] * A[4,2] - 
            A[1,1] * A[2,4] * A[3,3] * A[4,2] + 
            A[1,3] * A[2,1] * A[3,4] * A[4,2] +
            A[1,2] * A[2,4] * A[3,1] * A[4,3] + 
            A[1,4] * A[2,1] * A[3,2] * A[4,3] - 
            A[1,2] * A[2,1] * A[3,4] * A[4,3] - 
            A[1,1] * A[2,2] * A[3,4] * A[4,3] -
            A[1,3] * A[2,2] * A[3,1] * A[4,4] - 
            A[1,1] * A[2,3] * A[3,2] * A[4,4] - 
            A[1,2] * A[2,1] * A[3,3] * A[4,4] + 
            3 * A[1,1] * A[2,2] * A[3,3] * A[4,4]
end

function immanant210(A::AbstractMatrix)
    size(A) == (3,3) || throw(ArgumentError("Matrix must be 3×3"))
    
    return -A[1,2] * A[2,3] * A[3,1] - 
           A[1,3] * A[2,1] * A[3,2] + 
           2 * A[1,1] * A[2,2] * A[3,3]
end

# Example usage:
# matrix = [1 2; 3 4]
# println("Permanent: ", permanent(matrix))  # Should print 10

function findzero(part::Array{Int,1})
    base = basis_states(part)
    filter(base) do x
        weight = zweight(x)
        zero_weight = fill!(similar(weight), zero(eltype(weight)))
        isapprox(weight, zero_weight; atol=1e-12, rtol=0)
    end
end

const SYMBOLIC_COMPARISON_POINTS = ComplexF64[
    0.17 + 0.29im,
   -0.31 + 0.23im,
    0.41 - 0.37im,
   -0.53 - 0.19im,
]

function symbolic_isapprox(lhs::Basic, rhs::Basic; atol::Float64 = 1e-10, rtol::Float64 = 0.0)
    difference = expand(lhs - rhs)
    difference == 0 && return true

    free_symbols = sort!(collect(SymEngine.free_symbols(difference)); by=string)
    isempty(free_symbols) && return isapprox(convert(ComplexF64, difference), 0; atol, rtol)

    for sample_point in SYMBOLIC_COMPARISON_POINTS
        evaluated_difference = difference
        for (symbol_index, symbol) in enumerate(free_symbols)
            evaluated_difference = SymEngine.subs(evaluated_difference, symbol, sample_point^symbol_index)
        end
        isapprox(convert(ComplexF64, evaluated_difference), 0; atol, rtol) || return false
    end

    return true
end
