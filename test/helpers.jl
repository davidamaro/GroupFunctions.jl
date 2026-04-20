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

function es_cero(pol::Basic; ϵ = 10^(-5))
    monomios_lista = SymEngine.free_symbols(pol)

    while length(monomios_lista) > 0
      mono = pop!(monomios_lista)
      pol = SymEngine.subs(pol, mono, randn())
    end

    abs(pol)^2 < ϵ
end
