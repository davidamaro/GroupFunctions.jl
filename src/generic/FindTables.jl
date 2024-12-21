module FindTables

using JuMP
using HiGHS

export encontrar_prototablones

# Paste the function code here
"""
Find all possible solutions to the matrix balancing problem
"""
function encontrar_prototablones(A::Vector{Int}, B::Vector{Int})
    if length(A) != length(B) || sum(A) != sum(B)
        println("No solution exists: row sums must equal column sums")
        return []
    end
    
    n = length(A)
    solutions = []
    
    function solve_with_fixed_values(fixed_values)
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        
        # Variables
        @variable(model, M[1:n, 1:n] >= 0, Int)
        
        # Apply fixed values
        for ((i,j), val) in fixed_values
            @constraint(model, M[i,j] == val)
        end
        
        # Basic constraints
        for i in 1:n
            @constraint(model, sum(M[i,:]) == A[i])
        end
        for j in 1:n
            @constraint(model, sum(M[:,j]) == B[j])
        end
        
        optimize!(model)
        
        if termination_status(model) == OPTIMAL
            return round.(Int, value.(M))
        else
            return nothing
        end
    end
    
    function explore_solutions(fixed_values=Dict())
        sol = solve_with_fixed_values(fixed_values)
        if sol === nothing
            return
        end
        
        # If we have a complete solution, add it
        if length(fixed_values) == n*n-1
            if !any(all(s .== sol) for s in solutions)
                push!(solutions, copy(sol))
            end
            return
        end
        
        # Find next unfixed position
        next_i, next_j = 1, 1
        for i in 1:n, j in 1:n
            if !haskey(fixed_values, (i,j))
                next_i, next_j = i, j
                break
            end
        end
        
        # Try different values for the next position
        max_val = min(A[next_i], B[next_j])
        for val in 0:max_val
            new_fixed = copy(fixed_values)
            new_fixed[(next_i, next_j)] = val
            explore_solutions(new_fixed)
        end
    end
    
    explore_solutions()
    return solutions
end

end  # module
