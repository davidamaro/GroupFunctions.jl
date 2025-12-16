module FindTables

using JuMP
using HiGHS

export find_tablaeux_fillings

mutable struct FixedValues
    positions::Matrix{Int}
    values::Vector{Int}  # Changed to Vector since we only need 1D for values
    count::Int
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
        return []
    end
    
    n = length(A)
    solutions = []
    
    function solve_with_fixed_values(fixed_values)
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        
        max_possible = max(maximum(A), maximum(B))
        @variable(model, 0 <= M[1:n, 1:n] <= max_possible, Int)
        # @variable(model, M[1:n, 1:n] >= 0, Int)
        
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
    
    function init_fixed_values()
        # n*n positions possible in total
        FixedValues(zeros(Int, n*n, 2), zeros(Int, n*n), 0)
    end
    
    function add_value!(fv::FixedValues, i::Int, j::Int, val::Int)
        fv.count += 1
        fv.positions[fv.count, 1] = i
        fv.positions[fv.count, 2] = j
        fv.values[fv.count] = val
    end
    
    function remove_last!(fv::FixedValues)
        fv.count -= 1
    end
    
    function solution_hash(sol::Matrix{Int})
        hash(vec(sol))
    end
    
    solution_hashes = Set{UInt}()
    
    function explore_solutions(fixed::FixedValues)
        # Convert fixed values to format needed by solver
        fixed_dict = Dict((fixed.positions[i,1], fixed.positions[i,2]) => fixed.values[i] 
                         for i in 1:fixed.count)
        
        sol = solve_with_fixed_values(fixed_dict)
        if sol === nothing
            return
        end
        
        # If we have a complete solution, add it if unique
        if fixed.count == n*n-1
            sol_hash = solution_hash(sol)
            if !(sol_hash in solution_hashes)
                push!(solution_hashes, sol_hash)
                push!(solutions, copy(sol))
            end
            return
        end
        
        # Find next unfixed position
        next_i, next_j = 1, 1
        found = false
        for i in 1:n, j in 1:n
            if !any(k -> fixed.positions[k,1] == i && fixed.positions[k,2] == j, 1:fixed.count)
                next_i, next_j = i, j
                found = true
                break
            end
        end
        
        if !found
            return
        end
        
        # Try different values for the next position
        max_val = min(A[next_i], B[next_j])
        for val in 0:max_val
            add_value!(fixed, next_i, next_j, val)
            explore_solutions(fixed)
            remove_last!(fixed)
        end
    end
    
    explore_solutions(init_fixed_values())
    return solutions
end

end  # module
