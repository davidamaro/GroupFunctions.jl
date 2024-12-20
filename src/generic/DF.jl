import Combinatorics: permutations
# using IntervalArithmetic
# import IntervalConstraintProgramming: SubPaving, Contractor
# import ModelingToolkit: @variables, Variable
using JuMP
# using Gurobi
using HiGHS
#import SArray
import LinearAlgebra: transpose
import StaticArrays: SArray
export group_function, mma_to_julia, @mma
export zweight, pweight
export encontrar_prototablones

const Irrep = Array{T,1} where T <: Integer
const YTableau = AbstractAlgebra.Generic.YoungTableau{T} where T <: Integer
const Content = Array{T,1} where T <: Integer
const MapST2SST = Dict{T,T} where T <: Integer

# function pave(Cs, working::Vector{IntervalBox{N,T}}, ϵ, bisection_point=nothing) where {N,T}

#     boundary = SubPaving{N,T}()

#     while !isempty(working)

#         X = pop!(working)

#         # apply each contractor in a round-robin fashion:
#         for C in Cs
#             X = C(X)
            
#             if isempty(X)
#                 break
#             end
#         end
        
#         if isempty(X)
#             continue
#         end

#         if diam(X) < ϵ
#             push!(boundary, X)

#         else
#             if isnothing(bisection_point)
#                 push!(working, bisect(X)...)
#             else
#                 push!(working, bisect(X, bisection_point)...)
#             end

#         end

#     end

#     return boundary

# end


# function integerize(X::Interval)
#     a = ceil(X.lo)
#     b = floor(X.hi)
    
#     if a > b
#         return emptyinterval(X)
#     end

#     return Interval(a, b)
# end 

# integerize(X::IntervalBox) = integerize.(X)

# @doc Markdown.doc"""
# > Return the size of the vector which represents the partition.

# encontrar_prototablones(Array, Array)

# # Examples:

# ```
# julia > encontrar_prototablones([1,1,1], [1,2,0])
#  [1, 0, 0, 0, 1, 1, 0, 0, 0]
#  [0, 1, 0, 1, 0, 1, 0, 0, 0]
#  [0, 0, 1, 1, 1, 0, 0, 0, 0]
# ```
# """
# function encontrar_prototablones(μ::Array{Int64,1}, ν::Array{Int64,1})

#     # intervalo y contractors
#     n::Int64 = length(μ)
#     x::IntervalBox{n^2, Float64}, c::Array{Contractor,1} = generarcontractors(μ,ν)

#     contractors::Array{Function,1} = [X -> C(0..0, X) for C in c]

#     helper = pop!(contractors)
#     X::IntervalBox{n^2, Float64} = Base.invokelatest(helper, x)

#     for C in contractors
#         X = Base.invokelatest(C, X)
#     end

#     contractors = [contractors; integerize]
#     solutions::Array{IntervalBox{n^2,Float64},1} = Base.invokelatest(pave, contractors,[X], 1.0)


#     SArray{Tuple{n^2},Int64,1,n^2}[Int.(x) for x in mid.(solutions)]
# end

# function generarcontractors(μ::Array{Int64,1}, ν::Array{Int64,1})
#     n::Int64 = length(μ)

#     vars = (@variables y[1:n, 1:n])[1]

#     contractors::Array{Contractor,1} = Contractor[]

#     for i in 1:n
#         push!(contractors, Contractor(vars, sum(vars[i, 1:n]) - μ[i]))
#     end

#     for i in 1:n
#         push!(contractors, Contractor(vars, sum(vars[1:n, i]) - ν[i]))
#     end

#     X::IntervalBox = IntervalBox(0..n^2, n^2)
#     X, contractors
# end
function flatten(x)
    n, _ = x |> size
    reshape(x, (n^2,1))
#     reshape(x, (1,n^2))
end


#
# using JuMP
# using HiGHS

# function find_all_matrices(rowsum::Array{Int64,1}, colsum::Array{Int64,1}; verbose::Bool = false, num_sol::Int = 1000)
    # # Input validation
    # @assert length(rowsum) == length(colsum) "Row and column sums must have equal length"
    # @assert sum(rowsum) == sum(colsum) "Sum of row sums must equal sum of column sums"
    
    # N = length(rowsum)
    # INDEX = 1:N
    # Results = Array{Array{Int64,2},1}()
    
    # # Create initial model
    # function create_base_model()
        # model = Model(HiGHS.Optimizer)
        # if !verbose
            # set_silent(model)
        # end
        
        # # Variables and basic constraints
        # @variable(model, y[INDEX,INDEX] >= 0, Int)
        # @constraint(model, rowCons[i=INDEX], sum(y[i,j] for j in INDEX) == rowsum[i])
        # @constraint(model, colCons[j=INDEX], sum(y[i,j] for i in INDEX) == colsum[j])
        
        # return model, y
    # end
    
    # # Function to add constraint excluding previous solution
    # function add_exclusion_constraint(model, y, prev_sol)
        # # Sum of absolute differences must be at least 1
        # @constraint(model, sum(y[i,j] for i in INDEX, j in INDEX where prev_sol[i,j] == 1) + 
                          # sum(1 - y[i,j] for i in INDEX, j in INDEX where prev_sol[i,j] == 0) <= 
                          # N*N - 1)
    # end
    
    # # Main solution loop
    # solutions_found = 0
    # while solutions_found < num_sol
        # model, y = create_base_model()
        
        # # Add exclusion constraints for all previously found solutions
        # for prev_sol in Results
            # add_exclusion_constraint(model, y, prev_sol)
        # end
        
        # # Optimize
        # optimize!(model)
        
        # # Check if we found a solution
        # if termination_status(model) != MOI.OPTIMAL
            # if verbose && solutions_found == 0
                # println("No (more) solutions exist")
            # end
            # break
        # end
        
        # # Extract solution
        # sol = zeros(Int64, N, N)
        # for i in INDEX, j in INDEX
            # sol[i,j] = round(Int64, value(y[i,j]))
        # end
        
        # # Add to results
        # push!(Results, sol)
        # solutions_found += 1
        
        # if verbose
            # println("Found solution $solutions_found:")
            # display(sol)
            # println()
        # end
    # end
    
    # if solutions_found == num_sol
        # println("Warning! Maybe not all solutions were found. Try increasing num_sol.")
    # end
    
    # if verbose
        # println("Total solutions found: ", length(Results))
    # end
    
    # @show Results
    # return Results
# end

#david
# function encontrar_prototablones(rowsum::Vector{Int}, colsum::Vector{Int}; verbose::Bool = false, num_sol::Int = 1000)
    # ## Model:
    # model = Model(HiGHS.Optimizer)
    # @assert length(rowsum) == length(colsum)

    # N = length(rowsum)
    # INDEX = 1:N
    # @variable(model, y[INDEX, INDEX], Bin)  # Binary decision variables for faster computation

    # @constraint(model, rowCons[i in INDEX], sum(y[i, j] for j in INDEX) == rowsum[i])
    # @constraint(model, colCons[j in INDEX], sum(y[i, j] for i in INDEX) == colsum[j])

    # if verbose; println(model); end

    # ## Set verbosity for HiGHS solver
    # if !verbose
        # set_silent(model)
    # end

    # Results = Array{Matrix{Int64}, 1}()

    # ## Add an array to track excluded solutions
    # excluded_solutions = []

    # ## Optimization loop
    # for _ in 1:num_sol
        # optimize!(model)

        # if termination_status(model) != MOI.OPTIMAL
            # break
        # end

        # # Extract the solution
        # sol = Int.(JuMP.value.(y))

        # # Add to results
        # push!(Results, sol)

        # if verbose
            # println("Found solution:")
            # println(sol)
        # end

        # # Add exclusion constraint for the current solution
        # push!(excluded_solutions, sol)
        # @constraint(model, sum(y[i, j] != sol[i, j] for i in INDEX, j in INDEX) >= 1)
    # end

    # @show unique(Results)
    # return unique(Results)
# end

# function encontrar_prototablones(rowsum::Array{Int64,1}, colsum::Array{Int64,1}; verbose::Bool = false, num_sol::Int = 1000)
    # ## Model:
    # @show rowsum, colsum
    # model = Model(HiGHS.Optimizer)
    # @assert length(rowsum) == length(colsum)

    # N = length(rowsum)
    # INDEX = 1:N
    # @variable(model, y[INDEX,INDEX] >= 0, Int)

    # @constraint(model, rowCons[i=INDEX], sum(y[i,j] for j in INDEX) == rowsum[i])
    # @constraint(model, colCons[j=INDEX], sum(y[i,j] for i in INDEX) == colsum[j])

    # if verbose; println(model); end

    # ## Set verbosity for HiGHS solver
    # if !verbose
        # set_silent(model)
    # end

    # Results = Array{Array{Int64},1}()

    # ## Function to add a constraint to exclude a solution
    # function exclude_solution!(model, sol)
        # @constraint(model, sum((y[i,j] != sol[i,j]) for i in INDEX, j in INDEX) >= 1)
    # end

    # ## Find all solutions
    # while true
        # optimize!(model)

        # if termination_status(model) != MOI.OPTIMAL
            # break
        # end

        # sol = Array{Int64}(undef, N, N)
        # for i in INDEX
            # for j in INDEX
                # sol[i,j] = Int(value(y[i,j]))
            # end
        # end

        # push!(Results, sol)

        # if verbose
            # println("Found solution:")
            # println(sol)
        # end

        # exclude_solution!(model, sol)
    # end

    # @show Results
    # return Results
# end
#
# 
# function encontrar_prototablones(rowsum::Vector{Int}, colsum::Vector{Int}; verbose::Bool = false, num_sol::Int = 1000)
# function find_all_solutions(A::Vector{Int}, B::Vector{Int})
function encontrar_prototablones(A::Vector{Int}, B::Vector{Int})
    solutions = []
    n = length(A)
    
    # Create base model
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    # Variables
    @variable(model, M[1:n, 1:n] >= 0, Int)
    
    # Base constraints
    for i in 1:n
        @constraint(model, sum(M[i,:]) == A[i])
    end
    for j in 1:n
        @constraint(model, sum(M[:,j]) == B[j])
    end
    
    # Keep finding solutions until no more exist
    prev_size = -1
    current_size = 0
    
    while current_size > prev_size
        optimize!(model)
        
        if termination_status(model) != MOI.OPTIMAL
            break
        end
        
        # Get current solution
        current_sol = round.(Int, value.(M))
        push!(solutions, copy(current_sol))
        
        # Update sizes for termination check
        prev_size = current_size
        current_size = length(solutions)
        
        # Add constraint to exclude current solution
        constraint_expr = AffExpr(0.0)
        for i in 1:n, j in 1:n
            if current_sol[i,j] == 0
                # For zero elements, count how many become non-zero
                constraint_expr += M[i,j]
            else
                # For non-zero elements, count elements that decrease
                constraint_expr += (current_sol[i,j] - M[i,j])/(2*current_sol[i,j])
                # Count elements that increase
                constraint_expr += (M[i,j] - current_sol[i,j])/(2*current_sol[i,j])
            end
        end
        @constraint(model, constraint_expr >= 1)
    end
    
    println("Found $(length(solutions)) total solutions")
    return solutions
end


# function encontrar_prototablones(rowsum::Vector{Int}, colsum::Vector{Int}; verbose::Bool = false, num_sol::Int = 1000)
    # ## Model:
    # @show rowsum, colsum
    # model = Model(HiGHS.Optimizer)
    # @assert length(rowsum) == length(colsum)

    # N = length(rowsum)
    # INDEX = 1:N
    # @variable(model, y[INDEX, INDEX] >= 0, Int)

    # @constraint(model, rowCons[i=INDEX], sum(y[i, j] for j in INDEX) == rowsum[i])
    # @constraint(model, colCons[j=INDEX], sum(y[i, j] for i in INDEX) == colsum[j])

    # if verbose; println(model); end

    # ## Set verbosity for HiGHS solver
    # if !verbose
        # set_silent(model)
    # end

    # Results = []

    # ## Function to add a constraint to exclude a solution
    # function exclude_solution!(model, sol)
        # @constraint(model, sum(abs(y[i, j] - sol[i, j]) for i in INDEX, j in INDEX) >= 1)
    # end

    # ## Find all solutions
    # while length(Results) < num_sol
        # optimize!(model)

        # if termination_status(model) != MOI.OPTIMAL
            # break
        # end

        # sol = [Int(value(y[i, j])) for i in INDEX, j in INDEX]

        # push!(Results, sol)

        # if verbose
            # println("Found solution:")
            # println(sol)
        # end

        # exclude_solution!(model, sol)
    # end

    # @show length(Results), Results
    # return Results
# end

# function encontrar_prototablones(rowsum::Array{Int64,1}, colsum::Array{Int64,1}; verbose::Bool = false, num_sol::Int = 1000)
    # ## Model:
    # @show rowsum, colsum
    # model = Model(HiGHS.Optimizer)
    # @assert length(rowsum) == length(colsum)

    # N = length(rowsum)
    # INDEX = 1:N
    # @variable(model, y[INDEX,INDEX] >= 0, Int)

    # @constraint(model, rowCons[i=INDEX], sum(y[i,j] for j in INDEX) == rowsum[i])
    # @constraint(model, colCons[j=INDEX], sum(y[i,j] for i in INDEX) == colsum[j])

    # if verbose; println(model); end

    # ## Set verbosity for HiGHS solver
    # if !verbose
        # set_silent(model)
    # end

    # Results = Array{Array{Int64},1}()

    # ## Function to add a constraint to exclude a solution
    # function exclude_solution!(model, sol)
        # @constraint(model, sum(y[i,j] * sol[i,j] for i in INDEX, j in INDEX) <= sum(sol) - 1)
    # end

    # ## Find all solutions
    # while true
        # optimize!(model)

        # if termination_status(model) != MOI.OPTIMAL
            # break
        # end

        # sol = Array{Int64}(undef, N, N)
        # for i in INDEX
            # for j in INDEX
                # sol[i,j] = Int(value(y[i,j]))
            # end
        # end

        # push!(Results, sol)

        # if verbose
            # println("Found solution:")
            # println(sol)
        # end

        # exclude_solution!(model, sol)
    # end

    # @show Results
    # return Results
# end
#
# function encontrar_prototablones(rowsum::Array{Int64,1}, colsum::Array{Int64,1}; verbose::Bool = false, num_sol::Int = 1000)
    # ## Model:
    # model = Model(HiGHS.Optimizer)
    # @assert length(rowsum) == length(colsum)

    # N = length(rowsum)
    # INDEX = 1:N
    # @variable(model, y[INDEX,INDEX] >= 0, Int)

    # @constraint(model, rowCons[i=INDEX], sum(y[i,j] for j in INDEX) == rowsum[i])
    # @constraint(model, colCons[j=INDEX], sum(y[i,j] for i in INDEX) == colsum[j])

    # if verbose; println(model); end

    # ## Set verbosity for HiGHS solver
    # if !verbose
        # set_silent(model)
    # end

    # ## Optimize:
    # optimize!(model)

    # ## Results:
    # if termination_status(model) != MOI.OPTIMAL
        # error("Optimization was not successful.")
    # end

    # Results = Array{Array{Int64},1}()
    # sol = Array{Int64}(undef, N, N)
    # for i in INDEX
        # for j in INDEX
            # sol[i,j] = Int(value(y[i,j]))
        # end
    # end
    # push!(Results, sol)

    # if verbose
        # println("Solution:")
        # println(sol)
    # end
    # @show Results

    # return (x -> vcat(x...)).(transpose(Results))
# end

# function encontrar_prototablones(rowsum::Array{Int64,1}, colsum::Array{Int64,1}; verbose::Bool = false, num_sol::Int = 1000)
        # ## Model:
    # model = Model(Gurobi.Optimizer)
    # @assert length(rowsum) == length(colsum)

    # N = length(rowsum)
    # INDEX = 1:N
    # @variable(model, y[INDEX,INDEX] >= 0, Int)

    # @constraint(model, rowCons[i=INDEX], sum(y[i,j] for j in INDEX) == rowsum[i])
    # @constraint(model, colCons[j=INDEX], sum(y[i,j] for i in INDEX) == colsum[j])

    # if verbose; print(model); end

    # ## Gurobi parameters - see https://www.gurobi.com/documentation/9.0/refman/finding_multiple_solutions.html
    # JuMP.set_optimizer_attribute(model, "PoolSearchMode", 2)   # exhaustive search mode
    # JuMP.set_optimizer_attribute(model, "PoolSolutions", num_sol)  # num_sol is an arbitrary (large enough) whole number

    # if !verbose; JuMP.set_optimizer_attribute(model, "OutputFlag", 0); end

    # ## Optimize:
    # JuMP.optimize!(model)

    # ## Results:
    # num_results = result_count(model)
    # if verbose; println("Number of results: ", num_results, "\n"); end

    # Results = Array{Array{Int64},1}()  ## Note the conversion to Int64
    # for n in 1:num_results
        # sol = Array{Int64}(undef, N, N)
        # for i in INDEX
            # for j in INDEX
                # sol[i,j] = JuMP.value(y[i,j]; result=n)
            # end
        # end

        # push!(Results,sol)
        # if verbose
            # println(Results[n])
        # end
    # end

    # length(Results) == num_sol && println("Warning! Maybe I did not obtain all the solutions. Try increasing num_sol.")

    # return (x->vcat(x...)).(transpose(Results))
    # # return Results
# end

@doc Markdown.doc"""
  encajar(lista::Vector{T}, n::T) where T <: Integer
> Embebe `lista` en un intervalo de tamaño n

# Examples:
```
julia> l = [1,2,3]; encajar(l,5)
[1,2,3,4,5]
julia> l = [2,3,4]; encajar(l,5)
[1,2,3,4,5]
julia> l = [3,4,5]; encajar(l,5)
[1,2,3,4,5]
```
"""
function encajar(lista::Vector{T}, n::T) where T <: Integer
  valor = minimum(lista)
    if valor > 1
      lista = lista .- (valor - 1)
    end
    valor = maximum(lista)
    if valor < n
      lista = vcat([lista, collect(Int64, valor+1:n)]...)
    end
    lista
end

function encontrar_representativos(c_a::Content, c_b::Content)
    proto = encontrar_prototablones(c_a, c_b)
    lista_proto_perm = sortperm.(map(x -> calcula_proto_permutacion(x |> collect), proto))
    lista_proto_perm = map(x -> encajar(x,length(c_a)), lista_proto_perm)
    map(x -> Perm(collect(Int64,x) |> sortperm), lista_proto_perm)
end

function encontrar_representativos(t_a, t_b)
    c_a = content(t_a)
    c_b = content(t_b)
    proto = encontrar_prototablones(c_a, c_b)
    lista_proto_perm = map(x -> calcula_proto_permutacion(x |> collect), proto)
    map(x -> Perm(x |> collect), sortperm.(lista_proto_perm))
end
###############################################################################
#
#   Codigo para encontrar tablones
#
###############################################################################
@doc Markdown.doc"""
    monomio(f::MapST2SST,g::MapST2SST, per::Perm, n::Int)
> Calcula el monomio a partir de dos funciones y una permutacion.

# Examples:
```
julia> f = Dict{Int,Int}(1=>1, 2=>2, 3=>3);
julia> g = Dict{Int,Int}(1=>1, 2=>2, 3=>3);
julia> monomio(f,g,Perm([1,2,3]), 3)
```
"""
function monomio(f::MapST2SST, g::MapST2SST, per::Perm, n::Int64)
    fff = Int[]
    ggg = Int[]
    mon = 1
    for x in 1:n
        try
           push!(fff,f[per[x]])
       catch y
           if isa(y, KeyError)
                push!(fff,per[x])
           end
       end
    end
    
     for x in 1:n
        try
           push!(ggg,g[x])
       catch y
           if isa(y, KeyError)
                push!(ggg,x)
           end
       end
    end
    
    colleccion = zip(fff,ggg) |> collect
    
    for (a,b) in colleccion
        mon *= SymEngine.symbols("u_$(a)_$(b)")
    end
    
    mon
end

function monomio_n(f::MapST2SST, g::MapST2SST, per::Perm, n::Int64, mat::Array{Complex{Float64}, 2})
    fff = Int[]
    ggg = Int[]
    mon = 1
    for x in 1:n
        try
           push!(fff,f[per[x]])
       catch y
           if isa(y, KeyError)
                push!(fff,per[x])
           end
       end
    end
    
     for x in 1:n
        try
           push!(ggg,g[x])
       catch y
           if isa(y, KeyError)
                push!(ggg,x)
           end
       end
    end
    
    colleccion = zip(fff,ggg) |> collect
    
    for (a,b) in colleccion
      mon *= mat[a,b]
    end
    
    mon
end

@doc Markdown.doc"""
  double_coset(vector μ, vector ν)
> Return the double coset representatives

# Examples:
```
julia> c_a = [2,1,0]; c_b = [1,0,2];
julia> double_coset(c_a, b_b)
```
"""
function double_coset(μ::Content, ν::Content)
    between = encontrar_representativos(μ, ν)

    a,b = calcular_sα.([μ, ν])
    lista_coset = Array{Perm,1}[]
    for γ in between
        lista_permutaciones = Perm[]
        for x in a, z in b
            push!(lista_permutaciones,x*γ*z)
        end
        push!(lista_coset, unique!(lista_permutaciones))
    end
    
    
    between, lista_coset
end


@doc Markdown.doc"""
  group_function[lista tablones standard, tab semi U, tab semi V]
> Return the size of the vector which represents the partition.

<<<<<<< HEAD
# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
=======
# Example:
```julia
>>>>>>> e0f07f5... mejora documentacion group_function.
julia> t = YoungTableau([2,1]); fill!(t, [1,2,3]);
julia> group_function([2,1,0], t, t)
```
"""
function group_function(λ::Irrep, tab_u::YTableau, tab_v::YTableau; verbose = false) 
    tablones = StandardYoungTableaux(filter(x -> x > 0, λ))
    i = indice_tablon_semistandard(tab_u)
    j = indice_tablon_semistandard(tab_v)
    f = genera_funcion(tab_u,λ)
    g = genera_funcion(tab_v,λ)
    
    # probablemente se pueda sustituir con sum(λ)
    n = tab_u |> content |> length
    
    inversos = sqrt((1/ Θ(tab_u,λ) )*(1/ Θ(tab_v,λ) ))
    if verbose 
      @show inversos
    end

    (lista_gamas, lista_cosets) = double_coset(content(tab_u, λ), content(tab_v, λ))

    if verbose
      @show length(vcat(lista_cosets...) |> unique), length(lista_gamas)
    end
    
    pol::Basic = zero(SymEngine.symbols("x"))#0.0
    #total::Complex{Float64} = zero(Complex{Float64})
    total::Basic = zero(Basic)
    
    for ind in 1:length(lista_gamas)
        γ = lista_gamas[ind]
        cjto_σ = lista_cosets[ind]
        mon = monomio(f, g, inv(γ), n)
        total = zero(Basic)
        for σ in cjto_σ
            total += generar_matriz(tablones, σ,λ)[i,j]
        end
        pol += (total*mon)
    end
    
    pol*inversos
end

@doc Markdown.doc"""
    group_function(λ::Irrep, tu::GTPattern, tv::GTPattern)
> Return the _symbolic_ group function corresponding to irrep `λ` and GT patterns
> `tu` and `tv`.

# Example:
```julia
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2]);
julia> group_function([2,1,0], t, t)
```
"""
function group_function(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern; verbose = false) 
  tab_u = pat_u |> YoungTableau
  tab_v = pat_v |> YoungTableau
    tablones = StandardYoungTableaux(filter(x -> x > 0, λ))
    i = indice_tablon_semistandard(tab_u)
    j = indice_tablon_semistandard(tab_v)
    f = genera_funcion(tab_u,λ)
    g = genera_funcion(tab_v,λ)
    
    # probablemente se pueda sustituir con sum(λ)
    n = tab_u |> content |> length
    
    #inversos = sqrt((1/Θ(tab_u,λ))*(1/Θ(tab_v,λ)))
    inversos = sqrt((1/ Θ(tab_u,λ))*(1/ Θ(tab_v,λ) ))
    if verbose 
      @show inversos
    end

    (lista_gamas, lista_cosets) = double_coset(content(tab_u, λ), content(tab_v, λ))

    if verbose
      @show length(vcat(lista_cosets...) |> unique), length(lista_gamas)
    end
    
    pol::Basic = zero(SymEngine.symbols("x"))#0.0
    #total::Float64 = zero(Float64)
    total::Basic = zero(Basic)
    
    for ind in 1:length(lista_gamas)
        γ = lista_gamas[ind]
        cjto_σ = lista_cosets[ind]
        mon = monomio(f, g, inv(γ), n)
        total = zero(Basic)
        for σ in cjto_σ
            total += generar_matriz(tablones, σ,λ)[i,j]
        end
        pol += (total*mon)
    end
    
    pol*inversos
end

@doc Markdown.doc"""
    group_function(λ::Irrep, tu::GTPattern, tv::GTPattern, mat::Array{Complex{Float64},2})
> Return the _numeric_ group function, for an SU(n) member `mat`, corresponding to irrep `λ` and a pair of GT patterns
> `tu` and `tv`.

```julia
julia> using RandomMatrices
julia> mat = rand(Haar(2),3)
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2]);
julia> group_function([2,1,0], t, t, mat)
```
"""
function group_function(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern, mat::Array{Complex{Float64}, 2}; verbose = false) 
    tab_u = pat_u |> YoungTableau
    tab_v = pat_v |> YoungTableau
    tablones = StandardYoungTableaux(filter(x -> x > 0, λ))
    i = indice_tablon_semistandard(tab_u)
    j = indice_tablon_semistandard(tab_v)
    f = genera_funcion(tab_u,λ)
    g = genera_funcion(tab_v,λ)
    
    # probablemente se pueda sustituir con sum(λ)
    n = tab_u |> content |> length
    
    inversos = sqrt((1/Θn(tab_u,λ))*(1/Θn(tab_v,λ)))
    if verbose 
      @show inversos
    end

    (lista_gamas, lista_cosets) = double_coset(content(tab_u, λ), content(tab_v, λ))

    if verbose
      @show length(vcat(lista_cosets...) |> unique), length(lista_gamas)
    end
    
    pol::Complex{Float64} = zero(Complex{Float64})
    total::Float64 = zero(Float64)
    
    for ind in 1:length(lista_gamas)
        γ = lista_gamas[ind]
        cjto_σ = lista_cosets[ind]
        mon = monomio_n(f, g, inv(γ), n,mat)
        total = 0.0
        for σ in cjto_σ
            total += generar_matriz(tablones, σ,λ)[i,j]
        end
        pol += (total*mon)
    end
    
    pol*inversos
end

@doc Markdown.doc"""
    group_function(λ::Irrep, tu::GTPattern, tv::GTPattern, mat::Array{Complex{Float64},2})
> Return the _numeric_ group function, for an SU(n) member `mat`, corresponding to irrep `λ` and STYT
> `tu` and `tv`.

# Example:
```julia
julia> using RandomMatrices
julia> mat = rand(Haar(2),3)
julia> t = YoungTableau([2,1]); fill!(t, [1,2,3]);
julia> group_function([2,1,0], t, t, mat)
```
"""
function group_function(λ::Irrep, tab_u::YTableau, tab_v::YTableau, mat::Array{Complex{Float64}, 2}; verbose = false) 
    tablones = StandardYoungTableaux(filter(x -> x > 0, λ))
    i = indice_tablon_semistandard(tab_u)
    j = indice_tablon_semistandard(tab_v)
    f = genera_funcion(tab_u,λ)
    g = genera_funcion(tab_v,λ)
    
    # probablemente se pueda sustituir con sum(λ)
    n = tab_u |> content |> length
    
    inversos = sqrt((1/Θn(tab_u,λ))*(1/Θn(tab_v,λ)))
    if verbose 
      @show inversos
    end

    (lista_gamas, lista_cosets) = double_coset(content(tab_u, λ), content(tab_v, λ))

    if verbose
      @show length(vcat(lista_cosets...) |> unique), length(lista_gamas)
    end
    
    pol::Complex{Float64} = zero(Complex{Float64})
    total::Float64 = zero(Float64)
    
    for ind in 1:length(lista_gamas)
        γ = lista_gamas[ind]
        cjto_σ = lista_cosets[ind]
        mon = monomio_n(f, g, inv(γ), n,mat)
        total = 0.0
        for σ in cjto_σ
            total += generar_matriz(tablones, σ,λ)[i,j]
        end
        pol += (total*mon)
    end
    
    pol*inversos
end

macro mma_str(s)
  mma_to_julia(s)
end

@doc Markdown.doc"""
    mma_to_julia(s::String)
> Calcula el polinomio de MMA a uno de Julia usando SymEngine

# Example:

```
julia > mma_to_julia("x[1, 1])
u_1_1
```
"""
function mma_to_julia(edo::String)
    edo = replace(edo, ", " => "_")
    edo = replace(edo, "["  => "_")
    edo = replace(edo, "]"  => "")

    edo = replace(edo, r"x(_\d_\d)(\^\d{1})" => s"*(SymEngine.symbols(\"u\1\")\2)")
    edo = replace(edo, r"(x_\d_\d)"          => s"SymEngine.symbols(\"\g<1>\")")

    edo = replace(edo, " " => "")
    edo = replace(edo, "x" => "u")

    eval(Meta.parse(edo))
end

@doc Markdown.doc"""
> Computes _zweight_ of a GTPattern.
> This array, if applied to each state of the irrep, is commonly known as the _weight diagram_ of an SU(n) irrep.

  zweight(x::GTPattern)

# Examples:

```julia
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2])
julia> zweight(t)
>[0.5,0.5]

```
"""
function zweight(gt::GTPattern)
    l = zeros(Int, length(gt.filas) + 1)
    l[1] = 0
    l[2:end] = reverse!(sum.((gt.filas)))
    total::Array{Float64,1} = zeros(Float64, length(gt.filas) - 1)#Float64[]
    @simd for k in 2:length(gt.filas)
      @inbounds total[k - 1] = (l[k] - (1/2)*(l[k+1] + l[k-1]))
        #push!(total, (l[k] - (1/2)*(l[k+1] + l[k-1])))
    end
    total
end

@doc Markdown.doc"""
> Computes _pweight_ of a GTPattern.
> This array is related to the occupation number.

  pweight(x::GTPattern)

# Examples:

```julia
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2])
julia> pweight(t)
>[0,1,2]

```
"""
function pweight(gt::GTPattern)
    l::Array{Int64,1} = zeros(Int, length(gt.filas) + 1)
    #l = reverse!(sum.((gt.filas)))
    @simd for i in 1:length(gt.filas)
      @inbounds l[i] = sum(gt.filas[i])
    end
    total::Array{Int64,1} = zeros(Int, length(l) - 1)
    @simd for k in 1:length(total)
      @inbounds total[k] = l[k]  - l[k+1]
    end
    total
end
