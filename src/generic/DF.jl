import Combinatorics: permutations
using JuMP
using HiGHS
import LinearAlgebra: transpose
import StaticArrays: SArray
export group_function, mma_to_julia, @mma
export zweight, pweight
export encontrar_prototablones

const Irrep = Array{T,1} where T <: Integer
const YTableau = AbstractAlgebra.Generic.YoungTableau{T} where T <: Integer
const Content = Array{T,1} where T <: Integer
const MapST2SST = Dict{T,T} where T <: Integer

function flatten(x)
    n, _ = x |> size
    reshape(x, (n^2,1))
end

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
