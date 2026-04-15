export group_function, group_function_sym
export zweight, pweight
export find_tablaeux_fillings
export find_double_coset_representative_matrices
export occupation_number 

const Irrep = Array{T,1} where T <: Integer
const YTableau = YoungTableau
const Content = Array{T,1} where T <: Integer
const MapST2SST = Dict{T,T} where T <: Integer


include("AllSolutionsMatrix.jl")
using .AllSolutionsMatrix
include("FindTables.jl")
using .FindTables



@doc """
    adjust_permutation_list(lista::Vector{T}, n::T) where T <: Integer -> Vector{T}

Extend a permutation, represented by a vector if integers, to fit within range [1,n]
 by shifting values and filling gaps.

Arguments:
- `list::Vector{T}`: Input (possibly incomplete, without 1-cycles) permutation
- `n::T`: Target maximum value

Returns:
- `Vector{T}`: Adjusted permutation – a vector containing all integers from min to n

Notes:
- Shifts values down if minimum is greater than 1
- Fills missing values up to n if maximum is less than n
- Preserves input type T

# Examples:
```
julia> l = [1,2,3]; adjust_permutation_list(l,5)
[1,2,3,4,5]
julia> l = [2,3,4]; adjust_permutation_list(l,5)
[1,2,3,4,5]
julia> l = [3,4,5]; adjust_permutation_list(l,5)
[1,2,3,4,5]
julia> l = [3,2,1]; adjust_permutation_list(l,5)
[3,2,1,4,5]
```
"""
function adjust_permutation_list(list::Vector{T}, n::T) where T <: Integer
    #TODO: refactor the function – right now it doesn't perform checks if it is indeed a permutation
    # Return early if list is empty
    isempty(list) && return collect(T, 1:n)
    
    # Create a copy to avoid modifying input
    adjusted_list = copy(list)
    
    # Shift values down if needed
    min_value = minimum(adjusted_list)
    if min_value > 1
        adjusted_list .-= (min_value - 1)
    end
    
    # Extend list up to n if needed
    max_value = maximum(adjusted_list)
    if max_value < n
        append!(adjusted_list, (max_value + 1):n)
    end
    
    return adjusted_list
end

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
    find_double_coset_representatives(c_a::Content, c_b::Content) -> Vector{Perm}
was: encontrar_representativos
Find representative permutations between two content vectors.

Arguments:
- `c_a::Content`: First content vector
- `c_b::Content`: Second content vector

Returns:
- Vector{Perm}: List of representative permutations
"""
function find_double_coset_representatives(c_a::Content, c_b::Content)
    # Get proto-tableaux for the given contents
    proto_tableaux = enumerate_matrices(c_a, c_b)
    
    # Calculate and process permutations
    dimension = length(c_a)
    
    # Transform proto-tableaux into permutations
    representatives = map(proto_tableaux) do proto_rows
        # Convert to permutation and sort
        proto_matrix = rows_to_matrix(proto_rows)
        proto_perm = expand_frequency_matrix(proto_matrix) |> sortperm
        # Adjust permutation to fit dimension
        adjusted_perm = adjust_permutation_list(proto_perm, dimension)
        # Create final permutation
        Perm(sortperm(collect(Int64, adjusted_perm)))
    end
    
    return representatives
end

"""
    find_double_coset_representative_matrices(c_a::Content, c_b::Content) -> Vector{Matrix{Int}}

Return the frequency matrices (proto-tableaux) used to build the double coset
representative permutations.
"""
function find_double_coset_representative_matrices(c_a::Content, c_b::Content)
    return map(rows_to_matrix, enumerate_matrices(c_a, c_b))
end

"""
    find_double_coset_representatives(t_a, t_b) -> Vector{Perm}

Find representative permutations between two tableaux.

Arguments:
- `t_a::AbstractTableau`: First tableau
- `t_b::AbstractTableau`: Second tableau

Returns:
- Vector{Perm}: List of representative permutations
"""
function find_double_coset_representatives(t_a, t_b)
    # Extract content from tableaux
    content_a = content(t_a)
    content_b = content(t_b)
    #TODO: make this function just a wrapper using the implementation above (with c_a::Content, c_b::Content)

    # Get proto-tableaux for the contents
    proto_tableaux = enumerate_matrices(content_a, content_b)
    
    # Transform proto-tableaux into permutations
    representatives = map(proto_tableaux) do proto_rows
        proto_matrix = rows_to_matrix(proto_rows)
        permutation = expand_frequency_matrix(proto_matrix)
        Perm(collect(sortperm(permutation)))
    end
    
    return representatives
end

"""
    find_double_coset_representative_matrices(t_a, t_b) -> Vector{Matrix{Int}}

Return the frequency matrices (proto-tableaux) for two tableaux by using their
contents.
"""
function find_double_coset_representative_matrices(t_a, t_b)
    content_a = content(t_a)
    content_b = content(t_b)
    return map(rows_to_matrix, enumerate_matrices(content_a, content_b))
end

###############################################################################
#
#   Codigo para encontrar tablones
#
###############################################################################
@doc """

    monomial(f::MapST2SST, g::MapST2SST, per::Perm, n::Int64) -> Basic

    (was:monomio)

Compute a symbolic monomial by combining mapped indices into SymEngine variables.

Arguments:
- `f::MapST2SST`: First dict mapping indices to indices, {i => f(i)}
- `g::MapST2SST`: Second dict mapping indices to indices, {j => g(j)}
- `per::Perm`: Permutation vector
- `n::Int64`: Size of the permutation

Returns:
- `Basic`: SymEngine monomial, product of u_{f(per[k])}_{g(k)}} for k=1:n 


Notes:
- Uses default identity mapping when key is not found in f or g
- Creates symbolic variables in the form u_i_j where i,j are indices

# Examples:
```
julia> f = Dict{Int,Int}(1=>1, 2=>2, 3=>3);
julia> g = Dict{Int,Int}(1=>1, 2=>2, 3=>3);
julia> monomial(f,g,Perm([1,2,3]), 3)
```
"""
function monomial(f::MapST2SST, g::MapST2SST, per::Perm, n::Int64)
    
    # Pre-allocate arrays with known size
    mapped_indices1 = Vector{Int}(undef, n)
    mapped_indices2 = Vector{Int}(undef, n)
    
    # Process first mapping with permutation
    @inbounds for i in 1:n
        mapped_indices1[i] = get(f, per[i], per[i])
    end
    
    # Process second mapping
    @inbounds for i in 1:n
        mapped_indices2[i] = get(g, i, i)
    end
    
    # Initialize monomial product
    result = one(Basic)
    
    # Compute symbolic product
    @inbounds for i in 1:n
        result *= SymEngine.symbols("u_$(mapped_indices1[i])_$(mapped_indices2[i])")
    end
    
    return result
end

function standard_tableaux_for_irrep(irrep::Irrep)
    return StandardYoungTableaux(filter(x -> x > 0, irrep))
end

function semistandard_indices(tab_u::YTableau, tab_v::YTableau)
    return index_of_semistandard_tableau(tab_u), index_of_semistandard_tableau(tab_v)
end

function semistandard_maps(tab_u::YTableau, tab_v::YTableau, irrep::Irrep)
    return standard_to_semistandard_map(tab_u, irrep), standard_to_semistandard_map(tab_v, irrep)
end

function normalization_factor(tab_u::YTableau, tab_v::YTableau, irrep::Irrep,
                              ::Type{<:Union{AbstractFloat, Complex{<:AbstractFloat}}})
    return sqrt((1 / Θn(tab_u, irrep)) * (1 / Θn(tab_v, irrep)))
end

function normalization_factor(tab_u::YTableau, tab_v::YTableau, irrep::Irrep, ::Type=Basic)
    return sqrt((1 / Θ(tab_u, irrep)) * (1 / Θ(tab_v, irrep)))
end

function coset_debug_stats(coset_list)
    return length(unique(vcat(coset_list...)))
end

function symbolic_matrix(n::Int)
    return [SymEngine.symbols("u_$(i)_$(j)") for i in 1:n, j in 1:n]
end

struct DoubleCosetSumContext{TG,TC,TM,TT,TTabs,TIrr,P,T}
    gamma_list::TG
    coset_list::TC
    map_u::TM
    map_v::TM
    dim::Int
    mat::TT
    standard_tableaux::TTabs
    row_index::Int
    col_index::Int
    irrep::TIrr
    pol_zero::P
    total_zero::T
end

function sum_over_double_cosets(ctx::DoubleCosetSumContext{TG,TC,TM,TT,TTabs,TIrr,P,T}) where {TG,TC,TM,TT,TTabs,TIrr,P,T}
    pol::P = ctx.pol_zero
    @inbounds for (gamma, coset) in zip(ctx.gamma_list, ctx.coset_list)
        mon = compute_monomial(ctx.map_u, ctx.map_v, inv(gamma), ctx.dim, ctx.mat)
        total::T = ctx.total_zero
        @inbounds for sigma in coset
            total += convert(T, generate_matrix(ctx.standard_tableaux, sigma, ctx.irrep)[ctx.row_index, ctx.col_index])
        end
        pol += total * mon
    end
    return pol
end

"""
    compute_monomial(mapping1::MapST2SST, mapping2::MapST2SST, permutation::Perm, size::Int64, matrix::Matrix{ComplexF64}) -> ComplexF64

Compute a monomial term by combining elements from a complex matrix based on mappings and permutations.

Arguments:
- `mapping1::MapST2SST`: First mapping function that transforms indices through permutation
- `mapping2::MapST2SST`: Second mapping function that transforms original indices
- `permutation::Perm`: Vector representing the permutation of indices
- `size::Int64`: Size of the permutation/dimension of the problem
- `matrix::Matrix{ComplexF64}`: Complex matrix for element-wise calculations

Returns:
- `ComplexF64`: The computed monomial value by multiplying selected matrix elements

Notes:
- Uses identity mapping (original index) when a key is not found in mapping1 or mapping2
- Matrix dimensions must be compatible with the mappings and size parameter
"""
function compute_monomial(mapping1::MapST2SST, mapping2::MapST2SST, permutation::Perm, size::Int64, matrix::AbstractMatrix{T}) where T
    # Pre-allocate arrays for transformed indices
    transformed_indices1 = Vector{Int}(undef, size)
    transformed_indices2 = Vector{Int}(undef, size)

    # Apply first mapping with permutation
    @inbounds for index in 1:size
        transformed_indices1[index] = get(mapping1, permutation[index], permutation[index])
    end

    # Apply second mapping
    @inbounds for index in 1:size
        transformed_indices2[index] = get(mapping2, index, index)
    end

    # Calculate the monomial by multiplying corresponding matrix elements
    @inbounds monomial = prod(matrix[transformed_indices1[i], transformed_indices2[i]] for i in 1:size)

    return monomial
end

@doc """
  double_coset(vector μ, vector ν)
> Return the double coset representatives

# Examples:
```
julia> c_a = [2,1,0]; c_b = [1,0,2];
julia> double_coset(c_a, b_b)
```
    double_coset(μ::Content, ν::Content) -> Tuple{Vector{Perm}, Vector{Vector{Perm}}}

Calculate double coset representatives and their corresponding cosets for given contents.

Arguments:
- `μ::Content`: First content vector
- `ν::Content`: Second content vector

Returns:
- Tuple containing:
  - Vector{Perm}: Representatives between the contents
  - Vector{Vector{Perm}}: List of unique cosets for each representative

Notes:
- Uses stabilizer_permutations for stabilizer computations
- Performs unique factorization of permutations
"""
function double_coset(μ::Content, ν::Content)
    # Find representatives between contents
    representatives = find_double_coset_representatives(μ, ν)
    
    # Calculate stabilizers for both contents
    stabilizer_left, stabilizer_right = stabilizer_permutations.([μ, ν])

    # Precompute cartesian size to sizehint! and reuse a set for dedup
    reps_count = length(representatives)
    coset_list = Vector{Vector{Perm}}(undef, reps_count)

    @inbounds for (idx, representative) in enumerate(representatives)
        coset = Vector{Perm}()
        sizehint!(coset, length(stabilizer_left) * length(stabilizer_right))
        for left_perm in stabilizer_left
            for right_perm in stabilizer_right
                push!(coset, left_perm * representative * right_perm)
            end
        end
        unique!(coset)
        coset_list[idx] = coset
    end
    
    return representatives, coset_list
end


@doc """
    group_function(λ::Irrep, tu::YoungTableau, tv::YoungTableau)
> Returns the _symbolic_ group function corresponding to irrep `λ` and Young tableaux
> `tu` and `tv`.

# Example:
```julia
julia> t = YoungTableau([2,1]); fill!(t, [1,2,3]);
julia> group_function([2,1,0], t, t)
```

    group_function(λ::Irrep, tab_u::YTableau, tab_v::YTableau; verbose::Bool = false) -> Basic

Compute group theoretical function based on Young tableaux and irreducible representations.

Arguments:
- `λ::Irrep`: Irreducible representation
- `tab_u::YTableau`: First Young tableau
- `tab_v::YTableau`: Second Young tableau
- `verbose::Bool`: Flag for detailed output (default: false)

Returns:
- `Complex`: Group function evaluated

Notes:
- Uses SymEngine for symbolic computation
- Involves matrix operations and coset calculations
"""
function group_function(λ::Irrep, tab_u::YTableau, tab_v::YTableau; verbose::Bool = false)
    n = length(λ)  # group dimension SU(n), not number of boxes
    mat = symbolic_matrix(n)
    return group_function(λ, tab_u, tab_v, mat; verbose = verbose)
end

@doc """
    group_function(λ::Irrep, tu::GTPattern, tv::GTPattern)
> Returns the _symbolic_ group function corresponding to irrep `λ` and GT patterns
> `tu` and `tv`.

# Example:
```julia
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2]);
julia> group_function([2,1,0], t, t)
```
    group_function(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern; verbose::Bool = false) -> Basic

Compute group theoretical function based on Gelfand-Tsetlin patterns and irreducible representations.

Arguments:
- `λ::Irrep`: Irreducible representation
- `pat_u::GTPattern`: First Gelfand-Tsetlin pattern
- `pat_v::GTPattern`: Second Gelfand-Tsetlin pattern
- `verbose::Bool`: Flag for detailed output (default: false)

Returns:
- `Basic`: Computed polynomial expression in SymEngine format

Notes:
- Converts GT patterns to Young tableaux for calculations
- Uses SymEngine for symbolic computation
"""
function group_function(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern; verbose::Bool = false)
    return group_function(λ, YoungTableau(pat_u), YoungTableau(pat_v); verbose = verbose)
end

@doc """
    group_function(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern, mat::AbstractMatrix{T}; verbose=false) where T

Returns the group function for an SU(n) element `mat`, corresponding to irrep `λ` and a pair of GT patterns
`pat_u` and `pat_v`. Converts GT patterns to Young tableaux and delegates to the main implementation.

```julia
julia> using RandomMatrices
julia> mat = rand(Haar(2),3)
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2]);
julia> group_function([2,1,0], t, t, mat)
```
"""
function group_function(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern, mat::AbstractMatrix{T}; verbose = false) where T
    return group_function(λ, YoungTableau(pat_u), YoungTableau(pat_v), mat; verbose = verbose)
end

@doc """
    group_function_sym(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern, mat::AbstractMatrix{T}; verbose::Bool = false) where T

Compatibility wrapper for the symbolic API. It forwards to the generic
`group_function` implementation so symbolic matrices continue to work through
the legacy `group_function_sym` entry point.
"""
function group_function_sym(λ::Irrep, pat_u::GTPattern, pat_v::GTPattern, mat::AbstractMatrix{T}; verbose::Bool = false) where T
    return group_function(λ, pat_u, pat_v, mat; verbose = verbose)
end

@doc """
    group_function(λ::Irrep, tab_u::YTableau, tab_v::YTableau, mat::AbstractMatrix{T}; verbose=false) where T

Returns the group function for an SU(n) element `mat`, corresponding to irrep `λ` and
Young tableaux `tab_u` and `tab_v`. For numeric `mat` (real or complex floating point),
computes the group function numerically. Otherwise, uses exact coefficients internally,
and the result is a polynomial in the entries of `mat`.

# Example:
```julia
julia> using RandomMatrices
julia> mat = rand(Haar(2),3)
julia> t = YoungTableau([2,1]); fill!(t, [1,2,3]);
julia> group_function([2,1,0], t, t, mat)
```
"""
function group_function(λ::Irrep, tab_u::YTableau, tab_v::YTableau, mat::AbstractMatrix{T}; verbose = false) where T
    standard_tableaux = standard_tableaux_for_irrep(λ)
    row_index = index_of_semistandard_tableau(tab_u)
    col_index = index_of_semistandard_tableau(tab_v)
    map_u = standard_to_semistandard_map(tab_u, λ)
    map_v = standard_to_semistandard_map(tab_v, λ)

    # probablemente se pueda sustituir con sum(λ)
    dim = tab_u |> content |> length
    normalization = normalization_factor(tab_u, tab_v, λ, T)
    if verbose
      @show normalization
    end

    (gamma_list, coset_list) = double_coset(content(tab_u, λ), content(tab_v, λ))

    if verbose
      @show length(vcat(coset_list...) |> unique), length(gamma_list)
    end

    ctx = DoubleCosetSumContext(
        gamma_list,
        coset_list,
        map_u,
        map_v,
        dim,
        mat,
        standard_tableaux,
        row_index,
        col_index,
        λ,
        zero(eltype(mat)),
        zero(eltype(mat)),
    )
    pol = sum_over_double_cosets(ctx)

    pol * normalization
end

@doc """
    group_function(λ::Irrep; verbose::Bool = false) -> Tuple{Matrix{Basic}, Vector{GTPattern}}

Compute all symbolic group functions associated with the partition `λ`.
The routine builds every valid GT pattern for `λ` and evaluates the group
function for each pair.

Arguments:
- `λ::Irrep`: Partition describing the irrep
- `verbose::Bool`: Forwarded to the underlying pairwise `group_function`

Returns:
- `Tuple`: `(values, patterns)` where `values[i,j]` corresponds to
  `group_function(λ, patterns[i], patterns[j])` and `patterns` is the basis
  returned by `basis_states(λ)`
"""
function group_function(λ::Irrep; verbose::Bool = false)
    patterns = basis_states(λ)
    n_states = length(patterns)
    values = Matrix{Basic}(undef, n_states, n_states)

    @inbounds for (i, pat_u) in enumerate(patterns)
        for (j, pat_v) in enumerate(patterns)
            values[i, j] = group_function(λ, pat_u, pat_v; verbose = verbose)
        end
    end

    return values, patterns
end

@doc """
    group_function(λ::Irrep, mat::AbstractMatrix{T}; verbose::Bool = false) where T

Compute all group functions associated with the partition `λ` and a matrix `mat`.
Generates the GT patterns for `λ` and evaluates every pair using the provided matrix.
For numeric matrices, returns numeric values; for symbolic matrices, returns polynomials.

Arguments:
- `λ::Irrep`: Partition describing the irrep
- `mat::AbstractMatrix{T}`: Matrix representing the SU(n) element
- `verbose::Bool`: Forwarded to the underlying pairwise `group_function`

Returns:
- `Tuple`: `(values, patterns)` where `values[i,j]` corresponds to
  `group_function(λ, patterns[i], patterns[j], mat)` and `patterns` is the basis
  returned by `basis_states(λ)`
"""
function group_function(λ::Irrep, mat::AbstractMatrix{T}; verbose::Bool = false) where T
    patterns = basis_states(λ)
    n_states = length(patterns)
    first_value = group_function(λ, patterns[1], patterns[1], mat; verbose = verbose)
    values = Matrix{typeof(first_value)}(undef, n_states, n_states)

    values[1, 1] = first_value
    @inbounds for (i, pat_u) in enumerate(patterns)
        for (j, pat_v) in enumerate(patterns)
            if i == 1 && j == 1
                continue
            end
            values[i, j] = group_function(λ, pat_u, pat_v, mat; verbose = verbose)
        end
    end

    return values, patterns
end

# macro mma_str(s)
  # mma_to_julia(s)
# end


@doc """
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
    l = zeros(Int, length(gt.rows) + 1)
    l[1] = 0
    l[2:end] = reverse!(sum.((gt.rows)))
    total::Array{Float64,1} = zeros(Float64, length(gt.rows) - 1)#Float64[]
    @simd for k in 2:length(gt.rows)
      @inbounds total[k - 1] = (l[k] - (1/2)*(l[k+1] + l[k-1]))
        #push!(total, (l[k] - (1/2)*(l[k+1] + l[k-1])))
    end
    total
end

@doc """
    pweight(x::GTPattern)

Compute the `pweight` of a `GTPattern`.

This array is related to the occupation number.

# Examples
```jldoctest
julia> using GroupFunctions

julia> t = GTPattern([[2, 1, 0], [2, 1], [2]]);

julia> pweight(t)
3-element Vector{Int64}:
 0
 1
 2
```
"""
function pweight(gt::GTPattern)
    l::Array{Int64,1} = zeros(Int, length(gt.rows) + 1)
    #l = reverse!(sum.((gt.filas)))
    @simd for i in 1:length(gt.rows)
      @inbounds l[i] = sum(gt.rows[i])
    end
    total::Array{Int64,1} = zeros(Int, length(l) - 1)
    @simd for k in 1:length(total)
      @inbounds total[k] = l[k]  - l[k+1]
    end
    total
end


occupation_number = reverse ∘ pweight
import Combinatorics: permutations
