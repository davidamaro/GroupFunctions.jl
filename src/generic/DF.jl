import LinearAlgebra: transpose, det, eigvals
export group_function
export zweight, pweight
export find_tablaeux_fillings

const Irrep = Array{T,1} where T <: Integer
const YTableau = AbstractAlgebra.Generic.YoungTableau{T} where T <: Integer
const Content = Array{T,1} where T <: Integer
const MapST2SST = Dict{T,T} where T <: Integer


include("FindTables.jl")
using .FindTables



@doc Markdown.doc"""
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
    proto_tableaux = find_tablaeux_fillings(c_a, c_b) #find tableaux fillings
    
    # Calculate and process permutations
    dimension = length(c_a)
    
    # Transform proto-tableaux into permutations
    representatives = map(proto_tableaux) do proto
        # Convert to permutation and sort
        proto_perm = proto |> collect |> expand_frequency_matrix |> sortperm
        # Adjust permutation to fit dimension
        adjusted_perm = adjust_permutation_list(proto_perm, dimension)
        # Create final permutation
        Perm(sortperm(collect(Int64, adjusted_perm)))
    end
    
    return representatives
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
    proto_tableaux = find_tablaeux_fillings(content_a, content_b)
    
    # Transform proto-tableaux into permutations
    representatives = map(proto_tableaux) do proto
        proto_array = collect(proto)
        permutation = expand_frequency_matrix(proto_array)
        Perm(collect(sortperm(permutation)))
    end
    
    return representatives
end

###############################################################################
#
#   Codigo para encontrar tablones
#
###############################################################################
@doc Markdown.doc"""

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

function pad_partition(irrep::Irrep, n::Int)
    length(irrep) <= n || error("Partition length $(length(irrep)) exceeds matrix size $n")
    padded = Vector{Int}(irrep)
    if length(padded) < n
        append!(padded, zeros(Int, n - length(padded)))
    end
    return padded
end

@doc Markdown.doc"""
    character_weyl(λ::Irrep, mat::Array{Complex{Float64},2}) -> ComplexF64
> Compute the U(n) character using the Weyl/Schur determinant formula, without
> using the characteristic polynomial.

# Notes
- `mat` must be square; the length of `λ` must not exceed `size(mat,1)`.
- The formula evaluates the Schur polynomial at the eigenvalues of `mat`:
  `det(z_i^(λ_j+n-j)) / det(z_i^(n-j))`.
"""
function character_weyl(λ::Irrep, mat::Array{Complex{Float64},2})
    n = size(mat, 1)
    size(mat, 2) == n || error("Matrix must be square; got $(size(mat))")
    padded = pad_partition(λ, n)
    eigenvalues = eigvals(mat)

    numer = Matrix{ComplexF64}(undef, n, n)
    denom = Matrix{ComplexF64}(undef, n, n)
    @inbounds for i in 1:n
        z = eigenvalues[i]
        for j in 1:n
            numer[i, j] = z^(padded[j] + n - j)
            denom[i, j] = z^(n - j)
        end
    end

    return det(numer) / det(denom)
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

function normalization_factor(tab_u::YTableau, tab_v::YTableau, irrep::Irrep)
    return sqrt((1 / Θ(tab_u, irrep)) * (1 / Θ(tab_v, irrep)))
end

function coset_debug_stats(coset_list)
    return length(unique(vcat(coset_list...)))
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
function compute_monomial(mapping1::MapST2SST, mapping2::MapST2SST, permutation::Perm, size::Int64, matrix::Matrix{ComplexF64})
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

@doc Markdown.doc"""
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


@doc Markdown.doc"""
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
    standard_tableaux = standard_tableaux_for_irrep(λ)
    index_u, index_v = semistandard_indices(tab_u, tab_v)
    mapping_u, mapping_v = semistandard_maps(tab_u, tab_v, λ)
    dimension = length(content(tab_u))
    normalization = normalization_factor(tab_u, tab_v, λ)

    if verbose
        @show normalization
    end

    content_u = content(tab_u, λ)
    content_v = content(tab_v, λ)
    gamma_list, coset_list = double_coset(content_u, content_v)

    if verbose
        unique_elements = coset_debug_stats(coset_list)
        total_gammas = length(gamma_list)
        @show unique_elements, total_gammas
    end

    polynomial = zero(Basic)
    @inbounds for (gamma, coset) in zip(gamma_list, coset_list)
        monomial_term = monomial(mapping_u, mapping_v, inv(gamma), dimension)
        coset_sum = zero(Basic)
        @inbounds for σ in coset
            coset_sum += generate_matrix(standard_tableaux, σ, λ)[index_u, index_v]
        end
        polynomial += coset_sum * monomial_term
    end

    return polynomial * normalization
end

@doc Markdown.doc"""
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
    # Convert GT patterns to Young tableaux
    tableau_u = YoungTableau(pat_u)
    tableau_v = YoungTableau(pat_v)
    
    # Initialize standard tableaux and indices
    standard_tableaux = StandardYoungTableaux(filter(x -> x > 0, λ))
    index_u = index_of_semistandard_tableau(tableau_u)
    index_v = index_of_semistandard_tableau(tableau_v)
    
    # Generate mapping functions
    mapping_u = standard_to_semistandard_map(tableau_u, λ)
    mapping_v = standard_to_semistandard_map(tableau_v, λ)
    
    # Calculate dimension from content
    dimension = length(content(tableau_u))
    
    # Calculate normalization factor
    normalization = sqrt((1 / Θ(tableau_u, λ)) * (1 / Θ(tableau_v, λ)))
    if verbose
        @show normalization
    end
    
    # Generate double cosets
    content_u = content(tableau_u, λ)
    content_v = content(tableau_v, λ)
    gamma_list, coset_list = double_coset(content_u, content_v)
    
    if verbose
        unique_elements = length(unique(vcat(coset_list...)))
        total_gammas = length(gamma_list)
        @show unique_elements, total_gammas
    end
    
    # Initialize polynomial and accumulator
    polynomial::Basic = zero(Basic)
    
    # Process each gamma and corresponding coset
    @inbounds for (index, gamma) in enumerate(gamma_list)
        coset = coset_list[index]
        monomial_term = monomial(mapping_u, mapping_v, inv(gamma), dimension)
        
        # Sum over current coset
        coset_sum = zero(Basic)
        @inbounds for permutation in coset
            coset_sum += generate_matrix(standard_tableaux, permutation, λ)[index_u, index_v]
        end
        
        polynomial += coset_sum * monomial_term
    end
    
    return polynomial * normalization
end

@doc Markdown.doc"""
    group_function(λ::Irrep, tu::GTPattern, tv::GTPattern, mat::Array{Complex{Float64},2})
> Returns the _numeric_ group function, for an SU(n) member `mat`, corresponding to irrep `λ` and a pair of GT patterns
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
    i = index_of_semistandard_tableau(tab_u)
    j = index_of_semistandard_tableau(tab_v)
    f = standard_to_semistandard_map(tab_u,λ)
    g = standard_to_semistandard_map(tab_v,λ)
    
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
        mon = compute_monomial(f, g, inv(γ), n,mat)
        total = 0.0
        for σ in cjto_σ
            total += generate_matrix(tablones, σ,λ)[i,j]
        end
        pol += (total*mon)
    end
    
    pol*inversos
end

@doc Markdown.doc"""
    group_function(λ::Irrep, tu::GTPattern, tv::GTPattern, mat::Array{Complex{Float64},2})
> Returns the _numeric_ group function, for an SU(n) member `mat`, corresponding to irrep `λ` and STYT
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
    i = index_of_semistandard_tableau(tab_u)
    j = index_of_semistandard_tableau(tab_v)
    f = standard_to_semistandard_map(tab_u,λ)
    g = standard_to_semistandard_map(tab_v,λ)
    
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
        mon = compute_monomial(f, g, inv(γ), n,mat)
        total = 0.0
        for σ in cjto_σ
            total += generate_matrix(tablones, σ,λ)[i,j]
        end
        pol += (total*mon)
    end
    
    pol*inversos
end

@doc Markdown.doc"""
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

@doc Markdown.doc"""
    group_function(λ::Irrep, mat::Array{Complex{Float64}, 2}; verbose::Bool = false) -> Tuple{Matrix{ComplexF64}, Vector{GTPattern}}

Compute all numeric group functions associated with the partition `λ` and a
matrix `mat`. Generates the GT patterns for `λ` and evaluates every pair using
the provided matrix.

Arguments:
- `λ::Irrep`: Partition describing the irrep
- `mat::Array{Complex{Float64}, 2}`: Matrix representing the SU(n) element
- `verbose::Bool`: Forwarded to the underlying pairwise `group_function`

Returns:
- `Tuple`: `(values, patterns)` where `values[i,j]` corresponds to
  `group_function(λ, patterns[i], patterns[j], mat)` and `patterns` is the basis
  returned by `basis_states(λ)`
"""
function group_function(λ::Irrep, mat::Array{Complex{Float64}, 2}; verbose::Bool = false)
    patterns = basis_states(λ)
    n_states = length(patterns)
    values = Matrix{ComplexF64}(undef, n_states, n_states)

    @inbounds for (i, pat_u) in enumerate(patterns)
        for (j, pat_v) in enumerate(patterns)
            values[i, j] = group_function(λ, pat_u, pat_v, mat; verbose = verbose)
        end
    end

    return values, patterns
end

# macro mma_str(s)
  # mma_to_julia(s)
# end


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
import Combinatorics: permutations
