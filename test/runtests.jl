using GroupFunctions, Test, SymEngine, Combinatorics

# , Immanants
import RandomMatrices: Haar
import LinearAlgebra: norm, tr, det, eigvals

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
    filter(x -> norm(zweight(x)) ≈ 0.0, base)
end

function es_cero(pol::Basic; ϵ = 10^(-5))
    monomios_lista = SymEngine.free_symbols(pol)

    while length(monomios_lista) > 0
      mono = pop!(monomios_lista)
      pol = SymEngine.subs(pol, mono, randn())
    end

    abs(pol)^2 < ϵ
end

include("internals.jl")
include("double_coset_representatives.jl")
include("Aqua.jl")
include("mathematica_irrep_300.jl")
include("mathematica_additional_su3_irreps.jl")
include("mathematica_even_more_su3_irreps.jl")

@testset "Number of zero weight states." begin
  solutions = find_tableaux_fillings([1, 1, 2, 1, 0], [1, 1, 2, 1, 0])
  @test length(solutions) == 33
  solutions = find_tableaux_fillings([2,2,2,0,0,0], [2,2,2,0,0,0])
  @test length(solutions) == 21
end

@testset "group_function with partition input" begin
    λ = [1, 0]
    symbolic_values, patterns = group_function(λ)
    expected_patterns = basis_states(λ)

    @test map(p -> p.rows, patterns) == map(p -> p.rows, expected_patterns)
    @test size(symbolic_values) == (length(patterns), length(patterns))
    @test symbolic_values[1, 2] == group_function(λ, expected_patterns[1], expected_patterns[2])

    mat = ComplexF64[1 0; 0 1]
    numeric_values, numeric_patterns = group_function(λ, mat)

    @test map(p -> p.rows, numeric_patterns) == map(p -> p.rows, patterns)
    @test size(numeric_values) == (length(patterns), length(patterns))
    @test numeric_values[2, 1] ≈ group_function(λ, expected_patterns[2], expected_patterns[1], mat)
end

@testset "group_function with symbolic matrix input" begin
    λ = [2, 0]
    mat = su2_block_symbolic(2, 1)

    symbolic_values, patterns = group_function(λ, mat)
    expected_patterns = basis_states(λ)

    @test map(p -> p.rows, patterns) == map(p -> p.rows, expected_patterns)
    @test size(symbolic_values) == (length(patterns), length(patterns))
    @test eltype(symbolic_values) == Basic
    @test symbolic_values[1, 2] == group_function(λ, expected_patterns[1], expected_patterns[2], mat)
    @test GroupFunctions.character(λ, mat) == sum(symbolic_values[i, i] for i in axes(symbolic_values, 1))
end

@testset "character helper" begin
    λ = [2, 0]
    symbolic_character = GroupFunctions.character(λ)
    symbolic_rep, symbolic_patterns = group_function(λ)
    @test symbolic_character == sum(symbolic_rep[i, i] for i in axes(symbolic_rep, 1))
    @test symbolic_character isa Basic

    numeric_mat = ComplexF64[1 0; 0 1]
    numeric_rep, _ = group_function(λ, numeric_mat)
    @test @inferred(GroupFunctions.character(λ, numeric_mat)) ≈ tr(numeric_rep)
end

@testset "standard_to_semistandard_map is non-decreasing" begin
    t = YoungTableau([3]); fill!(t, [1,2,5])
    std = GroupFunctions.standard_tableau_from_semistandard(t)
    @test std.fill == [1,2,3]
    f = GroupFunctions.standard_to_semistandard_map(t, [3,0,0,0,0])
    fvals = [f[i] for i in 1:5]
    @test fvals == [1,2,5,5,5]
    @test issorted(fvals)

    t2 = YoungTableau([2,1]); fill!(t2, [1,1,2])
    f2 = GroupFunctions.standard_to_semistandard_map(t2)
    f2vals = [f2[i] for i in 1:3]
    @test issorted(f2vals)

    f2_ext = GroupFunctions.standard_to_semistandard_map(t2, [2,1,0,0])
    f2ext_vals = [f2_ext[i] for i in 1:4]
    @test f2ext_vals == [1,1,2,2]
    @test issorted(f2ext_vals)

    t3 = YoungTableau([2,1]); fill!(t3, [1,3,3])
    f3 = GroupFunctions.standard_to_semistandard_map(t3, [2,1,0,0,0])
    f3vals = [f3[i] for i in 1:5]
    @test issorted(f3vals)
    @test f3vals[4] == f3vals[3] == f3vals[5]
end

@testset "comparison of the character (2,1,0) of SU3" begin
    mat = rand(Haar(2), 3)
    nmat = mat / det(mat)^(1/3)
    lambdas = eigvals(nmat)
    ideal_character = lambdas[1]/lambdas[2] + lambdas[2]/lambdas[1] +
                      lambdas[1]/lambdas[3] + lambdas[3]/lambdas[1] +
                      lambdas[2]/lambdas[3] + lambdas[3]/lambdas[2] + 2
    @test @inferred(GroupFunctions.character([2,1,0], nmat)) ≈ ideal_character
end

@testset "comparison of the character (2,0)" begin
    mat = rand(Haar(2), 3)
    nmat = mat / det(mat)^(1/3)
    lambdas = eigvals(nmat)
    ideal_character = lambdas[1]^2 + lambdas[2]^2 + lambdas[3]^2 +
                      lambdas[1]*lambdas[2] + lambdas[1]*lambdas[3] +
                      lambdas[2]*lambdas[3]
    @test @inferred(GroupFunctions.character([2,0,0], nmat)) ≈ ideal_character
end

@testset "comparison of the character (3,0)" begin
    mat = rand(Haar(2), 3)
    nmat = mat / det(mat)^(1/3)
    lambdas = eigvals(nmat)
    ideal_character = lambdas[1]^3 + lambdas[2]^3 + lambdas[3]^3 +
                      lambdas[1]^2*lambdas[2] + lambdas[1]^2*lambdas[3] +
                      lambdas[2]^2*lambdas[1] + lambdas[2]^2*lambdas[3] +
                      lambdas[3]^2*lambdas[1] + lambdas[3]^2*lambdas[2] +
                      lambdas[1]*lambdas[2]*lambdas[3]
    @test @inferred(GroupFunctions.character([3,0,0], nmat)) ≈ ideal_character
end

@testset "irrep 221 de SU(4)" begin
    t_u = YoungTableau([2,2, 1])
    fill!(t_u, [1,2,3,3,4])
    t_v = YoungTableau([2,2, 1])
    fill!(t_v, [1,2,3,3,4])
    edo_mma = "0.5 x[1, 4] x[2, 3] x[3, 2] x[3, 3] x[4, 1] +  0.5 x[1, 3] x[2, 4] x[3, 2] x[3, 3] x[4, 1] -  0.5 x[1, 4] x[2, 2] x[3, 3]^2* x[4, 1] -  0.5 x[1, 2] x[2, 4] x[3, 3]^2* x[4, 1] -  1.0 x[1, 3] x[2, 3] x[3, 2] x[3, 4] x[4, 1] +  0.5 x[1, 3] x[2, 2] x[3, 3] x[3, 4] x[4, 1] +  0.5 x[1, 2] x[2, 3] x[3, 3] x[3, 4] x[4, 1] +  0.5 x[1, 4] x[2, 3] x[3, 1] x[3, 3] x[4, 2] +  0.5 x[1, 3] x[2, 4] x[3, 1] x[3, 3] x[4, 2] -  0.5 x[1, 4] x[2, 1] x[3, 3]^2* x[4, 2] -  0.5 x[1, 1] x[2, 4] x[3, 3]^2* x[4, 2] -  1.0 x[1, 3] x[2, 3] x[3, 1] x[3, 4] x[4, 2] +  0.5 x[1, 3] x[2, 1] x[3, 3] x[3, 4] x[4, 2] +  0.5 x[1, 1] x[2, 3] x[3, 3] x[3, 4] x[4, 2] -  1.0 x[1, 4] x[2, 3] x[3, 1] x[3, 2] x[4, 3] -  1.0 x[1, 3] x[2, 4] x[3, 1] x[3, 2] x[4, 3] +  0.5 x[1, 4] x[2, 2] x[3, 1] x[3, 3] x[4, 3] + 0.5 x[1, 2] x[2, 4] x[3, 1] x[3, 3] x[4, 3] +  0.5 x[1, 4] x[2, 1] x[3, 2] x[3, 3] x[4, 3] +  0.5 x[1, 1] x[2, 4] x[3, 2] x[3, 3] x[4, 3] +  0.5 x[1, 3] x[2, 2] x[3, 1] x[3, 4] x[4, 3] +  0.5 x[1, 2] x[2, 3] x[3, 1] x[3, 4] x[4, 3] +  0.5 x[1, 3] x[2, 1] x[3, 2] x[3, 4] x[4, 3] +  0.5 x[1, 1] x[2, 3] x[3, 2] x[3, 4] x[4, 3] -  1.0 x[1, 2] x[2, 1] x[3, 3] x[3, 4] x[4, 3] -  1.0 x[1, 1] x[2, 2] x[3, 3] x[3, 4] x[4, 3] +  2.0 x[1, 3] x[2, 3] x[3, 1] x[3, 2] x[4, 4] -  1.0 x[1, 3] x[2, 2] x[3, 1] x[3, 3] x[4, 4] -  1.0 x[1, 2] x[2, 3] x[3, 1] x[3, 3] x[4, 4] -  1.0 x[1, 3] x[2, 1] x[3, 2] x[3, 3] x[4, 4] -  1.0 x[1, 1] x[2, 3] x[3, 2] x[3, 3] x[4, 4] +  x[1, 2] x[2, 1] x[3, 3]^2* x[4, 4] + x[1, 1] x[2, 2] x[3, 3]^2* x[4, 4]" |> mma_to_julia
    edo_julia = group_function([2,2,1,0],t_u, t_v) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
end

@testset "SU(2) irrep homomorphism" begin
    u = rand(Haar(2), 2)
    u = u / det(u)^(1 / 2)
    v = rand(Haar(2), 2)
    v = v / det(v)^(1 / 2)
    uv = u * v

    for λ in ([2, 0], [3, 0])
        rep_u, _ = group_function(λ, u)
        rep_v, _ = group_function(λ, v)
        rep_uv, _ = group_function(λ, uv)
        @test isapprox(rep_u * rep_v, rep_uv; atol=1e-8, rtol=1e-8)
    end
end

@testset "SU(3) irrep homomorphism" begin
    u = rand(Haar(2), 3)
    u = u / det(u)^(1 / 3)
    v = rand(Haar(2), 3)
    v = v / det(v)^(1 / 3)
    uv = u * v

    for λ in ([2, 1, 0], [2, 2, 0])
        rep_u, _ = group_function(λ, u)
        rep_v, _ = group_function(λ, v)
        rep_uv, _ = group_function(λ, uv)
        @test isapprox(rep_u * rep_v, rep_uv; atol=1e-8, rtol=1e-8)
    end
end

@testset "Comparison between Legendre polynomails and D-function" begin
  α1,β1,γ1 = rand(Float64,3)
  xx=su2_block(3,1,(α1,β1,γ1))
  α2,β2 = rand(Float64,2)
  yy=su2_block(3,2,(α2,β2,α2))
  α3,β3,γ3 = rand(Float64,3)
  zz=su2_block(3,1,(α3,β3,γ3))

  uu = xx*yy*zz

  edo = GTPattern([[2,1,0],[1,1], [1]], [1])
  @test group_function([2,1,0],edo, edo, uu) ≈ (1/4)*(1+3*cos(β2))

  edo = GTPattern([[4,2,0],[2,2], [2]], [2])
  @test group_function([4,2,0],edo, edo, uu) ≈ (1/3^2)*(1+3*cos(β2)+5*(1/2)*(3*cos(β2)^2 - 1))
end

@testset "irrep 11 de SU(3)" begin
    t_1 = YoungTableau([1, 1])# edo 21
    fill!(t_1, [2,3])
    t_2 = YoungTableau([1, 1])# edo 69
    fill!(t_2, [1,3])
    t_3 = YoungTableau([1, 1])# edo 69
    fill!(t_3, [1,2])
  ##ya
    edo_mma = "-1.0 x[2, 3] x[3, 2] + x[2, 2] x[3, 3]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_1, t_1) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##ya
    edo_mma = "-1.0 x[2, 3] x[3, 1] + x[2, 1] x[3, 3]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_1, t_2) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##ya
    edo_mma = "-1.0 x[2, 2] x[3, 1] + x[2, 1] x[3, 2]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_1, t_3) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##ya
    edo_mma = "-1.0 x[1, 3] x[3, 2] + x[1, 2] x[3, 3]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_2, t_1) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##ya
    edo_mma = "-1.0 x[1, 3] x[3, 1] + x[1, 1] x[3, 3]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_2, t_2) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##ya
    edo_mma = "-1.0 x[1, 2] x[3, 1] + x[1, 1] x[3, 2]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_2, t_3) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##ya
    edo_mma = "-1.0 x[1, 3] x[2, 2] + x[1, 2] x[2, 3]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_3, t_1) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##
    edo_mma = "-1.0 x[1, 3] x[2, 1] + x[1, 1] x[2, 3]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_3, t_2) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
  ##
    edo_mma = "-1.0 x[1, 2] x[2, 1] + x[1, 1] x[2, 2]" |> mma_to_julia
    edo_julia = group_function([1,1,0],t_3, t_3) |> expand
    edo_final = expand(edo_mma - edo_julia)
    @test es_cero(edo_final)
end
@testset "irrep 21 de SU(3)" begin
    t_1 = YoungTableau([2,1])# edo 21
    fill!(t_1, [2,3,3])
    t_2 = YoungTableau([2,1])# edo 69
    fill!(t_2, [1,3,3])
    t_3 = YoungTableau([2,1])# edo 69
    fill!(t_3, [1,3,2])
    t_4 = YoungTableau([2,1])# edo 69
    fill!(t_4, [2,2,3])
    t_5 = YoungTableau([2,1])# edo 21
    fill!(t_5, [1,2,3])
    t_6 = YoungTableau([2,1])# edo 69
    fill!(t_6, [1,1,3])
    t_7 = YoungTableau([2,1])# edo 69
    fill!(t_7, [1,2,2])
    t_8 = YoungTableau([2,1])# edo 69
    fill!(t_8, [1,1,2])

      mma_states = [
      [
      "-1.0 x[2, 3] x[3, 2] x[3, 3] + x[2, 2] x[3, 3]^2",
      "-1.0 x[2, 3] x[3, 1] x[3, 3] + x[2, 1] x[3, 3]^2",
      "-1.224744871391589 x[2, 2] x[3, 1] x[3, 3] + 1.224744871391589 x[2, 1] x[3, 2] x[3, 3]",
      "-1.0 x[2, 3] x[3, 2]^2 + x[2, 2] x[3, 2] x[3, 3]",
      "-1.4142135623730951 x[2, 3] x[3, 1] x[3, 2] + 0.7071067811865475 x[2, 2] x[3, 1] x[3, 3] + 0.7071067811865475 x[2, 1] x[3, 2] x[3, 3]",
      "-1.0 x[2, 3] x[3, 1]^2 + x[2, 1] x[3, 1] x[3, 3]",
      "-1.0 x[2, 2] x[3, 1] x[3, 2] + x[2, 1] x[3, 2]^2",
      "-1.0 x[2, 2] x[3, 1]^2 + x[2, 1] x[3, 1] x[3, 2]"
      ],
      [
      "-1.0 x[1, 3] x[3, 2] x[3, 3] + x[1, 2] x[3, 3]^2",
      "-1.0 x[1, 3] x[3, 1] x[3, 3] + x[1, 1] x[3, 3]^2",
      "-1.224744871391589 x[1, 2] x[3, 1] x[3, 3] + 1.224744871391589 x[1, 1] x[3, 2] x[3, 3]",
      "-1.0 x[1, 3] x[3, 2]^2 + x[1, 2] x[3, 2] x[3, 3]",
      "-1.4142135623730951 x[1, 3] x[3, 1] x[3, 2] + 0.7071067811865475 x[1, 2] x[3, 1] x[3, 3] + 0.7071067811865475 x[1, 1] x[3, 2] x[3, 3]",
      "-1.0 x[1, 3] x[3, 1]^2 + x[1, 1] x[3, 1] x[3, 3]",
      "-1.0 x[1, 2] x[3, 1] x[3, 2] + x[1, 1] x[3, 2]^2",
      "-1.0 x[1, 2] x[3, 1]^2 + x[1, 1] x[3, 1] x[3, 2]"
      ],
      [
      "-1.224744871391589 x[1, 3] x[2, 2] x[3, 3] + 1.224744871391589 x[1, 2] x[2, 3] x[3, 3]",
      "-1.224744871391589 x[1, 3] x[2, 1] x[3, 3] + 1.224744871391589 x[1, 1] x[2, 3] x[3, 3]",
      "0.5 x[1, 3] x[2, 2] x[3, 1] - 0.5 x[1, 2] x[2, 3] x[3, 1] - 0.5 x[1, 3] x[2, 1] x[3, 2] + 0.5 x[1, 1] x[2, 3] x[3, 2] - 1.0 x[1, 2] x[2, 1] x[3, 3] + x[1, 1] x[2, 2] x[3, 3]",
      "-1.224744871391589 x[1, 3] x[2, 2] x[3, 2] + 1.224744871391589 x[1, 2] x[2, 3] x[3, 2]",
      "-0.8660254037844386 x[1, 3] x[2, 2] x[3, 1] + 0.8660254037844386 x[1, 2] x[2, 3] x[3, 1] - 0.8660254037844386 x[1, 3] x[2, 1] x[3, 2] + 0.8660254037844386 x[1, 1] x[2, 3] x[3, 2]",
      "-1.224744871391589 x[1, 3] x[2, 1] x[3, 1] + 1.224744871391589 x[1, 1] x[2, 3] x[3, 1]",
      "-1.224744871391589 x[1, 2] x[2, 1] x[3, 2] + 1.224744871391589 x[1, 1] x[2, 2] x[3, 2]",
      "-1.224744871391589 x[1, 2] x[2, 1] x[3, 1] + 1.224744871391589 x[1, 1] x[2, 2] x[3, 1]"
      ],
      [
      "-1.0 x[2, 3]^2 x[3, 2] + x[2, 2] x[2, 3] x[3, 3]",
      "-1.0 x[2, 3]^2 x[3, 1] + x[2, 1] x[2, 3] x[3, 3]",
      "-1.224744871391589 x[2, 2] x[2, 3] x[3, 1] + 1.224744871391589 x[2, 1] x[2, 3] x[3, 2]",
      "-1.0 x[2, 2] x[2, 3] x[3, 2] + x[2, 2]^2 x[3, 3]",
      "-0.7071067811865475 x[2, 2] x[2, 3] x[3, 1] - 0.7071067811865475 x[2, 1] x[2, 3] x[3, 2] + 1.4142135623730951 x[2, 1] x[2, 2] x[3, 3]",
      "-1.0 x[2, 1] x[2, 3] x[3, 1] + x[2, 1]^2 x[3, 3]",
      "-1.0 x[2, 2]^2 x[3, 1] + x[2, 1] x[2, 2] x[3, 2]",
      "-1.0 x[2, 1] x[2, 2] x[3, 1] + x[2, 1]^2 x[3, 2]"
      ],
      [
      "-1.4142135623730951 x[1, 3] x[2, 3] x[3, 2] + 0.7071067811865475 x[1, 3] x[2, 2] x[3, 3] + 0.7071067811865475 x[1, 2] x[2, 3] x[3, 3]",
      "-1.4142135623730951 x[1, 3] x[2, 3] x[3, 1] + 0.7071067811865475 x[1, 3] x[2, 1] x[3, 3] + 0.7071067811865475 x[1, 1] x[2, 3] x[3, 3]",
      "-0.8660254037844386 x[1, 3] x[2, 2] x[3, 1] - 0.8660254037844386 x[1, 2] x[2, 3] x[3, 1] + 0.8660254037844386 x[1, 3] x[2, 1] x[3, 2] + 0.8660254037844386 x[1, 1] x[2, 3] x[3, 2]",
      "-0.7071067811865475 x[1, 3] x[2, 2] x[3, 2] - 0.7071067811865475 x[1, 2] x[2, 3] x[3, 2] + 1.4142135623730951 x[1, 2] x[2, 2] x[3, 3]",
      "-0.5 x[1, 3] x[2, 2] x[3, 1] - 0.5 x[1, 2] x[2, 3] x[3, 1] - 0.5 x[1, 3] x[2, 1] x[3, 2] - 0.5 x[1, 1] x[2, 3] x[3, 2] + x[1, 2] x[2, 1] x[3, 3] + x[1, 1] x[2, 2] x[3, 3]",
      "-0.7071067811865475 x[1, 3] x[2, 1] x[3, 1] - 0.7071067811865475 x[1, 1] x[2, 3] x[3, 1] + 1.4142135623730951 x[1, 1] x[2, 1] x[3, 3]",
      "-1.4142135623730951 x[1, 2] x[2, 2] x[3, 1] + 0.7071067811865475 x[1, 2] x[2, 1] x[3, 2] + 0.7071067811865475 x[1, 1] x[2, 2] x[3, 2]",
      "-0.7071067811865475 x[1, 2] x[2, 1] x[3, 1] - 0.7071067811865475 x[1, 1] x[2, 2] x[3, 1] + 1.4142135623730951 x[1, 1] x[2, 1] x[3, 2]"
      ],
      [
      "-1.0 x[1, 3]^2 x[3, 2] + x[1, 2] x[1, 3] x[3, 3]",
      "-1.0 x[1, 3]^2 x[3, 1] + x[1, 1] x[1, 3] x[3, 3]",
      "-1.224744871391589 x[1, 2] x[1, 3] x[3, 1] + 1.224744871391589 x[1, 1] x[1, 3] x[3, 2]",
      "-1.0 x[1, 2] x[1, 3] x[3, 2] + x[1, 2]^2 x[3, 3]",
      "-0.7071067811865475 x[1, 2] x[1, 3] x[3, 1] - 0.7071067811865475 x[1, 1] x[1, 3] x[3, 2] + 1.4142135623730951 x[1, 1] x[1, 2] x[3, 3]",
      "-1.0 x[1, 1] x[1, 3] x[3, 1] + x[1, 1]^2 x[3, 3]",
      "-1.0 x[1, 2]^2 x[3, 1] + x[1, 1] x[1, 2] x[3, 2]",
      "-1.0 x[1, 1] x[1, 2] x[3, 1] + x[1, 1]^2 x[3, 2]"
      ],
      [
      "-1.0 x[1, 3] x[2, 2] x[2, 3] + x[1, 2] x[2, 3]^2",
      "-1.0 x[1, 3] x[2, 1] x[2, 3] + x[1, 1] x[2, 3]^2",
      "-1.224744871391589 x[1, 2] x[2, 1] x[2, 3] + 1.224744871391589 x[1, 1] x[2, 2] x[2, 3]",
      "-1.0 x[1, 3] x[2, 2]^2 + x[1, 2] x[2, 2] x[2, 3]",
      "-1.4142135623730951 x[1, 3] x[2, 1] x[2, 2] + 0.7071067811865475 x[1, 2] x[2, 1] x[2, 3] + 0.7071067811865475 x[1, 1] x[2, 2] x[2, 3]",
      "-1.0 x[1, 3] x[2, 1]^2 + x[1, 1] x[2, 1] x[2, 3]",
      "-1.0 x[1, 2] x[2, 1] x[2, 2] + x[1, 1] x[2, 2]^2",
      "-1.0 x[1, 2] x[2, 1]^2 + x[1, 1] x[2, 1] x[2, 2]"
      ],
      [
      "-1.0 x[1, 3]^2 x[2, 2] + x[1, 2] x[1, 3] x[2, 3]",
      "-1.0 x[1, 3]^2 x[2, 1] + x[1, 1] x[1, 3] x[2, 3]",
      "-1.224744871391589 x[1, 2] x[1, 3] x[2, 1] + 1.224744871391589 x[1, 1] x[1, 3] x[2, 2]",
      "-1.0 x[1, 2] x[1, 3] x[2, 2] + x[1, 2]^2 x[2, 3]",
      "-0.7071067811865475 x[1, 2] x[1, 3] x[2, 1] - 0.7071067811865475 x[1, 1] x[1, 3] x[2, 2] + 1.4142135623730951 x[1, 1] x[1, 2] x[2, 3]",
      "-1.0 x[1, 1] x[1, 3] x[2, 1] + x[1, 1]^2 x[2, 3]",
      "-1.0 x[1, 2]^2 x[2, 1] + x[1, 1] x[1, 2] x[2, 2]",
      "-1.0 x[1, 1] x[1, 2] x[2, 1] + x[1, 1]^2 x[2, 2]"
      ]
      ];

    bolsa_estados = [ t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8 ]
    output_test= []
    for (idx,edox) in enumerate(bolsa_estados), (idy,edoy) in enumerate(bolsa_estados)
      edo_mma = mma_states[idx][idy] |> mma_to_julia
      edo_julia = group_function([2,1,0],edox, edoy) |> expand
      edo_final = expand(edo_mma - edo_julia)
      if !es_cero(edo_final)
        @show edo_mma, edo_julia
      end
      push!(output_test, es_cero(edo_final))
    end
    @test all(output_test)
end
@testset "irrep 21 de SU(4)" begin
    # placeholder
    t_1 = YoungTableau([2,1])# edo 21
    fill!(t_1, [3, 4, 4])
    t_2 = YoungTableau([2,1])# edo 69
    fill!(t_2, [2, 4, 4])
    t_3 = YoungTableau([2,1])# edo 69
    fill!(t_3, [1, 4, 4])
    t_4 = YoungTableau([2,1])# edo 69
    fill!(t_4, [2, 4, 3])
    t_5 = YoungTableau([2,1])# edo 21
    fill!(t_5, [1, 4, 3])
    t_6 = YoungTableau([2,1])# edo 69
    fill!(t_6, [1, 4, 2])
    t_7 = YoungTableau([2,1])# edo 69
    fill!(t_7, [3, 3, 4])
    t_8 = YoungTableau([2,1])# edo 69
    fill!(t_8, [2, 3, 4])

    basura_ejemplo = [
    [
    "-1.0 x[3, 4] x[4, 3] x[4, 4] + x[3, 3] x[4, 4]^2",
    "-1.0 x[3, 4] x[4, 2] x[4, 4] + x[3, 2] x[4, 4]^2",
    "-1.0 x[3, 4] x[4, 1] x[4, 4] + x[3, 1] x[4, 4]^2",
    "-1.224744871391589 x[3, 3] x[4, 2] x[4, 4] + 1.224744871391589 x[3, 2] x[4, 3] x[4, 4]",
    "-1.224744871391589 x[3, 3] x[4, 1] x[4, 4] + 1.224744871391589 x[3, 1] x[4, 3] x[4, 4]",
    "-1.224744871391589 x[3, 2] x[4, 1] x[4, 4] + 1.224744871391589 x[3, 1] x[4, 2] x[4, 4]",
    "-1.0 x[3, 4] x[4, 3]^2 + x[3, 3] x[4, 3] x[4, 4]",
    "-1.4142135623730951 x[3, 4] x[4, 2] x[4, 3] + 0.7071067811865475 x[3, 3] x[4, 2] x[4, 4] + 0.7071067811865475 x[3, 2] x[4, 3] x[4, 4]"
    ],
    [
    "-1.0 x[2, 4] x[4, 3] x[4, 4] + x[2, 3] x[4, 4]^2",
    "-1.0 x[2, 4] x[4, 2] x[4, 4] + x[2, 2] x[4, 4]^2",
    "-1.0 x[2, 4] x[4, 1] x[4, 4] + x[2, 1] x[4, 4]^2",
    "-1.224744871391589 x[2, 3] x[4, 2] x[4, 4] + 1.224744871391589 x[2, 2] x[4, 3] x[4, 4]",
    "-1.224744871391589 x[2, 3] x[4, 1] x[4, 4] + 1.224744871391589 x[2, 1] x[4, 3] x[4, 4]",
    "-1.224744871391589 x[2, 2] x[4, 1] x[4, 4] + 1.224744871391589 x[2, 1] x[4, 2] x[4, 4]",
    "-1.0 x[2, 4] x[4, 3]^2 + x[2, 3] x[4, 3] x[4, 4]",
    "-1.4142135623730951 x[2, 4] x[4, 2] x[4, 3] + 0.7071067811865475 x[2, 3] x[4, 2] x[4, 4] + 0.7071067811865475 x[2, 2] x[4, 3] x[4, 4]"
    ],
    [
    "-1.0 x[1, 4] x[4, 3] x[4, 4] + x[1, 3] x[4, 4]^2",
    "-1.0 x[1, 4] x[4, 2] x[4, 4] + x[1, 2] x[4, 4]^2",
    "-1.0 x[1, 4] x[4, 1] x[4, 4] + x[1, 1] x[4, 4]^2",
    "-1.224744871391589 x[1, 3] x[4, 2] x[4, 4] + 1.224744871391589 x[1, 2] x[4, 3] x[4, 4]",
    "-1.224744871391589 x[1, 3] x[4, 1] x[4, 4] + 1.224744871391589 x[1, 1] x[4, 3] x[4, 4]",
    "-1.224744871391589 x[1, 2] x[4, 1] x[4, 4] + 1.224744871391589 x[1, 1] x[4, 2] x[4, 4]",
    "-1.0 x[1, 4] x[4, 3]^2 + x[1, 3] x[4, 3] x[4, 4]",
    "-1.4142135623730951 x[1, 4] x[4, 2] x[4, 3] + 0.7071067811865475 x[1, 3] x[4, 2] x[4, 4] + 0.7071067811865475 x[1, 2] x[4, 3] x[4, 4]"
    ],
    [
    "-1.224744871391589 x[2, 4] x[3, 3] x[4, 4] + 1.224744871391589 x[2, 3] x[3, 4] x[4, 4]",
    "-1.224744871391589 x[2, 4] x[3, 2] x[4, 4] + 1.224744871391589 x[2, 2] x[3, 4] x[4, 4]",
    "-1.224744871391589 x[2, 4] x[3, 1] x[4, 4] + 1.224744871391589 x[2, 1] x[3, 4] x[4, 4]",
    "0.5 x[2, 4] x[3, 3] x[4, 2] - 0.5 x[2, 3] x[3, 4] x[4, 2] - 0.5 x[2, 4] x[3, 2] x[4, 3] + 0.5 x[2, 2] x[3, 4] x[4, 3] - 1.0 x[2, 3] x[3, 2] x[4, 4] + x[2, 2] x[3, 3] x[4, 4]",
    "0.5 x[2, 4] x[3, 3] x[4, 1] - 0.5 x[2, 3] x[3, 4] x[4, 1] - 0.5 x[2, 4] x[3, 1] x[4, 3] + 0.5 x[2, 1] x[3, 4] x[4, 3] - 1.0 x[2, 3] x[3, 1] x[4, 4] + x[2, 1] x[3, 3] x[4, 4]",
    "0.5 x[2, 4] x[3, 2] x[4, 1] - 0.5 x[2, 2] x[3, 4] x[4, 1] - 0.5 x[2, 4] x[3, 1] x[4, 2] + 0.5 x[2, 1] x[3, 4] x[4, 2] - 1.0 x[2, 2] x[3, 1] x[4, 4] + x[2, 1] x[3, 2] x[4, 4]",
    "-1.224744871391589 x[2, 4] x[3, 3] x[4, 3] + 1.224744871391589 x[2, 3] x[3, 4] x[4, 3]",
    "-0.8660254037844386 x[2, 4] x[3, 3] x[4, 2] + 0.8660254037844386 x[2, 3] x[3, 4] x[4, 2] - 0.8660254037844386 x[2, 4] x[3, 2] x[4, 3] + 0.8660254037844386 x[2, 2] x[3, 4] x[4, 3]"
    ],
    [
    "-1.224744871391589 x[1, 4] x[3, 3] x[4, 4] + 1.224744871391589 x[1, 3] x[3, 4] x[4, 4]",
    "-1.224744871391589 x[1, 4] x[3, 2] x[4, 4] + 1.224744871391589 x[1, 2] x[3, 4] x[4, 4]",
    "-1.224744871391589 x[1, 4] x[3, 1] x[4, 4] + 1.224744871391589 x[1, 1] x[3, 4] x[4, 4]",
    "0.5 x[1, 4] x[3, 3] x[4, 2] - 0.5 x[1, 3] x[3, 4] x[4, 2] - 0.5 x[1, 4] x[3, 2] x[4, 3] + 0.5 x[1, 2] x[3, 4] x[4, 3] - 1.0 x[1, 3] x[3, 2] x[4, 4] + x[1, 2] x[3, 3] x[4, 4]",
    "0.5 x[1, 4] x[3, 3] x[4, 1] - 0.5 x[1, 3] x[3, 4] x[4, 1] - 0.5 x[1, 4] x[3, 1] x[4, 3] + 0.5 x[1, 1] x[3, 4] x[4, 3] - 1.0 x[1, 3] x[3, 1] x[4, 4] + x[1, 1] x[3, 3] x[4, 4]",
    "0.5 x[1, 4] x[3, 2] x[4, 1] - 0.5 x[1, 2] x[3, 4] x[4, 1] - 0.5 x[1, 4] x[3, 1] x[4, 2] + 0.5 x[1, 1] x[3, 4] x[4, 2] - 1.0 x[1, 2] x[3, 1] x[4, 4] + x[1, 1] x[3, 2] x[4, 4]",
    "-1.224744871391589 x[1, 4] x[3, 3] x[4, 3] + 1.224744871391589 x[1, 3] x[3, 4] x[4, 3]",
    "-0.8660254037844386 x[1, 4] x[3, 3] x[4, 2] + 0.8660254037844386 x[1, 3] x[3, 4] x[4, 2] - 0.8660254037844386 x[1, 4] x[3, 2] x[4, 3] + 0.8660254037844386 x[1, 2] x[3, 4] x[4, 3]"
    ],
    [
    "-1.224744871391589 x[1, 4] x[2, 3] x[4, 4] + 1.224744871391589 x[1, 3] x[2, 4] x[4, 4]",
    "-1.224744871391589 x[1, 4] x[2, 2] x[4, 4] + 1.224744871391589 x[1, 2] x[2, 4] x[4, 4]",
    "-1.224744871391589 x[1, 4] x[2, 1] x[4, 4] + 1.224744871391589 x[1, 1] x[2, 4] x[4, 4]",
    "0.5 x[1, 4] x[2, 3] x[4, 2] - 0.5 x[1, 3] x[2, 4] x[4, 2] - 0.5 x[1, 4] x[2, 2] x[4, 3] + 0.5 x[1, 2] x[2, 4] x[4, 3] - 1.0 x[1, 3] x[2, 2] x[4, 4] + x[1, 2] x[2, 3] x[4, 4]",
    "0.5 x[1, 4] x[2, 3] x[4, 1] - 0.5 x[1, 3] x[2, 4] x[4, 1] - 0.5 x[1, 4] x[2, 1] x[4, 3] + 0.5 x[1, 1] x[2, 4] x[4, 3] - 1.0 x[1, 3] x[2, 1] x[4, 4] + x[1, 1] x[2, 3] x[4, 4]",
    "0.5 x[1, 4] x[2, 2] x[4, 1] - 0.5 x[1, 2] x[2, 4] x[4, 1] - 0.5 x[1, 4] x[2, 1] x[4, 2] + 0.5 x[1, 1] x[2, 4] x[4, 2] - 1.0 x[1, 2] x[2, 1] x[4, 4] + x[1, 1] x[2, 2] x[4, 4]",
    "-1.224744871391589 x[1, 4] x[2, 3] x[4, 3] + 1.224744871391589 x[1, 3] x[2, 4] x[4, 3]",
    "-0.8660254037844386 x[1, 4] x[2, 3] x[4, 2] + 0.8660254037844386 x[1, 3] x[2, 4] x[4, 2] - 0.8660254037844386 x[1, 4] x[2, 2] x[4, 3] + 0.8660254037844386 x[1, 2] x[2, 4] x[4, 3]"
    ],
    [
    "-1.0 x[3, 4]^2 x[4, 3] + x[3, 3] x[3, 4] x[4, 4]",
    "-1.0 x[3, 4]^2 x[4, 2] + x[3, 2] x[3, 4] x[4, 4]",
    "-1.0 x[3, 4]^2 x[4, 1] + x[3, 1] x[3, 4] x[4, 4]",
    "-1.224744871391589 x[3, 3] x[3, 4] x[4, 2] + 1.224744871391589 x[3, 2] x[3, 4] x[4, 3]",
    "-1.224744871391589 x[3, 3] x[3, 4] x[4, 1] + 1.224744871391589 x[3, 1] x[3, 4] x[4, 3]",
    "-1.224744871391589 x[3, 2] x[3, 4] x[4, 1] + 1.224744871391589 x[3, 1] x[3, 4] x[4, 2]",
    "-1.0 x[3, 3] x[3, 4] x[4, 3] + x[3, 3]^2 x[4, 4]",
    "-0.7071067811865475 x[3, 3] x[3, 4] x[4, 2] - 0.7071067811865475 x[3, 2] x[3, 4] x[4, 3] + 1.4142135623730951 x[3, 2] x[3, 3] x[4, 4]"
    ],
    [
    "-1.4142135623730951 x[2, 4] x[3, 4] x[4, 3] + 0.7071067811865475 x[2, 4] x[3, 3] x[4, 4] + 0.7071067811865475 x[2, 3] x[3, 4] x[4, 4]",
    "-1.4142135623730951 x[2, 4] x[3, 4] x[4, 2] + 0.7071067811865475 x[2, 4] x[3, 2] x[4, 4] + 0.7071067811865475 x[2, 2] x[3, 4] x[4, 4]",
    "-1.4142135623730951 x[2, 4] x[3, 4] x[4, 1] + 0.7071067811865475 x[2, 4] x[3, 1] x[4, 4] + 0.7071067811865475 x[2, 1] x[3, 4] x[4, 4]",
    "-0.8660254037844386 x[2, 4] x[3, 3] x[4, 2] - 0.8660254037844386 x[2, 3] x[3, 4] x[4, 2] + 0.8660254037844386 x[2, 4] x[3, 2] x[4, 3] + 0.8660254037844386 x[2, 2] x[3, 4] x[4, 3]",
    "-0.8660254037844386 x[2, 4] x[3, 3] x[4, 1] - 0.8660254037844386 x[2, 3] x[3, 4] x[4, 1] + 0.8660254037844386 x[2, 4] x[3, 1] x[4, 3] + 0.8660254037844386 x[2, 1] x[3, 4] x[4, 3]",
    "-0.8660254037844386 x[2, 4] x[3, 2] x[4, 1] - 0.8660254037844386 x[2, 2] x[3, 4] x[4, 1] + 0.8660254037844386 x[2, 4] x[3, 1] x[4, 2] + 0.8660254037844386 x[2, 1] x[3, 4] x[4, 2]",
    "-0.7071067811865475 x[2, 4] x[3, 3] x[4, 3] - 0.7071067811865475 x[2, 3] x[3, 4] x[4, 3] + 1.4142135623730951 x[2, 3] x[3, 3] x[4, 4]",
    "-0.5 x[2, 4] x[3, 3] x[4, 2] - 0.5 x[2, 3] x[3, 4] x[4, 2] - 0.5 x[2, 4] x[3, 2] x[4, 3] - 0.5 x[2, 2] x[3, 4] x[4, 3] + x[2, 3] x[3, 2] x[4, 4] + x[2, 2] x[3, 3] x[4, 4]"
    ]
    ];
    bolsa_estados = [ t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8 ]
    output_test = []
    for (idx,edox) in enumerate(bolsa_estados), (idy,edoy) in enumerate(bolsa_estados)
      edo_mma = basura_ejemplo[idx][idy] |> mma_to_julia
      edo_julia = group_function([2,1,0,0],edox, edoy) |> expand
      edo_final = expand(edo_mma - edo_julia)
      if !es_cero(edo_final)
        @show edo_mma, edo_julia
      end
      push!(output_test, es_cero(edo_final) )
    end
    @test all(output_test)
end

@testset "Comparison with immanant 210" begin
    mat = rand(Haar(2), 3)
    pt_1 = GTPattern([[2,1,0], [2,0],[1]],[1])
    pt_2 = GTPattern([[2,1,0], [1,1],[1]],[1])
    suma::Complex{Float64} = group_function([2,1,0], pt_1, pt_1, mat)+ group_function([2,1,0], pt_2, pt_2, mat)
    @test suma ≈ immanant210(mat)
end

@testset "Comparison with immanant 2110" begin
    part = [2,1,1,0]
    zeroweightstates = findzero(part)

    mat = rand(Haar(2), sum(part))

    total::Complex{Float64} = sum(group_function(part, p, p, mat) for p in zeroweightstates)

    @test isapprox(abs(immanant2110(mat) -total)^2, 0.0, atol = 10^(-6))
end

@testset "Sum rules 3x3" begin
    welcome = basis_states([2,0,0]);

    α1,β1,γ1 = rand(Float64,3)
    xx=su2_block(3,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=su2_block(3,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=su2_block(3,1,(α3,β3,γ3))

    mat = xx*yy*zz;
    rate1 = abs( group_function([2,1,0], welcome[5], welcome[3], mat) )^2 + abs( group_function([2,1,0], welcome[5], welcome[5], mat) )^2
    matc1 = yy*zz;
    rate2 = abs( group_function([2,1,0], welcome[5], welcome[3], matc1) )^2 + abs( group_function([2,1,0], welcome[5], welcome[5], matc1) )^2
    @test rate1 ≈ rate2
end

@testset "Sum rules 4x4" begin
    α1,β1,γ1 = rand(Float64,3)
    xx=su2_block(4,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=su2_block(4,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=su2_block(4,1,(α3,β3,γ3))
    α4,β4 = rand(Float64,3)
    xx2=su2_block(4,3,(α4,β4,α4))
    α5,β5 = rand(Float64,2)
    yy2=su2_block(4,2,(α5,β5,α5))
    α6,β6,γ6 = rand(Float64,3)
    zz2=su2_block(4,1,(α6,β6,γ6))

    welcome = basis_states([2,0,0,0])

    mat4 = xx*yy*zz
    rate1 = abs( group_function([2,0,0,0], welcome[9], welcome[9], mat4) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[4], mat4) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[7], mat4) )^2
    mat4c1 = xx*yy*zz*xx2*yy2
    rate2 = abs( group_function([2,0,0,0], welcome[9], welcome[9], mat4c1) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[4], mat4c1) )^2 + abs( group_function([2,0,0,0], welcome[9], welcome[7], mat4c1) )^2
    mat4 = xx*yy*zz*xx2*yy2*zz2
    @test rate1 ≈ rate2
end

@testset "Comparison with permanent of a submatrix" begin
    α1,β1,γ1 = rand(Float64,3)
    xx=su2_block(4,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=su2_block(4,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=su2_block(4,1,(α3,β3,γ3))
    α4,β4 = rand(Float64,3)
    xx2=su2_block(4,3,(α4,β4,α4))
    α5,β5 = rand(Float64,2)
    yy2=su2_block(4,2,(α5,β5,α5))
    α6,β6,γ6 = rand(Float64,3)
    zz2=su2_block(4,1,(α6,β6,γ6))

    mat4 = xx*yy*zz*xx2*yy2*zz2
    welcome = basis_states([3,0,0,0])

    edox = filter(x -> pweight(x) == [0,1,1,1] , welcome)[1]
    edoy = filter(x -> pweight(x) == [1,0,2,0], welcome)[1]

    @test group_function([3,0,0,0], edox, edoy, mat4) ≈ permanent(mat4[[1,2,3], [2,2,4]])/sqrt(2)
end
    # @test group_function([3,0,0,0], edox, edoy, mat4) ≈ immanant(Partition([3]), mat4[[1,2,3], [2,2,4]])/sqrt(2)

@testset "Testing labelling for construction of unitary matrices" begin
  α1,β1,γ1 = rand(Float64,3)
  xx=su2_block(4,1,(α1,β1,γ1))
  α2,β2 = rand(Float64,2)
  yy=su2_block(4,2,(α2,β2,α2))
  α3,β3,γ3 = rand(Float64,3)
  zz=su2_block(4,1,(α3,β3,γ3))
  α4,β4 = rand(Float64,3)
  xx2=su2_block(4,3,(α4,β4,α4))
  α5,β5 = rand(Float64,2)
  yy2=su2_block(4,2,(α5,β5,α5))
  α6,β6,γ6 = rand(Float64,3)
  zz2=su2_block(4,1,(α6,β6,γ6))

  matsimple = sud_from_angles([α1,β1,γ1,α2,β2,α3,β3,γ3,α4,β4,α5,β5,α6,β6,γ6], 4)
  matsimple_quotient = sud_from_angles([α1,β1,γ1,α2,β2,α3,β3,γ3,α4,β4,α5,β5,α6,β6,γ6], 4; quotient = true)
  mat =  zz2*yy2*xx2*zz*yy*xx 
  mat_quotient =  zz2*yy2*xx2#*zz*yy*xx 

  @test norm(matsimple - mat) < 10^(-11)
  @test norm(matsimple_quotient - mat_quotient) < 10^(-11)
end

@testset "Testing labelling for construction of unitary matrices" begin
    α1,β1,γ1 = rand(Float64,3)
    xx=su2_block(3,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=su2_block(3,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=su2_block(3,1,(α3,β3,γ3))
    mat = sud_from_angles([α1,β1,γ1, α2,β2, α3,β3,γ3 ], 3)
    mat2 = zz*yy*xx;
    @test norm(mat-mat2) < 10^(-5)

end

@testset "Identity irrep 100 for U(3)" begin
  mat = rand(Haar(2),3)
  kak = basis_states([1,0,0])
  @test all(isapprox.(0, [group_function([1,0,0], x, y, mat) - mat[end-i+1, end-j+1] for (i,x) in enumerate(kak), (j,y) in enumerate(kak)], atol=1e-10))
end

@testset "Identity irrep 1000 for U(4)" begin
  mat = rand(Haar(2),4)
  irrep = [1,0,0,0]
  kak = basis_states(irrep)
  @test all(isapprox.(0, [group_function(irrep, x, y, mat) - mat[end-i+1, end-j+1] for (i,x) in enumerate(kak), (j,y) in enumerate(kak)], atol=1e-10))
end

@testset "Symmetry with respect to hermitian conjugate 210" begin
  mat = rand(Haar(2),3)
  invmat = inv(mat)
  irrep = [2,1,0]
  kak = basis_states(irrep)
  @test all(isapprox.(0, [group_function(irrep, x, y, mat) - conj(group_function(irrep, y, x, invmat)) for (i,x) in enumerate(kak), (j,y) in enumerate(kak)], atol=1e-10))
end

@testset "Symmetry with respect to hermitian conjugate 221" begin
  mat = rand(Haar(2),3)
  invmat = inv(mat)
  irrep = [2,2,1]
  kak = basis_states(irrep)
  @test all(isapprox.(0, [group_function(irrep, x, y, mat) - conj(group_function(irrep, y, x, invmat)) for (i,x) in enumerate(kak), (j,y) in enumerate(kak)], atol=1e-10))
end

@testset "group_function variant consistency" begin
    # Test 1: All entry points produce same symbolic result
    @testset "Symbolic variants consistency" begin
        λ = [2, 1, 0]  # 8-dim irrep of SU(3)
        patterns = basis_states(λ)
        pat_u, pat_v = patterns[1], patterns[2]
        tab_u, tab_v = YoungTableau(pat_u), YoungTableau(pat_v)

        # All symbolic variants should give same polynomial
        result_tab = group_function(λ, tab_u, tab_v)           # YTableau, no mat
        result_pat = group_function(λ, pat_u, pat_v)           # GTPattern, no mat

        @test es_cero(expand(result_tab - result_pat))
    end

    # Test 2: Symbolic with explicit symbolic_matrix equals no-matrix version
    @testset "Symbolic matrix vs implicit symbolic" begin
        λ = [1, 1, 0]  # 3-dim irrep of SU(3)
        patterns = basis_states(λ)
        pat_u, pat_v = patterns[1], patterns[2]
        tab_u, tab_v = YoungTableau(pat_u), YoungTableau(pat_v)

        sym_mat = GroupFunctions.symbolic_matrix(3)
        result_implicit = group_function(λ, tab_u, tab_v)
        result_explicit = group_function(λ, tab_u, tab_v, sym_mat)

        @test es_cero(expand(result_implicit - result_explicit))
    end

    # Test 3: Numeric evaluation matches symbolic substituted
    @testset "Numeric vs substituted symbolic" begin
        λ = [2, 0]  # 3-dim irrep of SU(2)
        mat = rand(Haar(2), 2)
        mat = mat / det(mat)^(1/2)

        patterns = basis_states(λ)
        pat_u, pat_v = patterns[1], patterns[2]

        # Numeric directly
        numeric_result = group_function(λ, pat_u, pat_v, mat)

        # Symbolic then substitute
        symbolic_result = group_function(λ, pat_u, pat_v)
        substituted = symbolic_result
        for i in 1:2, j in 1:2
            substituted = SymEngine.subs(substituted,
                SymEngine.symbols("u_$(i)_$(j)"), mat[i,j])
        end
        #@show numeric_result
        #@show ComplexF64(N(substituted))
        @test numeric_result ≈ convert(ComplexF64,substituted)
    end

    # Test 4: Batch vs individual calls
    @testset "Batch vs individual consistency" begin
        λ = [2, 1, 0]
        mat = rand(Haar(2), 3)
        mat = mat / det(mat)^(1/3)

        batch_result, patterns = group_function(λ, mat)

        for (i, pat_u) in enumerate(patterns)
            for (j, pat_v) in enumerate(patterns)
                individual = group_function(λ, pat_u, pat_v, mat)
                @test batch_result[i, j] ≈ individual
            end
        end
    end
end
