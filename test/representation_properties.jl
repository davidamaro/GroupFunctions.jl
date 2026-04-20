@testset "comparison of the character (2,1,0) of SU3" begin
    mat = rand(Haar(2), 3)
    special_unitary_matrix = mat / det(mat)^(1/3)
    lambdas = eigvals(special_unitary_matrix)
    ideal_character = lambdas[1]/lambdas[2] + lambdas[2]/lambdas[1] +
                      lambdas[1]/lambdas[3] + lambdas[3]/lambdas[1] +
                      lambdas[2]/lambdas[3] + lambdas[3]/lambdas[2] + 2
    @test @inferred(GroupFunctions.character([2,1,0], special_unitary_matrix)) ≈ ideal_character
end

@testset "comparison of the character (2,0)" begin
    mat = rand(Haar(2), 3)
    special_unitary_matrix = mat / det(mat)^(1/3)
    lambdas = eigvals(special_unitary_matrix)
    ideal_character = lambdas[1]^2 + lambdas[2]^2 + lambdas[3]^2 +
                      lambdas[1]*lambdas[2] + lambdas[1]*lambdas[3] +
                      lambdas[2]*lambdas[3]
    @test @inferred(GroupFunctions.character([2,0,0], special_unitary_matrix)) ≈ ideal_character
end

@testset "comparison of the character (3,0)" begin
    mat = rand(Haar(2), 3)
    special_unitary_matrix = mat / det(mat)^(1/3)
    lambdas = eigvals(special_unitary_matrix)
    ideal_character = lambdas[1]^3 + lambdas[2]^3 + lambdas[3]^3 +
                      lambdas[1]^2*lambdas[2] + lambdas[1]^2*lambdas[3] +
                      lambdas[2]^2*lambdas[1] + lambdas[2]^2*lambdas[3] +
                      lambdas[3]^2*lambdas[1] + lambdas[3]^2*lambdas[2] +
                      lambdas[1]*lambdas[2]*lambdas[3]
    @test @inferred(GroupFunctions.character([3,0,0], special_unitary_matrix)) ≈ ideal_character
end

@testset "SU(2) irrep homomorphism" begin
    left_factor = rand(Haar(2), 2)
    left_factor = left_factor / det(left_factor)^(1 / 2)
    right_factor = rand(Haar(2), 2)
    right_factor = right_factor / det(right_factor)^(1 / 2)
    product_factor = left_factor * right_factor

    for λ in ([2, 0], [3, 0])
        left_representation, _ = group_function(λ, left_factor)
        right_representation, _ = group_function(λ, right_factor)
        product_representation, _ = group_function(λ, product_factor)
        @test isapprox(left_representation * right_representation, product_representation; atol=1e-8, rtol=1e-8)
    end
end

@testset "SU(3) irrep homomorphism" begin
    left_factor = rand(Haar(2), 3)
    left_factor = left_factor / det(left_factor)^(1 / 3)
    right_factor = rand(Haar(2), 3)
    right_factor = right_factor / det(right_factor)^(1 / 3)
    product_factor = left_factor * right_factor

    for λ in ([2, 1, 0], [2, 2, 0])
        left_representation, _ = group_function(λ, left_factor)
        right_representation, _ = group_function(λ, right_factor)
        product_representation, _ = group_function(λ, product_factor)
        @test isapprox(left_representation * right_representation, product_representation; atol=1e-8, rtol=1e-8)
    end
end

@testset "Comparison between Legendre polynomails and D-function" begin
  α1,β1,γ1 = rand(Float64,3)
  leading_mix12 = su2_block(3,1,(α1,β1,γ1))
  α2,β2 = rand(Float64,2)
  middle_mix23 = su2_block(3,2,(α2,β2,α2))
  α3,β3,γ3 = rand(Float64,3)
  trailing_mix12 = su2_block(3,1,(α3,β3,γ3))

  composed_matrix = leading_mix12 * middle_mix23 * trailing_mix12

  linear_legendre_pattern = GTPattern([[2,1,0],[1,1], [1]], [1])
  @test group_function([2,1,0], linear_legendre_pattern, linear_legendre_pattern, composed_matrix) ≈ (1/4)*(1+3*cos(β2))

  quadratic_legendre_pattern = GTPattern([[4,2,0],[2,2], [2]], [2])
  @test group_function([4,2,0], quadratic_legendre_pattern, quadratic_legendre_pattern, composed_matrix) ≈ (1/3^2)*(1+3*cos(β2)+5*(1/2)*(3*cos(β2)^2 - 1))
end

@testset "Identity irrep 100 for U(3)" begin
  mat = rand(Haar(2),3)
  basis_patterns = basis_states([1,0,0])
  actual = [group_function([1,0,0], x, y, mat) for x in basis_patterns, y in basis_patterns]
  expected = [mat[end-i+1, end-j+1] for (i, _) in enumerate(basis_patterns), (j, _) in enumerate(basis_patterns)]
  @test isapprox(actual, expected; atol=1e-10, rtol=0)
end

@testset "Identity irrep 1000 for U(4)" begin
  mat = rand(Haar(2),4)
  irrep = [1,0,0,0]
  basis_patterns = basis_states(irrep)
  actual = [group_function(irrep, x, y, mat) for x in basis_patterns, y in basis_patterns]
  expected = [mat[end-i+1, end-j+1] for (i, _) in enumerate(basis_patterns), (j, _) in enumerate(basis_patterns)]
  @test isapprox(actual, expected; atol=1e-10, rtol=0)
end

@testset "Symmetry with respect to hermitian conjugate 210" begin
  mat = rand(Haar(2),3)
  invmat = inv(mat)
  irrep = [2,1,0]
  basis_patterns = basis_states(irrep)
  direct_values = [group_function(irrep, x, y, mat) for x in basis_patterns, y in basis_patterns]
  adjoint_symmetry_values = [conj(group_function(irrep, y, x, invmat)) for x in basis_patterns, y in basis_patterns]
  @test isapprox(direct_values, adjoint_symmetry_values; atol=1e-10, rtol=0)
end

@testset "Symmetry with respect to hermitian conjugate 221" begin
  mat = rand(Haar(2),3)
  invmat = inv(mat)
  irrep = [2,2,1]
  basis_patterns = basis_states(irrep)
  direct_values = [group_function(irrep, x, y, mat) for x in basis_patterns, y in basis_patterns]
  adjoint_symmetry_values = [conj(group_function(irrep, y, x, invmat)) for x in basis_patterns, y in basis_patterns]
  @test isapprox(direct_values, adjoint_symmetry_values; atol=1e-10, rtol=0)
end
