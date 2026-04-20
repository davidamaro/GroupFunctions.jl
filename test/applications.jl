@testset "Comparison with immanant 210" begin
    mat = rand(Haar(2), 3)
    first_zero_weight_pattern = GTPattern([[2,1,0], [2,0],[1]],[1])
    second_zero_weight_pattern = GTPattern([[2,1,0], [1,1],[1]],[1])
    diagonal_sum::Complex{Float64} = group_function([2,1,0], first_zero_weight_pattern, first_zero_weight_pattern, mat) +
                                     group_function([2,1,0], second_zero_weight_pattern, second_zero_weight_pattern, mat)
    @test diagonal_sum ≈ immanant210(mat)
end

@testset "Comparison with immanant 2110" begin
    irrep = [2,1,1,0]
    zero_weight_states = findzero(irrep)

    mat = rand(Haar(2), sum(irrep))

    total::Complex{Float64} = sum(group_function(irrep, pattern, pattern, mat) for pattern in zero_weight_states)

    @test isapprox(immanant2110(mat), total; atol=1e-6, rtol=0)
end

@testset "Sum rules 3x3" begin
    basis_patterns = basis_states([2,0,0])

    α1,β1,γ1 = rand(Float64,3)
    leading_mix12 = su2_block(3,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    middle_mix23 = su2_block(3,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    trailing_mix12 = su2_block(3,1,(α3,β3,γ3))

    reference_pattern = basis_patterns[5]
    first_selected_pattern = basis_patterns[3]
    second_selected_pattern = basis_patterns[5]

    full_matrix = leading_mix12 * middle_mix23 * trailing_mix12
    full_rate = abs(group_function([2,1,0], reference_pattern, first_selected_pattern, full_matrix))^2 +
                abs(group_function([2,1,0], reference_pattern, second_selected_pattern, full_matrix))^2
    truncated_matrix = middle_mix23 * trailing_mix12
    truncated_rate = abs(group_function([2,1,0], reference_pattern, first_selected_pattern, truncated_matrix))^2 +
                     abs(group_function([2,1,0], reference_pattern, second_selected_pattern, truncated_matrix))^2
    @test full_rate ≈ truncated_rate
end

@testset "Sum rules 4x4" begin
    α1,β1,γ1 = rand(Float64,3)
    first_mix12 = su2_block(4,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    first_mix23 = su2_block(4,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    second_mix12 = su2_block(4,1,(α3,β3,γ3))
    α4,β4 = rand(Float64,3)
    mix34 = su2_block(4,3,(α4,β4,α4))
    α5,β5 = rand(Float64,2)
    second_mix23 = su2_block(4,2,(α5,β5,α5))
    α6,β6,γ6 = rand(Float64,3)
    third_mix12 = su2_block(4,1,(α6,β6,γ6))

    basis_patterns = basis_states([2,0,0,0])
    reference_pattern = basis_patterns[9]
    selected_patterns = basis_patterns[[9, 4, 7]]

    initial_matrix = first_mix12 * first_mix23 * second_mix12
    initial_rate = sum(abs(group_function([2,0,0,0], reference_pattern, selected_pattern, initial_matrix))^2
                       for selected_pattern in selected_patterns)
    truncated_matrix = first_mix12 * first_mix23 * second_mix12 * mix34 * second_mix23
    truncated_rate = sum(abs(group_function([2,0,0,0], reference_pattern, selected_pattern, truncated_matrix))^2
                         for selected_pattern in selected_patterns)
    full_matrix = truncated_matrix * third_mix12
    @test initial_rate ≈ truncated_rate
    @test size(full_matrix) == (4, 4)
end

@testset "Comparison with permanent of a submatrix" begin
    α1,β1,γ1 = rand(Float64,3)
    first_mix12 = su2_block(4,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    first_mix23 = su2_block(4,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    second_mix12 = su2_block(4,1,(α3,β3,γ3))
    α4,β4 = rand(Float64,3)
    mix34 = su2_block(4,3,(α4,β4,α4))
    α5,β5 = rand(Float64,2)
    second_mix23 = su2_block(4,2,(α5,β5,α5))
    α6,β6,γ6 = rand(Float64,3)
    third_mix12 = su2_block(4,1,(α6,β6,γ6))

    unitary_matrix = first_mix12 * first_mix23 * second_mix12 * mix34 * second_mix23 * third_mix12
    basis_patterns = basis_states([3,0,0,0])

    first_pattern = filter(pattern -> pweight(pattern) == [0,1,1,1], basis_patterns)[1]
    second_pattern = filter(pattern -> pweight(pattern) == [1,0,2,0], basis_patterns)[1]

    @test group_function([3,0,0,0], first_pattern, second_pattern, unitary_matrix) ≈ permanent(unitary_matrix[[1,2,3], [2,2,4]])/sqrt(2)
end
    # @test group_function([3,0,0,0], first_pattern, second_pattern, unitary_matrix) ≈ immanant(Partition([3]), unitary_matrix[[1,2,3], [2,2,4]])/sqrt(2)

@testset "Testing labelling for construction of unitary matrices" begin
  α1,β1,γ1 = rand(Float64,3)
  first_mix12 = su2_block(4,1,(α1,β1,γ1))
  α2,β2 = rand(Float64,2)
  first_mix23 = su2_block(4,2,(α2,β2,α2))
  α3,β3,γ3 = rand(Float64,3)
  second_mix12 = su2_block(4,1,(α3,β3,γ3))
  α4,β4 = rand(Float64,3)
  mix34 = su2_block(4,3,(α4,β4,α4))
  α5,β5 = rand(Float64,2)
  second_mix23 = su2_block(4,2,(α5,β5,α5))
  α6,β6,γ6 = rand(Float64,3)
  third_mix12 = su2_block(4,1,(α6,β6,γ6))

  constructed_matrix = sud_from_angles([α1,β1,γ1,α2,β2,α3,β3,γ3,α4,β4,α5,β5,α6,β6,γ6], 4)
  constructed_quotient_matrix = sud_from_angles([α1,β1,γ1,α2,β2,α3,β3,γ3,α4,β4,α5,β5,α6,β6,γ6], 4; quotient = true)
  reference_matrix = third_mix12 * second_mix23 * mix34 * second_mix12 * first_mix23 * first_mix12
  reference_quotient_matrix = third_mix12 * second_mix23 * mix34

  @test norm(constructed_matrix - reference_matrix) < 10^(-11)
  @test norm(constructed_quotient_matrix - reference_quotient_matrix) < 10^(-11)
end

@testset "Testing labelling for construction of unitary matrices" begin
    α1,β1,γ1 = rand(Float64,3)
    first_mix12 = su2_block(3,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    mix23 = su2_block(3,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    second_mix12 = su2_block(3,1,(α3,β3,γ3))
    constructed_matrix = sud_from_angles([α1,β1,γ1, α2,β2, α3,β3,γ3 ], 3)
    reference_matrix = second_mix12 * mix23 * first_mix12
    @test norm(constructed_matrix - reference_matrix) < 10^(-5)

end
