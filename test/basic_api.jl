@testset "Number of zero weight states." begin
  solutions = find_tablaeux_fillings([1, 1, 2, 1, 0], [1, 1, 2, 1, 0])
  @test length(solutions) == 33
  solutions = find_tablaeux_fillings([2,2,2,0,0,0], [2,2,2,0,0,0])
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

@testset "group_function variant consistency" begin
    # Test 1: All entry points produce same symbolic result
    @testset "Symbolic variants consistency" begin
        λ = [2, 1, 0]  # 8-dim irrep of SU(3)
        patterns = basis_states(λ)
        first_pattern, second_pattern = patterns[1], patterns[2]
        first_tableau, second_tableau = YoungTableau(first_pattern), YoungTableau(second_pattern)

        # All symbolic variants should give same polynomial
        result_tab = group_function(λ, first_tableau, second_tableau)   # YTableau, no mat
        result_pat = group_function(λ, first_pattern, second_pattern)   # GTPattern, no mat

        @test symbolic_isapprox(result_tab, result_pat)
    end

    # Test 2: Symbolic with explicit symbolic_matrix equals no-matrix version
    @testset "Symbolic matrix vs implicit symbolic" begin
        λ = [1, 1, 0]  # 3-dim irrep of SU(3)
        patterns = basis_states(λ)
        first_pattern, second_pattern = patterns[1], patterns[2]
        first_tableau, second_tableau = YoungTableau(first_pattern), YoungTableau(second_pattern)

        symbolic_matrix = GroupFunctions.symbolic_matrix(3)
        result_implicit = group_function(λ, first_tableau, second_tableau)
        result_explicit = group_function(λ, first_tableau, second_tableau, symbolic_matrix)

        @test symbolic_isapprox(result_implicit, result_explicit)
    end

    # Test 3: Numeric evaluation matches symbolic substituted
    @testset "Numeric vs substituted symbolic" begin
        λ = [2, 0]  # 3-dim irrep of SU(2)
        mat = rand(Haar(2), 2)
        mat = mat / det(mat)^(1/2)

        patterns = basis_states(λ)
        first_pattern, second_pattern = patterns[1], patterns[2]

        # Numeric directly
        numeric_result = group_function(λ, first_pattern, second_pattern, mat)

        # Symbolic then substitute
        symbolic_result = group_function(λ, first_pattern, second_pattern)
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

        for (i, first_pattern) in enumerate(patterns)
            for (j, second_pattern) in enumerate(patterns)
                individual = group_function(λ, first_pattern, second_pattern, mat)
                @test batch_result[i, j] ≈ individual
            end
        end
    end
end
