function legacy_character(irrep, matrix = nothing)
    patterns = basis_states(irrep)
    if matrix === nothing
        value = group_function(irrep, patterns[1], patterns[1])
        for index in 2:length(patterns)
            value += group_function(irrep, patterns[index], patterns[index])
        end
        return value
    end

    value = group_function(irrep, patterns[1], patterns[1], matrix)
    for index in 2:length(patterns)
        value += group_function(irrep, patterns[index], patterns[index], matrix)
    end
    return value
end

@testset "direct character" begin
    @testset "symbolic SU(3) and SU(4) irreps" begin
        irreps = (
            [1, 0, 0],
            [2, 0, 0],
            [1, 1, 0],
            [2, 1, 0],
            [1, 0, 0, 0],
            [2, 0, 0, 0],
            [1, 1, 0, 0],
            [2, 1, 0, 0],
        )

        for irrep in irreps
            direct = GroupFunctions.character(irrep)
            legacy = legacy_character(irrep)
            @test expand(direct - legacy) == 0
            @test direct isa Basic
            @test all(
                occursin(r"^u_\d+_\d+$", string(symbol))
                for symbol in SymEngine.free_symbols(direct)
            )
        end
    end

    @testset "explicit symbolic and numeric matrices" begin
        irrep = [2, 1, 0]
        symbolic = [SymEngine.symbols("v_$(i)_$(j)") for i in 1:3, j in 1:3]
        @test expand(GroupFunctions.character(irrep, symbolic) - legacy_character(irrep, symbolic)) == 0

        numeric = ComplexF64[
            1 + 2im  2 - 1im  0.5im
            -1 + im  3 + 0im  2 + 2im
            0.25     -2im     1 - im
        ]
        @test @inferred(GroupFunctions.character(irrep, numeric)) ≈ legacy_character(irrep, numeric)
    end

    @testset "input checks" begin
        @test_throws DimensionMismatch GroupFunctions.character([2, 1, 0], zeros(ComplexF64, 2, 2))
        @test_throws ArgumentError GroupFunctions.character([1, 2, 0])
        @test_throws ArgumentError GroupFunctions.character([1, -1])
    end
end
