@testset "schur_polynomial" begin
    x1 = SymEngine.symbols("x_1")
    x2 = SymEngine.symbols("x_2")
    x3 = SymEngine.symbols("x_3")

    @test schur_polynomial([0, 0]) == one(Basic)
    @test symbolic_isapprox(
        schur_polynomial([2, 0]),
        x1^2 + x1 * x2 + x2^2,
    )

    expected_210 = x1^2 * x2 + x1^2 * x3 + x1 * x2^2 + 2 * x1 * x2 * x3 +
                   x1 * x3^2 + x2^2 * x3 + x2 * x3^2

    @test symbolic_isapprox(schur_polynomial([2, 1, 0]), expected_210)
    @test schur_polynomial([2, 1, 0], [2, 3, 5]) == 280

    y = [SymEngine.symbols("y_1"), SymEngine.symbols("y_2")]
    @test symbolic_isapprox(schur_polynomial([1, 0], y), y[1] + y[2])
    @test symbolic_isapprox(schur_polynomial([3], [x1]), x1^3)
    @test symbolic_isapprox(schur_polynomial([3], variables = [x1, x2]), x1^3 + x1^2 * x2 + x1 * x2^2 + x2^3)

    @test_throws ArgumentError schur_polynomial(Int[])
    @test_throws ArgumentError schur_polynomial([1, 2])
    @test_throws ArgumentError schur_polynomial([1, -1])
    @test symbolic_isapprox(schur_polynomial([1, 0], [x1]), x1)
    @test schur_polynomial([1, 1], [x1]) == 0
end
