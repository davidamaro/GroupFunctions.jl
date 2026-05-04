function _poly_equal(lhs, rhs)
    return expand(lhs - rhs) == 0
end

macro test_poly_equal(lhs, rhs)
    return :(@test _poly_equal($(esc(lhs)), $(esc(rhs))))
end

function substitute(expr, replacements::Dict)
    result = expr
    for (old, new) in replacements
        result = SymEngine.subs(result, old, new)
    end

    return expand(result)
end

@testset "test_s31_two_variables" begin
    z1, z2 = SymEngine.symbols("z1 z2")

    lhs = schur_polynomial([3, 1], variables = [z1, z2])
    rhs = z1^3 * z2 + z1^2 * z2^2 + z1 * z2^3
    @test_poly_equal lhs rhs

    normalized = schur_polynomial([3, 1, 0, 0], variables = [z1, z2])
    @test_poly_equal normalized lhs
end

@testset "test_s31_three_variables" begin
    z1, z2, z3 = SymEngine.symbols("z1 z2 z3")

    lhs = schur_polynomial([3, 1], variables = [z1, z2, z3])
    rhs = z1^3 * (z2 + z3) +
          z1^2 * (z2 + z3)^2 +
          z1 * (z2^3 + 2 * z2^2 * z3 + 2 * z2 * z3^2 + z3^3) +
          z2 * z3 * (z2^2 + z2 * z3 + z3^2)
    expanded_rhs = z1^3 * z2 + z1^3 * z3 +
                   z1^2 * z2^2 + 2 * z1^2 * z2 * z3 + z1^2 * z3^2 +
                   z1 * z2^3 + 2 * z1 * z2^2 * z3 + 2 * z1 * z2 * z3^2 + z1 * z3^3 +
                   z2^3 * z3 + z2^2 * z3^2 + z2 * z3^3

    @test_poly_equal lhs rhs
    @test_poly_equal lhs expanded_rhs
    @test_poly_equal schur_polynomial([3, 1, 0], variables = [z1, z2, z3]) lhs
end

@testset "test_s3_two_variables" begin
    z1, z2 = SymEngine.symbols("z1 z2")

    lhs = schur_polynomial([3], variables = [z1, z2])
    rhs = z1^3 + z1^2 * z2 + z1 * z2^2 + z2^3
    @test_poly_equal lhs rhs
end

@testset "test_s21_three_variables" begin
    z1, z2, z3 = SymEngine.symbols("z1 z2 z3")

    lhs = schur_polynomial([2, 1], variables = [z1, z2, z3])
    rhs = (z1 + z2) * (z1 + z3) * (z2 + z3)
    expanded_rhs = z1^2 * z2 + z1^2 * z3 +
                   z1 * z2^2 + 2 * z1 * z2 * z3 + z1 * z3^2 +
                   z2^2 * z3 + z2 * z3^2

    @test_poly_equal lhs rhs
    @test_poly_equal lhs expanded_rhs
end

@testset "test_s42_three_variables" begin
    z1, z2, z3 = SymEngine.symbols("z1 z2 z3")

    lhs = schur_polynomial([4, 2], variables = [z1, z2, z3])
    rhs = (z1^2 + z1 * z2 + z2^2) *
          (z1^2 + z1 * z3 + z3^2) *
          (z2^2 + z2 * z3 + z3^2)
    @test_poly_equal lhs rhs
end

@testset "test_s21_square_decomposition_three_variables" begin
    z1, z2, z3 = SymEngine.symbols("z1 z2 z3")
    z = [z1, z2, z3]

    lhs = schur_polynomial([2, 1], variables = z)^2
    rhs = schur_polynomial([4, 2], variables = z) +
          schur_polynomial([4, 1, 1], variables = z) +
          schur_polynomial([3, 3], variables = z) +
          2 * schur_polynomial([3, 2, 1], variables = z) +
          schur_polynomial([2, 2, 2], variables = z)
    @test_poly_equal lhs rhs
end

@testset "test_s1_three_variables" begin
    z1, z2, z3 = SymEngine.symbols("z1 z2 z3")

    lhs = schur_polynomial([1], variables = [z1, z2, z3])
    rhs = z1 + z2 + z3
    @test_poly_equal lhs rhs
end

@testset "test_s2_two_variables" begin
    x1, x2 = SymEngine.symbols("x1 x2")

    lhs = schur_polynomial([2], variables = [x1, x2])
    rhs = x1^2 + x1 * x2 + x2^2
    @test_poly_equal lhs rhs
end

@testset "test_plethysm_u3_to_u2" begin
    z1, z2, z3 = SymEngine.symbols("z1 z2 z3")
    x1, x2 = SymEngine.symbols("x1 x2")

    lhs_z = schur_polynomial([2], variables = [z1, z2, z3])
    lhs_x = substitute(lhs_z, Dict(z1 => x1^2, z2 => x2^2, z3 => x1 * x2))

    rhs_x = schur_polynomial([4], variables = [x1, x2]) +
            schur_polynomial([2, 2], variables = [x1, x2])
    expected_lhs = x1^4 + x1^3 * x2 + 2 * x1^2 * x2^2 + x1 * x2^3 + x2^4
    s4 = x1^4 + x1^3 * x2 + x1^2 * x2^2 + x1 * x2^3 + x2^4
    s22 = x1^2 * x2^2

    @test_poly_equal lhs_x rhs_x
    @test_poly_equal lhs_x expected_lhs
    @test_poly_equal rhs_x s4 + s22
end
