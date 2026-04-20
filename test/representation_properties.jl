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

@testset "Identity irrep 100 for U(3)" begin
  mat = rand(Haar(2),3)
  kak = basis_states([1,0,0])
  actual = [group_function([1,0,0], x, y, mat) for x in kak, y in kak]
  expected = [mat[end-i+1, end-j+1] for (i, _) in enumerate(kak), (j, _) in enumerate(kak)]
  @test isapprox(actual, expected; atol=1e-10, rtol=0)
end

@testset "Identity irrep 1000 for U(4)" begin
  mat = rand(Haar(2),4)
  irrep = [1,0,0,0]
  kak = basis_states(irrep)
  actual = [group_function(irrep, x, y, mat) for x in kak, y in kak]
  expected = [mat[end-i+1, end-j+1] for (i, _) in enumerate(kak), (j, _) in enumerate(kak)]
  @test isapprox(actual, expected; atol=1e-10, rtol=0)
end

@testset "Symmetry with respect to hermitian conjugate 210" begin
  mat = rand(Haar(2),3)
  invmat = inv(mat)
  irrep = [2,1,0]
  kak = basis_states(irrep)
  direct = [group_function(irrep, x, y, mat) for x in kak, y in kak]
  conjugate_inverse = [conj(group_function(irrep, y, x, invmat)) for x in kak, y in kak]
  @test isapprox(direct, conjugate_inverse; atol=1e-10, rtol=0)
end

@testset "Symmetry with respect to hermitian conjugate 221" begin
  mat = rand(Haar(2),3)
  invmat = inv(mat)
  irrep = [2,2,1]
  kak = basis_states(irrep)
  direct = [group_function(irrep, x, y, mat) for x in kak, y in kak]
  conjugate_inverse = [conj(group_function(irrep, y, x, invmat)) for x in kak, y in kak]
  @test isapprox(direct, conjugate_inverse; atol=1e-10, rtol=0)
end
