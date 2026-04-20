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

    @test isapprox(immanant2110(mat), total; atol=1e-6, rtol=0)
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
