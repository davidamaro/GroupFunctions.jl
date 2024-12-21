using GroupFunctions, Test, SymEngine, Immanants
import RandomMatrices: Haar
import LinearAlgebra: norm

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
@testset "Number of zero weight states." begin
  solutions = encontrar_prototablones([1, 1, 2, 1, 0], [1, 1, 2, 1, 0])
  @test length(solutions) == 33
  solutions = encontrar_prototablones([2,2,2,0,0,0], [2,2,2,0,0,0])
  @test length(solutions) == 21
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

@testset "Comparison between Legendre polynomails and D-function" begin
  α1,β1,γ1 = rand(Float64,3)
  xx=bloquesun(3,1,(α1,β1,γ1))
  α2,β2 = rand(Float64,2)
  yy=bloquesun(3,2,(α2,β2,α2))
  α3,β3,γ3 = rand(Float64,3)
  zz=bloquesun(3,1,(α3,β3,γ3))

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

@testset "Comparisong with immanant" begin
    mat = rand(Haar(2), 3)
    pt_1 = GTPattern([[2,1,0], [2,0],[1]],[1])
    pt_2 = GTPattern([[2,1,0], [1,1],[1]],[1])
    suma::Complex{Float64} = group_function([2,1,0], pt_1, pt_1, mat)+ group_function([2,1,0], pt_2, pt_2, mat)
    @test suma ≈ immanant(Partition([2,1]), mat)
end

@testset "Comparisong with immanant" begin
    part = [2,1,1,0]
    zeroweightstates = findzero(part)

    mat = rand(Haar(2), sum(part))

    total::Complex{Float64} = sum(group_function(part, p, p, mat) for p in zeroweightstates)

    @test isapprox(abs(immanant(Partition(filter(x -> x> 0, part)), mat) -total)^2, 0.0, atol = 10^(-6))
end

@testset "Sum rules 3x3" begin
    welcome = basis_states([2,0,0]);

    α1,β1,γ1 = rand(Float64,3)
    xx=bloquesun(3,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=bloquesun(3,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=bloquesun(3,1,(α3,β3,γ3))

    mat = xx*yy*zz;
    rate1 = abs( group_function([2,1,0], welcome[5], welcome[3], mat) )^2 + abs( group_function([2,1,0], welcome[5], welcome[5], mat) )^2
    matc1 = yy*zz;
    rate2 = abs( group_function([2,1,0], welcome[5], welcome[3], matc1) )^2 + abs( group_function([2,1,0], welcome[5], welcome[5], matc1) )^2
    @test rate1 ≈ rate2
end

@testset "Sum rules 4x4" begin
    α1,β1,γ1 = rand(Float64,3)
    xx=bloquesun(4,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=bloquesun(4,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=bloquesun(4,1,(α3,β3,γ3))
    α4,β4 = rand(Float64,3)
    xx2=bloquesun(4,3,(α4,β4,α4))
    α5,β5 = rand(Float64,2)
    yy2=bloquesun(4,2,(α5,β5,α5))
    α6,β6,γ6 = rand(Float64,3)
    zz2=bloquesun(4,1,(α6,β6,γ6))

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
    xx=bloquesun(4,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=bloquesun(4,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=bloquesun(4,1,(α3,β3,γ3))
    α4,β4 = rand(Float64,3)
    xx2=bloquesun(4,3,(α4,β4,α4))
    α5,β5 = rand(Float64,2)
    yy2=bloquesun(4,2,(α5,β5,α5))
    α6,β6,γ6 = rand(Float64,3)
    zz2=bloquesun(4,1,(α6,β6,γ6))

    mat4 = xx*yy*zz*xx2*yy2*zz2
    welcome = basis_states([3,0,0,0])

    edox = filter(x -> pweight(x) == [0,1,1,1] , welcome)[1]
    edoy = filter(x -> pweight(x) == [1,0,2,0], welcome)[1]

    @test group_function([3,0,0,0], edox, edoy, mat4) ≈ immanant(Partition([3]), mat4[[1,2,3], [2,2,4]])/sqrt(2)
end

@testset "Testing labelling for construction of unitary matrices" begin
  α1,β1,γ1 = rand(Float64,3)
  xx=bloquesun(4,1,(α1,β1,γ1))
  α2,β2 = rand(Float64,2)
  yy=bloquesun(4,2,(α2,β2,α2))
  α3,β3,γ3 = rand(Float64,3)
  zz=bloquesun(4,1,(α3,β3,γ3))
  α4,β4 = rand(Float64,3)
  xx2=bloquesun(4,3,(α4,β4,α4))
  α5,β5 = rand(Float64,2)
  yy2=bloquesun(4,2,(α5,β5,α5))
  α6,β6,γ6 = rand(Float64,3)
  zz2=bloquesun(4,1,(α6,β6,γ6))

  matsimple = simple([α1,β1,γ1,α2,β2,α3,β3,γ3,α4,β4,α5,β5,α6,β6,γ6], 4)
  matsimple_quotient = simple([α1,β1,γ1,α2,β2,α3,β3,γ3,α4,β4,α5,β5,α6,β6,γ6], 4; quotient = true)
  mat =  zz2*yy2*xx2*zz*yy*xx 
  mat_quotient =  zz2*yy2*xx2#*zz*yy*xx 

  @test norm(matsimple - mat) < 10^(-11)
  @test norm(matsimple_quotient - mat_quotient) < 10^(-11)
end

@testset "Testing labelling for construction of unitary matrices" begin
    α1,β1,γ1 = rand(Float64,3)
    xx=bloquesun(3,1,(α1,β1,γ1))
    α2,β2 = rand(Float64,2)
    yy=bloquesun(3,2,(α2,β2,α2))
    α3,β3,γ3 = rand(Float64,3)
    zz=bloquesun(3,1,(α3,β3,γ3))
    mat = simple([α1,β1,γ1, α2,β2, α3,β3,γ3 ], 3)
    mat2 = zz*yy*xx;
    @test norm(mat-mat2) < 10^(-5)

end

