@testset "descomposicion de una transposicion" begin
    lista::Array{Tuple{Int64,Int64},1} = [(1,1)]
    GroupFunctions.individual!((2,5), lista)
    @test lista == [(1,1), (4,5), (3,4), (2,3), (3,4), (4,5)]
end

