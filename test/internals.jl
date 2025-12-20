@testset "descomposicion de una transposicion" begin
    lista::Array{Tuple{Int64,Int64},1} = [(1,1)]
    GroupFunctions.expand_transposition!((2,5), lista)
    @test lista == [(1,1), (4,5), (3,4), (2,3), (3,4), (4,5)]
end

@testset "swap_adjacent_entries" begin
    base = YoungTableau([2,1])
    fill!(base, [1,2,3])
    swapped = GroupFunctions.swap_adjacent_entries(base, 3, [2,1])
    @test swapped.fill == [1,3,2]
    @test base.fill == [1,2,3]

    single = YoungTableau([2,1])
    fill!(single, [1,1,2])
    updated = GroupFunctions.swap_adjacent_entries(single, 3, [3])
    @test updated.fill == [1,1,3]
    @test single.fill == [1,1,2]
end
