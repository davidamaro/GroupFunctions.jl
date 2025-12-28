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


@testset "expand_frequency_matrix" begin
    # Test 1: Simple 2×2 case
    proto = [2 0; 1 1]  # 2×2 matrix
    result = GroupFunctions.expand_frequency_matrix(proto)
    @test result == [1, 1, 1, 2]  # 2 ones from row 1, then 1 one and 1 two from row 2

    #Test 3/2: 3×3 case
    proto = [1 1 1; 0 2 0; 1 0 1]  # 3×3
    result = GroupFunctions.expand_frequency_matrix(proto)
    @test result == [1, 2, 3, 2, 2, 1, 3]  

    # Test 2: 1D input (should reshape)
    proto_1d = [2, 0, 1, 1]  # Will become 2×2
    result_1d = GroupFunctions.expand_frequency_matrix(proto_1d)
    @test result_1d == [1, 1, 2, 2]  
    # 2 ones (index 1), 1 two (index 3), 0 ones (index 2), and 1 two (index four)



end
