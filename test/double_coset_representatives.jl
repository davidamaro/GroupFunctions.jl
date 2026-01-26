using Test
using GroupFunctions

function reps_as_tuples(reps)
    return sort([Tuple(r.d) for r in reps])
end

@testset "double coset representatives (contents)" begin
    cases = [
        ([1, 1], [1, 1], [(1, 2), (2, 1)]),
        ([2, 1, 0], [1, 2, 0], [(1, 2, 3), (2, 3, 1)]),
        ([2, 1, 1], [1, 2, 1], [
            (1, 2, 3, 4),
            (1, 2, 4, 3),
            (1, 4, 2, 3),
            (2, 3, 1, 4),
            (2, 3, 4, 1),
            (2, 4, 1, 3),
            (2, 4, 3, 1),
        ]),
    ]

    for (A, B, expected) in cases
        reps = GroupFunctions.find_double_coset_representatives(A, B)
        @test reps_as_tuples(reps) == sort(expected)
    end
end
