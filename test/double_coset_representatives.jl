using Test
using GroupFunctions

function reps_as_tuples(reps)
    return sort([Tuple(r.d) for r in reps])
end

function matrices_as_tuples(mats)
    return sort([Tuple(vec(m)) for m in mats])
end

@testset "double coset representative matrices (contents)" begin
    # Expected matrices are the proto-tableaux/frequency matrices whose
    # row sums are content_a and column sums are content_b.
    cases = [
        ([1, 1], [1, 1], [
            [1 0; 0 1],
            [0 1; 1 0],
        ]),
        ([2, 1, 0], [1, 2, 0], [
            [0 2 0; 1 0 0; 0 0 0],
            [1 1 0; 0 1 0; 0 0 0],
        ]),
        ([2, 1, 1], [1, 2, 1], [
            [0 1 1; 1 0 0; 0 1 0],
            [0 2 0; 1 0 0; 0 0 1],
            [0 1 1; 0 1 0; 1 0 0],
            [1 0 1; 0 1 0; 0 1 0],
            [1 1 0; 0 1 0; 0 0 1],
            [0 2 0; 0 0 1; 1 0 0],
            [1 1 0; 0 0 1; 0 1 0],
        ]),
    ]

    for (content_a, content_b, expected_mats) in cases
        mats = GroupFunctions.find_double_coset_representative_matrices(content_a, content_b)
        @test matrices_as_tuples(mats) == matrices_as_tuples(expected_mats)
    end
end

@testset "double coset representatives (contents)" begin
    # Each case is (content_a, content_b, expected_reps).
    # content vectors are counts of labels 1..n (see content(::YoungTableau)),
    # so their length n fixes the permutation size. In these examples, the
    # tuples are one-line notation for permutations of 1..n.
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

    for (content_a, content_b, expected_reps) in cases
        reps = GroupFunctions.find_double_coset_representatives(content_a, content_b)
        @test reps_as_tuples(reps) == sort(expected_reps)
    end
end
