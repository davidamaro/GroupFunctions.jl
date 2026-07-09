using Test
using GroupFunctions

function reps_as_tuples(reps)
    return sort([Tuple(r.d) for r in reps])
end

function matrices_as_tuples(mats)
    return sort([Tuple(vec(m)) for m in mats])
end

function cosets_as_tuples(cosets)
    return [sort([Tuple(perm.d) for perm in coset]) for coset in cosets]
end

function brute_force_double_coset(content_a, content_b)
    reps = GroupFunctions.find_double_coset_representatives(content_a, content_b)
    left_stabilizer, right_stabilizer = GroupFunctions.stabilizer_permutations.([content_a, content_b])

    cosets = map(reps) do representative
        unique([left_perm * representative * right_perm
                for left_perm in left_stabilizer
                for right_perm in right_stabilizer])
    end

    return reps, cosets
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

@testset "double coset elements match brute force" begin
    # These contents follow the expanded convention used by group_function:
    # length(content) equals the underlying permutation size.
    cases = [
        ([1, 1], [1, 1]),
        ([2, 1, 0], [1, 2, 0]),
        ([2, 1, 1, 0], [1, 2, 1, 0]),
        ([2, 0, 2, 0], [1, 1, 2, 0]),
    ]

    for (content_a, content_b) in cases
        brute_reps, brute_cosets = brute_force_double_coset(content_a, content_b)
        reps, cosets = GroupFunctions.double_coset(content_a, content_b)

        @test reps_as_tuples(reps) == reps_as_tuples(brute_reps)
        @test cosets_as_tuples(cosets) == cosets_as_tuples(brute_cosets)
    end
end

@testset "double coset avoids stabilizer-product duplicates" begin
    _, cosets = GroupFunctions.double_coset([0, 4, 0, 0], [0, 4, 0, 0])

    @test length(cosets) == 1
    @test length(cosets[1]) == factorial(4)
    @test length(unique(cosets[1])) == factorial(4)
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
