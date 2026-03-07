module AllSolutionsMatrix

export enumerate_matrices

"""
    sum_ints(xs)

Input:
- `xs`: vector of integers.

Output:
- `total`: sum of elements in `xs`.

Examples:
julia> sum_ints([1, 1, 3])
5
"""
function sum_ints(xs::AbstractVector{T}) where {T<:Integer}
    total = zero(T)
    @inbounds for x in xs
        total += x
    end
    return total
end

"""
    all_nonnegative(xs)

Input:
- `xs`: vector of integers.

Output:
- `ok`: `true` if every entry in `xs` is nonnegative, otherwise `false`.

Examples:
julia> all_nonnegative([0, 2, 3])
true
"""
function all_nonnegative(xs::AbstractVector{T}) where {T<:Integer}
    @inbounds for x in xs
        if x < 0
            return false
        end
    end
    return true
end

"""
    any_negative(xs)

Input:
- `xs`: vector of integers.

Output:
- `has_negative`: `true` if any entry in `xs` is negative, otherwise `false`.

Examples:
julia> any_negative([2, 0, 3])
false
"""
function any_negative(xs::AbstractVector{T}) where {T<:Integer}
    @inbounds for x in xs
        if x < 0
            return true
        end
    end
    return false
end

"""
    suffix_sums(lambda)

Input:
- `lambda`: vector of row sums.

Output:
- `suffix`: vector of length `length(lambda) + 1` where
  `suffix[i] == sum(lambda[i:end])` and `suffix[end] == 0`.

Examples:
julia> suffix_sums([1, 1, 3])
4-element Vector{Int64}:
 5
 4
 3
 0
"""
function suffix_sums(lambda::AbstractVector{T}) where {T<:Integer}
    r = length(lambda)
    suffix = Vector{T}(undef, r + 1)
    suffix[r + 1] = zero(T)
    @inbounds for i in r:-1:1
        suffix[i] = suffix[i + 1] + lambda[i]
    end
    return suffix
end

"""
    within_bounds(value, cap)

Input:
- `value`: integer to test.
- `cap`: upper bound.

Output:
- `ok`: `true` if `0 <= value <= cap`, otherwise `false`.

Examples:
julia> within_bounds(2, 3)
true
"""
function within_bounds(value::T, cap::T) where {T<:Integer}
    return value >= zero(T) && value <= cap
end

"""
    vector_subtract(capacities, row)

Input:
- `capacities`: vector of column capacities.
- `row`: vector to subtract elementwise.

Output:
- `remaining`: new vector where `remaining[j] = capacities[j] - row[j]`.

Examples:
julia> vector_subtract([3, 1], [1, 1])
2-element Vector{Int64}:
 2
 0
"""
function vector_subtract(
    capacities::AbstractVector{T},
    row::AbstractVector{T},
) where {T<:Integer}
    s = length(capacities)
    remaining = Vector{T}(undef, s)
    @inbounds for j in 1:s
        remaining[j] = capacities[j] - row[j]
    end
    return remaining
end

"""
    prepend_value(value, tail)

Input:
- `value`: integer to place first.
- `tail`: vector to append after `value`.

Output:
- `row`: new vector where `row[1] == value` and `row[2:end] == tail`.

Examples:
julia> prepend_value(2, [3, 4])
3-element Vector{Int64}:
 2
 3
 4
"""
function prepend_value(value::T, tail::Vector{T}) where {T<:Integer}
    row = Vector{T}(undef, length(tail) + 1)
    row[1] = value
    @inbounds for i in 1:length(tail)
        row[i + 1] = tail[i]
    end
    return row
end

"""
    append_row(rows, row)

Input:
- `rows`: existing list of rows.
- `row`: row to append.

Output:
- `new_rows`: new list where `row` is appended to `rows`.

Examples:
julia> append_row([[1, 0]], [0, 1])
2-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
"""
function append_row(rows::Vector{Vector{T}}, row::Vector{T}) where {T<:Integer}
    new_rows = Vector{Vector{T}}(undef, length(rows) + 1)
    @inbounds for i in 1:length(rows)
        new_rows[i] = rows[i]
    end
    new_rows[end] = row
    return new_rows
end

"""
    empty_matrices(::Type{T})

Input:
- `T`: integer element type.

Output:
- `matrices`: an empty vector of matrices (each matrix is a vector of rows).

Examples:
julia> empty_matrices(Int)
Vector{Vector{Vector{Int64}}}[]
"""
function empty_matrices(::Type{T}) where {T<:Integer}
    return Vector{Vector{Vector{T}}}()
end

"""
    empty_rows(::Type{T})

Input:
- `T`: integer element type.

Output:
- `rows`: an empty vector of row vectors.

Examples:
julia> empty_rows(Int)
Vector{Vector{Int64}}[]
"""
function empty_rows(::Type{T}) where {T<:Integer}
    return Vector{Vector{T}}()
end

"""
    margins_compatible(lambda, mu)

Input:
- `lambda`: vector of row sums.
- `mu`: vector of column sums.

Output:
- `ok`: `true` if both vectors are nonnegative and have the same total sum.

Examples:
julia> margins_compatible([2, 2], [3, 1])
true
"""
function margins_compatible(
    lambda::AbstractVector{T},
    mu::AbstractVector{T},
) where {T<:Integer}
    return all_nonnegative(lambda) &&
           all_nonnegative(mu) &&
           sum_ints(lambda) == sum_ints(mu)
end

"""
    is_pruned(capacities, suffix, i)

Input:
- `capacities`: remaining column capacities.
- `suffix`: suffix sums of the row-sum vector.
- `i`: current row index (1-based).

Output:
- `pruned`: `true` if the state violates nonnegativity or the remaining total.

Examples:
julia> is_pruned([3, 1], [4, 2, 0], 1)
false
"""
function is_pruned(
    capacities::AbstractVector{T},
    suffix::AbstractVector{T},
    i::Int,
) where {T<:Integer}
    return any_negative(capacities) || sum_ints(capacities) != suffix[i]
end

"""
    gen_bounded_compositions(target, caps)

Input:
- `target`: target row sum.
- `caps`: remaining column capacities.

Output:
- `rows`: all row vectors `row` with `sum(row) == target` and
  `0 <= row[j] <= caps[j]`.

Examples:
julia> gen_bounded_compositions(2, [3, 1])
2-element Vector{Vector{Int64}}:
 [1, 1]
 [2, 0]
"""
function gen_bounded_compositions(
    target::T,
    caps::AbstractVector{T},
) where {T<:Integer}
    s = length(caps)
    if s == 0
        return target == zero(T) ? Vector{Vector{T}}([T[]]) : empty_rows(T)
    end

    function go(col::Int, remaining::T)
        if col == s
            return within_bounds(remaining, caps[s]) ? Vector{Vector{T}}([T[remaining]]) :
                   Vector{Vector{T}}()
        end

        out = Vector{Vector{T}}()
        maxval = min(remaining, caps[col])
        if maxval < zero(T)
            return out
        end
        for val in zero(T):maxval
            tails = go(col + 1, remaining - val)
            for t in tails
                push!(out, prepend_value(val, t))
            end
        end
        return out
    end

    return go(1, target)
end

"""
    finish_last_row(partial_rows, last_row, required_sum)

Input:
- `partial_rows`: rows already fixed.
- `last_row`: candidate last row (remaining capacities).
- `required_sum`: target sum for the last row.

Output:
- `matrices`: one matrix if `sum(last_row) == required_sum`, otherwise empty.

Examples:
julia> finish_last_row([[1, 1]], [2, 0], 2)
1-element Vector{Vector{Vector{Int64}}}:
 [[1, 1], [2, 0]]
"""
function finish_last_row(
    partial_rows::Vector{Vector{T}},
    last_row::Vector{T},
    required_sum::T,
) where {T<:Integer}
    if sum_ints(last_row) == required_sum
        return Vector{Vector{Vector{T}}}([append_row(partial_rows, last_row)])
    end
    return empty_matrices(T)
end

"""
    extend_partial(i, partial_rows, capacities, lambda, suffix)

Input:
- `i`: current row index (1-based).
- `partial_rows`: list of completed rows.
- `capacities`: remaining column capacities.
- `lambda`: full row-sum vector.
- `suffix`: suffix sums of `lambda`.

Output:
- `matrices`: all matrices extending the current partial state.

Examples:
julia> extend_partial(1, Vector{Vector{Int}}(), [3, 1], [2, 2], suffix_sums([2, 2]))
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 1], [2, 0]]
 [[2, 0], [1, 1]]
"""
function extend_partial(
    i::Int,
    partial_rows::Vector{Vector{T}},
    capacities::Vector{T},
    lambda::Vector{T},
    suffix::Vector{T},
) where {T<:Integer}
    r = length(lambda)
    if is_pruned(capacities, suffix, i)
        return empty_matrices(T)
    end

    if i == r
        return finish_last_row(partial_rows, capacities, lambda[r])
    end

    row_choices = gen_bounded_compositions(lambda[i], capacities)
    out = Vector{Vector{Vector{T}}}()
    for row in row_choices
        next_capacities = vector_subtract(capacities, row)
        next_rows = append_row(partial_rows, row)
        append!(out, extend_partial(i + 1, next_rows, next_capacities, lambda, suffix))
    end
    return out
end

"""
    enumerate_matrices(lambda, mu)

Input:
- `lambda`: vector of row sums.
- `mu`: vector of column sums.

Output:
- `matrices`: all nonnegative integer matrices (as lists of rows) with row sums
  `lambda` and column sums `mu`. If the margins are incompatible, the output
  is an empty list. If input element types differ, the output uses the promoted
  integer type.

Examples:
julia> enumerate_matrices([2, 2], [3, 1])
2-element Vector{Vector{Vector{Int64}}}:
 [[1, 1], [2, 0]]
 [[2, 0], [1, 1]]
"""
function enumerate_matrices(
    lambda::AbstractVector{<:Integer},
    mu::AbstractVector{<:Integer},
)
    T = promote_type(eltype(lambda), eltype(mu))
    lambda_vec = Vector{T}(lambda)
    mu_vec = Vector{T}(mu)
    if !margins_compatible(lambda_vec, mu_vec)
        return empty_matrices(T)
    end
    suffix = suffix_sums(lambda_vec)

    return extend_partial(1, Vector{Vector{T}}(), mu_vec, lambda_vec, suffix)
end

end
