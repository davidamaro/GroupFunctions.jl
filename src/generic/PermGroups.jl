import Base: *, ==, getindex, hash, inv, length, show

"""
    Perm(v)

Minimal permutation type used by GroupFunctions. The permutation is stored as
its image vector, so `v[i]` is the image of `i`.
"""
struct Perm{T<:Integer}
    d::Vector{T}

    function Perm{T}(v::AbstractVector{T}, check::Bool = true) where {T<:Integer}
        data = collect(v)
        if check
            n = length(data)
            sort(data) == collect(T, 1:n) || error("Unable to coerce to permutation: non-unique elements in array")
        end
        new{T}(data)
    end
end

Perm(n::T) where {T<:Integer} = Perm{T}(collect(T, 1:n), false)
Perm(v::AbstractVector{<:Integer}, check::Bool = true) = Perm{Int64}(Int64.(v), check)

length(p::Perm) = length(p.d)
getindex(p::Perm, n::Integer) = p.d[n]
==(lhs::Perm, rhs::Perm) = lhs.d == rhs.d
hash(p::Perm, h::UInt) = hash(p.d, h)

function show(io::IO, p::Perm)
    print(io, "Perm(", p.d, ")")
end

function *(g::Perm{S}, h::Perm{T}) where {S<:Integer, T<:Integer}
    length(g) == length(h) || throw(ArgumentError("incompatible permutations"))
    U = promote_type(S, T)
    Perm{U}(U.(h.d[g.d]), false)
end

function inv(g::Perm{T}) where {T<:Integer}
    inverse_data = similar(g.d)
    @inbounds for idx in eachindex(g.d)
        inverse_data[g[idx]] = idx
    end
    Perm{T}(inverse_data, false)
end

function cycles(g::Perm{T}) where {T<:Integer}
    seen = falses(length(g))
    cycle_list = Vector{Vector{T}}()
    sizehint!(cycle_list, length(g))

    for start in 1:length(g)
        seen[start] && continue
        cycle = T[]
        current = start
        while !seen[current]
            push!(cycle, current)
            seen[current] = true
            current = g[current]
        end
        push!(cycle_list, cycle)
    end

    cycle_list
end

@doc """
    adjacent_transpositions(g::Perm)
Return the decomposition of `g` into adjacent transpositions.
"""
function adjacent_transpositions(original::Perm)
    cycles_list = cycles(original)

    total = 0
    for cyc in cycles_list
        len = length(cyc)
        len <= 1 && continue
        @inbounds for i in 1:len-1
            a = cyc[i]
            b = cyc[i + 1]
            if b < a
                a, b = b, a
            end
            total += 2 * (b - a) - 1
        end
    end

    output = Vector{NTuple{2,Int}}(undef, total)
    idx = 1
    for cyc in cycles_list
        len = length(cyc)
        len <= 1 && continue
        @inbounds for i in 1:len-1
            idx = expand_transposition!(cyc[i], cyc[i + 1], output, idx)
        end
    end

    output
end

@inline function expand_transposition!(j::Int, k::Int, output::Vector{NTuple{2,Int}}, idx::Int)
    if k < j
        j, k = k, j
    end

    @inbounds for i in k:-1:j+1
        output[idx] = (i - 1, i)
        idx += 1
    end
    @inbounds for i in j+1:k-1
        output[idx] = (i, i + 1)
        idx += 1
    end
    idx
end

function expand_transposition!(transposition::NTuple{2,Int}, output::Vector{NTuple{2,Int}})
    j, k = transposition
    if k < j
        j, k = k, j
    end

    needed = 2 * (k - j) - 1
    start = length(output) + 1
    resize!(output, length(output) + needed)

    expand_transposition!(j, k, output, start)
    output
end

const descomp_total = adjacent_transpositions
const individual! = expand_transposition!
