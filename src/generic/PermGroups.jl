@doc Markdown.doc"""
    adjacent_transpositions(g::Perm)
Return the decomposition of `g` into adjacent transpositions.
"""
function adjacent_transpositions(original::Perm)
    cycles_list = cycles(original)

    # First pass: compute exact size
    total = 0
    for cyc in cycles_list
        len = length(cyc)
        len <= 1 && continue
        @inbounds for i in 1:len-1
            a = cyc[i]; b = cyc[i+1]
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
            idx = expand_transposition!(cyc[i], cyc[i+1], output, idx)
        end
    end

    return output
end

function decompose_cycle!(cycle::AbstractVector{<:Integer}, output::Vector{NTuple{2,Int}})
    len = length(cycle)
    len <= 1 && return output

    @inbounds for i in 1:len-1
        expand_transposition!((cycle[i], cycle[i+1]), output)
    end
    return output
end

@inline function expand_transposition!(j::Int, k::Int, output::Vector{NTuple{2,Int}}, idx::Int)
    if k < j
        j, k = k, j
    end

    @inbounds for i in k:-1:j+1
        output[idx] = (i-1, i)
        idx += 1
    end
    @inbounds for i in j+1:k-1
        output[idx] = (i, i+1)
        idx += 1
    end
    return idx
end

function expand_transposition!(transposition::NTuple{2,Int}, output::Vector{NTuple{2,Int}})
    j, k = transposition
    if k < j
        j, k = k, j
    end

    needed = 2 * (k - j) - 1
    start = length(output) + 1
    resize!(output, length(output) + needed)

    idx = expand_transposition!(j, k, output, start)
    return output
end

# Backward-compatible aliases
const descomp_total = adjacent_transpositions
const descomponer_ciclo! = decompose_cycle!
const individual! = expand_transposition!
