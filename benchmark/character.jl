using GroupFunctions

function legacy_character(irrep)
    patterns = basis_states(irrep)
    value = group_function(irrep, patterns[1], patterns[1])
    for index in 2:length(patterns)
        value += group_function(irrep, patterns[index], patterns[index])
    end
    return value
end

function median_sample(f; samples = 25)
    f() # compile and warm caches before measuring
    times = Vector{Float64}(undef, samples)
    bytes = Vector{Int}(undef, samples)
    for sample in 1:samples
        GC.gc()
        times[sample] = @elapsed f()
        bytes[sample] = @allocated f()
    end
    sort!(times)
    sort!(bytes)
    middle = (samples + 1) ÷ 2
    return times[middle], bytes[middle]
end

irreps = (
    [1, 0, 0],
    [2, 0, 0],
    [1, 1, 0],
    [2, 1, 0],
    [1, 0, 0, 0],
    [2, 0, 0, 0],
    [1, 1, 0, 0],
    [2, 1, 0, 0],
)

println("irrep | states | direct ms | legacy ms | speedup | direct bytes | legacy bytes")
for irrep in irreps
    direct_time, direct_bytes = median_sample(() -> character(irrep))
    legacy_time, legacy_bytes = median_sample(() -> legacy_character(irrep))
    println(
        irrep, " | ", length(basis_states(irrep)), " | ",
        round(1_000 * direct_time; digits = 3), " | ",
        round(1_000 * legacy_time; digits = 3), " | ",
        round(legacy_time / direct_time; digits = 2), "x | ",
        direct_bytes, " | ", legacy_bytes,
    )
end
