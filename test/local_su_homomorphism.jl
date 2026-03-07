using GroupFunctions
import LinearAlgebra: det, qr, Diagonal, diag
import Random: seed!
using Random

function random_su(n::Int)
    x = randn(ComplexF64, n, n)
    q, r = qr(x)
    q = Matrix(q)
    d = diag(r)
    phases = d ./ abs.(d)
    q = q * Diagonal(conj.(phases))
    return q / det(q)^(1 / n)
end

function check_homomorphism(λ::Vector{Int}, n::Int, trials::Int, tol::Float64)
    for _ in 1:trials
        u = random_su(n)
        v = random_su(n)
        uv = u * v
        rep_u, _ = group_function(λ, u)
        rep_v, _ = group_function(λ, v)
        rep_uv, _ = group_function(λ, uv)
        @assert isapprox(rep_u * rep_v, rep_uv; atol=tol, rtol=tol)
    end
end

trials = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 3
tol = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1e-8
seed = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0
seed != 0 && seed!(seed)

check_homomorphism([2, 0], 2, trials, tol)
check_homomorphism([3, 0], 2, trials, tol)
check_homomorphism([2, 1, 0], 3, trials, tol)
check_homomorphism([2, 2, 0], 3, trials, tol)
check_homomorphism([2, 1, 1, 0, 0], 5, trials, tol)

println("ok")
