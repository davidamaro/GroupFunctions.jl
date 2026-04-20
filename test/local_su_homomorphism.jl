using GroupFunctions
import LinearAlgebra: det, qr, Diagonal, diag
import Random: seed!
using Random

function random_su(n::Int)
    gaussian_matrix = randn(ComplexF64, n, n)
    q, r = qr(gaussian_matrix)
    unitary_factor = Matrix(q)
    diagonal_entries = diag(r)
    phases = diagonal_entries ./ abs.(diagonal_entries)
    unitary_factor = unitary_factor * Diagonal(conj.(phases))
    return unitary_factor / det(unitary_factor)^(1 / n)
end

function check_homomorphism(λ::Vector{Int}, n::Int, trials::Int, tol::Float64)
    for _ in 1:trials
        left_factor = random_su(n)
        right_factor = random_su(n)
        product_factor = left_factor * right_factor
        left_representation, _ = group_function(λ, left_factor)
        right_representation, _ = group_function(λ, right_factor)
        product_representation, _ = group_function(λ, product_factor)
        @assert isapprox(left_representation * right_representation, product_representation; atol=tol, rtol=tol)
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
