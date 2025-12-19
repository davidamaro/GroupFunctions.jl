using LinearAlgebra: I, mul!

export bloquesun, simplefactorization, simple
export julia_to_mma, mma_to_julia

function julia_to_mma(symbolic_expression::SymEngine.Basic)
    expr = string(symbolic_expression)
    expr = replace(expr, r"u_(\d+)_(\d+)" => s"u[\1,\2]")
    expr = replace(expr, " " => "")
    expr = replace(expr, "*" => " * ")
    expr = replace(expr, r"\)([+-])" => s") \1")
    expr = replace(expr, r"([^\s])\+" => s"\1 +")
    expr = replace(expr, r"\s+" => " ")
    return strip(expr)
end



@doc Markdown.doc"""
    mma_to_julia(s::String)
> Calcula el polinomio de MMA a uno de Julia usando SymEngine

# Example:

```
julia > mma_to_julia("x[1, 1])
u_1_1
```
"""
function mma_to_julia(edo::String)
    edo = replace(edo, ", " => "_", "[" => "_", "]" => "")
    edo = replace(edo, r"x(_\d_\d)(\^\d{1})" => s"*(SymEngine.symbols(\"u\1\")\2)")
    edo = replace(edo, r"(x_\d_\d)"          => s"SymEngine.symbols(\"\g<1>\")")
    edo = replace(edo, " " => "", "x" => "u")

    eval(Meta.parse(edo))
end


@doc Markdown.doc"""
    bloquesun(size::Int64, position::Int64, angles::Tuple{Float64,Float64,Float64})
> Computes the SU(2) matrix using `angles` an embeds it into a `size` times `size` matrix.

# Example:

```julia
α5,β5 = rand(Float64,2)
yy2=bloquesun(4,2,(α5,β5,α5))
```
"""
function bloquesun(tamano::Int, pos::Int, angulos::NTuple{3,Float64})
    @assert pos < tamano && pos > 0
    @assert pos + 1 <= tamano

    base = Matrix{ComplexF64}(I, tamano, tamano)
    α,β,γ = angulos
    cosβ = cos(β/2)
    sinβ = sin(β/2)
    base[pos, pos] = exp(-im*(α+γ)) * cosβ
    base[pos, pos+1] = -exp(-im*(α-γ)) * sinβ
    base[pos+1, pos] = exp(im*(α-γ)) * sinβ
    base[pos+1, pos+1] = exp(im*(α+γ)) * cosβ
    return base
end

function simplefactorization(size::Int, cociente::Int)
    @assert size > cociente

    mat = Matrix{ComplexF64}(I, size, size)
    tmp = similar(mat)
    @inbounds for tope in 1:size-(1 + cociente)
        @inbounds for k in tope:-1:1
            if k == 1
                angulos = rand(Float64,3) .* π
                blk = bloquesun(size, k, (angulos[1], angulos[2], angulos[3]))
            else
                angulos = rand(Float64,2) .* π
                blk = bloquesun(size, k, (angulos[1], angulos[2], angulos[1]))
            end
            mul!(tmp, mat, blk)
            mat, tmp = tmp, mat
        end
    end
    return mat
end

function simplefactorization(size::Int)
    mat = Matrix{ComplexF64}(I, size, size)
    tmp = similar(mat)
    @inbounds for tope in 1:size-1
        @inbounds for k in tope:-1:1
            if k == 1
                angulos = rand(Float64,3) .* π
                blk = bloquesun(size, k, (angulos[1], angulos[2], angulos[3]))
            else
                angulos = rand(Float64,2) .* π
                blk = bloquesun(size, k, (angulos[1], angulos[2], angulos[1]))
            end
            mul!(tmp, mat, blk)
            mat, tmp = tmp, mat
        end
    end
    return mat
end

function simple(valores::Vector{Float64}, n::Int; quotient::Bool = false)
    @assert length(valores) == n^2 - 1

    angulos = copy(valores)
    if quotient
        cutoff = 3*(n-2) + (n-3)*(n-2)
        cutoff > 0 && fill!(view(angulos, 1:cutoff), 0.0)
    end

    mat = Matrix{ComplexF64}(I, n, n)
    tmp = similar(mat)
    idx = 1

    @inbounds for tope in 1:n-1
        @inbounds for k in tope:-1:1
            if k == 1
                ang1, ang2, ang3 = angulos[idx], angulos[idx+1], angulos[idx+2]
                idx += 3
                blk = bloquesun(n, k, (ang1, ang2, ang3))
                mul!(tmp, blk, mat)
            else
                ang1, ang2 = angulos[idx], angulos[idx+1]
                idx += 2
                blk = bloquesun(n, k, (ang1, ang2, ang1))
                mul!(tmp, blk, mat)
            end
            mat, tmp = tmp, mat
        end
    end

    return mat
end
