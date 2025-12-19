using LinearAlgebra: I, mul!

export bloquesun, simplefactorization, simple
export julia_to_mma, mma_to_julia

"""
    julia_to_mma(expr::SymEngine.Basic)

Convert a SymEngine expression using `u_i_j` symbols into a Mathematica-style
string using `u[i,j]`, adding light operator spacing for readability.
"""
function julia_to_mma(expr::SymEngine.Basic)
    s = string(expr)
    s = replace(s, r"u_(\d+)_(\d+)" => s"u[\1,\2]")
    s = replace(s, " " => "")
    s = replace(s, "*" => " * ")
    s = replace(s, r"\)([+-])" => s") \1")
    s = replace(s, r"([^\s])\+" => s"\1 +")
    s = replace(s, r"\s+" => " ")
    return strip(s)
end



@doc Markdown.doc"""
    mma_to_julia(s::String)

Convert a Mathematica-style string using `x[i,j]` into a SymEngine expression
with `u_i_j` symbols.

Example:

```julia
julia> mma_to_julia("x[1, 1]")
u_1_1
```
"""
function mma_to_julia(expr_str::String)
    expr_str = replace(expr_str, ", " => "_", "[" => "_", "]" => "")
    expr_str = replace(expr_str, r"x(_\d_\d)(\^\d{1})" => s"*(SymEngine.symbols(\"u\1\")\2)")
    expr_str = replace(expr_str, r"(x_\d_\d)"          => s"SymEngine.symbols(\"\g<1>\")")
    expr_str = replace(expr_str, " " => "", "x" => "u")

    eval(Meta.parse(expr_str))
end


@doc Markdown.doc"""
    bloquesun(size::Int, position::Int, angles::NTuple{3,Float64})

Embed a 2×2 SU(2) rotation defined by Euler angles `(α, β, γ)` into an
`size × size` identity matrix, acting on rows/cols `position` and `position+1`.
"""
function bloquesun(size::Int, position::Int, angles::NTuple{3,Float64})
    @assert 0 < position < size
    @assert position + 1 <= size

    base = Matrix{ComplexF64}(I, size, size)
    α, β, γ = angles
    cosβ = cos(β/2)
    sinβ = sin(β/2)
    base[position, position] = exp(-im*(α+γ)) * cosβ
    base[position, position+1] = -exp(-im*(α-γ)) * sinβ
    base[position+1, position] = exp(im*(α-γ)) * sinβ
    base[position+1, position+1] = exp(im*(α+γ)) * cosβ
    return base
end

"""
    simplefactorization(size::Int, quotient::Int)

Build a random `size × size` unitary as a product of SU(2) blocks; the
`quotient` argument skips the last `quotient` layers of factors.
"""
function simplefactorization(size::Int, quotient::Int)
    @assert size > quotient

    mat = Matrix{ComplexF64}(I, size, size)
    tmp = similar(mat)
    @inbounds for layer in 1:size-(1 + quotient)
        @inbounds for position in layer:-1:1
            if position == 1
                angles = rand(Float64,3) .* π
                blk = bloquesun(size, position, (angles[1], angles[2], angles[3]))
            else
                angles = rand(Float64,2) .* π
                blk = bloquesun(size, position, (angles[1], angles[2], angles[1]))
            end
            mul!(tmp, mat, blk)
            mat, tmp = tmp, mat
        end
    end
    return mat
end

"""
    simplefactorization(size::Int)

Build a random `size × size` unitary as a full product of SU(2) blocks.
"""
function simplefactorization(size::Int)
    mat = Matrix{ComplexF64}(I, size, size)
    tmp = similar(mat)
    @inbounds for layer in 1:size-1
        @inbounds for position in layer:-1:1
            if position == 1
                angles = rand(Float64,3) .* π
                blk = bloquesun(size, position, (angles[1], angles[2], angles[3]))
            else
                angles = rand(Float64,2) .* π
                blk = bloquesun(size, position, (angles[1], angles[2], angles[1]))
            end
            mul!(tmp, mat, blk)
            mat, tmp = tmp, mat
        end
    end
    return mat
end

"""
    simple(angles::Vector{Float64}, size::Int; quotient::Bool = false)

Construct a unitary matrix from a supplied list of Euler angles, optionally
zeroing the earliest layers when `quotient` is true.
"""
function simple(angles::Vector{Float64}, size::Int; quotient::Bool = false)
    @assert length(angles) == size^2 - 1

    angles_buf = copy(angles)
    if quotient
        cutoff = 3*(size-2) + (size-3)*(size-2)
        cutoff > 0 && fill!(view(angles_buf, 1:cutoff), 0.0)
    end

    mat = Matrix{ComplexF64}(I, size, size)
    tmp = similar(mat)
    idx = 1

    @inbounds for layer in 1:size-1
        @inbounds for position in layer:-1:1
            if position == 1
                ang1, ang2, ang3 = angles_buf[idx], angles_buf[idx+1], angles_buf[idx+2]
                idx += 3
                blk = bloquesun(size, position, (ang1, ang2, ang3))
                mul!(tmp, blk, mat)
            else
                ang1, ang2 = angles_buf[idx], angles_buf[idx+1]
                idx += 2
                blk = bloquesun(size, position, (ang1, ang2, ang1))
                mul!(tmp, blk, mat)
            end
            mat, tmp = tmp, mat
        end
    end

    return mat
end
