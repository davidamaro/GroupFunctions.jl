using LinearAlgebra: I, mul!

export su2_block, bloquesun, su2_factorization, simplefactorization, simple, sud_from_angles
export su2_block_symbolic
export bs_block, swap_block, bs_block_symbolic, swap_block_symbolic
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
    expr_str = replace(expr_str, ", " => "_")
    expr_str = replace(expr_str, "[" => "_")
    expr_str = replace(expr_str, "]" => "")
    expr_str = replace(expr_str, r"x(_\d_\d)(\^\d{1})" => s"*(SymEngine.symbols(\"u\1\")\2)")
    expr_str = replace(expr_str, r"(x_\d_\d)"          => s"SymEngine.symbols(\"\g<1>\")")
    expr_str = replace(expr_str, " " => "")
    expr_str = replace(expr_str, "x" => "u")

    eval(Meta.parse(expr_str))
end


@doc Markdown.doc"""
    su2_block(size::Int, position::Int, angles::NTuple{3,Float64})
    su2_block(size::Int, indices::NTuple{2,Int}, angles::NTuple{3,Float64})

Embed a 2×2 SU(2) rotation defined by Euler angles `(α, β, γ)` into a
`size × size` identity matrix.

# Method 1: Adjacent indices
`su2_block(size, position, angles)` acts on rows/cols `position` and `position+1`.

# Method 2: Arbitrary indices
`su2_block(size, (i, j), angles)` acts on rows/cols `i` and `j` where `1 ≤ i < j ≤ size`.
This allows creating SU(2) blocks at non-adjacent positions.

# Examples
```julia
# Adjacent indices (traditional usage)
mat = su2_block(4, 2, (0.1, 0.2, 0.3))  # acts on rows/cols (2, 3)

# Non-adjacent indices
mat = su2_block(5, (1, 4), (0.1, 0.2, 0.3))  # acts on rows/cols (1, 4)
```
"""
function su2_block(size::Int, position::Int, angles::NTuple{3,Float64})
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

function su2_block(size::Int, indices::NTuple{2,Int}, angles::NTuple{3,Float64})
    i, j = indices
    @assert 1 <= i < j <= size "indices must satisfy 1 ≤ i < j ≤ size"

    base = Matrix{ComplexF64}(I, size, size)
    α, β, γ = angles
    cosβ = cos(β/2)
    sinβ = sin(β/2)

    # Fill the 2x2 block at positions (i,i), (i,j), (j,i), (j,j)
    base[i, i] = exp(-im*(α+γ)) * cosβ
    base[i, j] = -exp(-im*(α-γ)) * sinβ
    base[j, i] = exp(im*(α-γ)) * sinβ
    base[j, j] = exp(im*(α+γ)) * cosβ

    return base
end

@doc Markdown.doc"""
    su2_block_symbolic(size::Int, position::Int; prefix::String = "v")
    su2_block_symbolic(size::Int, indices::NTuple{2,Int}; prefix::String = "v")

Create a symbolic `size × size` matrix with a generic 2×2 block represented by
symbolic variables. The rest of the matrix is an identity.

This is the symbolic analog of `su2_block`, where instead of numeric angles, the 2×2
block elements are symbolic variables. This is useful for symbolic computations similar
to how `group_function` creates symbolic expressions in terms of matrix elements `u_i_j`.

# Method 1: Adjacent indices
`su2_block_symbolic(size, position; prefix)` creates a 2×2 block at `position` and `position+1`.

# Method 2: Arbitrary indices
`su2_block_symbolic(size, (i, j); prefix)` creates a 2×2 block at rows/cols `i` and `j`
where `1 ≤ i < j ≤ size`. This allows creating symbolic blocks at non-adjacent positions.

# Arguments
- `size::Int`: Dimension of the square matrix
- `position::Int`: Starting row/column for the 2×2 block (must satisfy `0 < position < size`)
- `indices::NTuple{2,Int}`: Tuple `(i, j)` specifying arbitrary row/column indices
- `prefix::String`: Prefix for symbolic variable names (default: "v")

# Returns
- `Matrix{Basic}`: A symbolic matrix with identity everywhere except a 2×2 symbolic block

# Examples
```julia
julia> using GroupFunctions, SymEngine

# Adjacent indices (traditional usage)
julia> mat = su2_block_symbolic(3, 1)
3×3 Matrix{Basic}:
 v_1_1  v_1_2  0
 v_2_1  v_2_2  0
     0      0  1

# Non-adjacent indices
julia> mat = su2_block_symbolic(5, (1, 4))
5×5 Matrix{Basic}:
 v_1_1  0  0  v_1_4  0
     0  1  0      0  0
     0  0  1      0  0
 v_4_1  0  0  v_4_4  0
     0  0  0      0  1

julia> mat = su2_block_symbolic(4, 2, prefix="w")
4×4 Matrix{Basic}:
 1      0      0      0
 0  w_2_2  w_2_3      0
 0  w_3_2  w_3_3      0
 0      0      0      1
```

# Notes
- The symbolic variables are named `\$(prefix)_i_j` where i,j are the row and column indices
- Unlike `su2_block`, this does not enforce SU(2) constraints (unitarity, determinant = 1)
- The resulting matrix can be used in symbolic group function computations
"""
function su2_block_symbolic(size::Int, position::Int; prefix::String = "v")
    @assert 0 < position < size "position must satisfy 0 < position < size"
    @assert position + 1 <= size "position + 1 must be ≤ size"

    # Create identity matrix with symbolic type
    base = Matrix{Basic}(I, size, size)

    # Create symbolic variables for the 2x2 block elements
    a11 = SymEngine.symbols("$(prefix)_$(position)_$(position)")
    a12 = SymEngine.symbols("$(prefix)_$(position)_$(position+1)")
    a21 = SymEngine.symbols("$(prefix)_$(position+1)_$(position)")
    a22 = SymEngine.symbols("$(prefix)_$(position+1)_$(position+1)")

    # Fill the 2x2 block
    base[position, position] = a11
    base[position, position+1] = a12
    base[position+1, position] = a21
    base[position+1, position+1] = a22

    return base
end

function su2_block_symbolic(size::Int, indices::NTuple{2,Int}; prefix::String = "v")
    i, j = indices
    @assert 1 <= i < j <= size "indices must satisfy 1 ≤ i < j ≤ size"

    # Create identity matrix with symbolic type
    base = Matrix{Basic}(I, size, size)

    # Create symbolic variables using the actual indices i, j
    a_ii = SymEngine.symbols("$(prefix)_$(i)_$(i)")
    a_ij = SymEngine.symbols("$(prefix)_$(i)_$(j)")
    a_ji = SymEngine.symbols("$(prefix)_$(j)_$(i)")
    a_jj = SymEngine.symbols("$(prefix)_$(j)_$(j)")

    # Fill the 2x2 block
    base[i, i] = a_ii
    base[i, j] = a_ij
    base[j, i] = a_ji
    base[j, j] = a_jj

    return base
end

@doc Markdown.doc"""
    bs_block(size::Int, position::Int)
    bs_block(size::Int, indices::NTuple{2,Int})

Embed a 2×2 beamsplitter matrix `[[1, -1], [1, 1]] / √2` into a `size × size` identity matrix.

This is a U(2) operation commonly used in quantum optics for modeling 50:50 beamsplitters.

# Method 1: Adjacent indices
`bs_block(size, position)` acts on rows/cols `position` and `position+1`.

# Method 2: Arbitrary indices
`bs_block(size, (i, j))` acts on rows/cols `i` and `j` where `1 ≤ i < j ≤ size`.

# Examples
```julia
# Adjacent modes
BS = bs_block(4, 2)  # acts on modes (2, 3)

# Non-adjacent modes
BS = bs_block(5, (1, 4))  # acts on modes (1, 4)
```

# Notes
- The beamsplitter matrix has determinant 1 and is unitary
- Common in quantum optics: implements a 50:50 beamsplitter
"""
function bs_block(size::Int, position::Int)
    @assert 0 < position < size "position must satisfy 0 < position < size"
    @assert position + 1 <= size "position + 1 must be ≤ size"

    base = Matrix{ComplexF64}(I, size, size)
    s = 1.0 / sqrt(2.0)

    base[position, position] = s
    base[position, position+1] = -s
    base[position+1, position] = s
    base[position+1, position+1] = s

    return base
end

function bs_block(size::Int, indices::NTuple{2,Int})
    i, j = indices
    @assert 1 <= i < j <= size "indices must satisfy 1 ≤ i < j ≤ size"

    base = Matrix{ComplexF64}(I, size, size)
    s = 1.0 / sqrt(2.0)

    base[i, i] = s
    base[i, j] = -s
    base[j, i] = s
    base[j, j] = s

    return base
end

@doc Markdown.doc"""
    bs_block_symbolic(size::Int, position::Int)
    bs_block_symbolic(size::Int, indices::NTuple{2,Int})

Embed a symbolic 2×2 beamsplitter matrix into a `size × size` identity matrix.

This creates a beamsplitter with the standard structure `[[1, -1], [1, 1]] / √2`
but using symbolic Basic types for further symbolic manipulation.

# Examples
```julia
# Adjacent indices
BS_sym = bs_block_symbolic(3, 1)

# Non-adjacent indices
BS_sym = bs_block_symbolic(5, (1, 4))
```
"""
function bs_block_symbolic(size::Int, position::Int)
    @assert 0 < position < size "position must satisfy 0 < position < size"
    @assert position + 1 <= size "position + 1 must be ≤ size"

    base = Matrix{Basic}(I, size, size)
    s = Basic(1) / sqrt(Basic(2))

    base[position, position] = s
    base[position, position+1] = -s
    base[position+1, position] = s
    base[position+1, position+1] = s

    return base
end

function bs_block_symbolic(size::Int, indices::NTuple{2,Int})
    i, j = indices
    @assert 1 <= i < j <= size "indices must satisfy 1 ≤ i < j ≤ size"

    base = Matrix{Basic}(I, size, size)
    s = Basic(1) / sqrt(Basic(2))

    base[i, i] = s
    base[i, j] = -s
    base[j, i] = s
    base[j, j] = s

    return base
end

@doc Markdown.doc"""
    swap_block(size::Int, position::Int)
    swap_block(size::Int, indices::NTuple{2,Int})

Embed a 2×2 swap matrix `[[0, 1], [1, 0]]` into a `size × size` identity matrix.

This is a U(2) operation (not SU(2), as det = -1) that swaps two modes or particles.

# Method 1: Adjacent indices
`swap_block(size, position)` acts on rows/cols `position` and `position+1`.

# Method 2: Arbitrary indices
`swap_block(size, (i, j))` acts on rows/cols `i` and `j` where `1 ≤ i < j ≤ size`.

# Examples
```julia
# Adjacent modes
SWAP = swap_block(4, 2)  # swaps modes (2, 3)

# Non-adjacent modes
SWAP = swap_block(5, (1, 4))  # swaps modes (1, 4)
```

# Notes
- The swap matrix has determinant -1 (U(2) but not SU(2))
- The swap matrix is unitary and self-inverse: SWAP² = I
"""
function swap_block(size::Int, position::Int)
    @assert 0 < position < size "position must satisfy 0 < position < size"
    @assert position + 1 <= size "position + 1 must be ≤ size"

    base = Matrix{ComplexF64}(I, size, size)

    base[position, position] = 0.0
    base[position, position+1] = 1.0
    base[position+1, position] = 1.0
    base[position+1, position+1] = 0.0

    return base
end

function swap_block(size::Int, indices::NTuple{2,Int})
    i, j = indices
    @assert 1 <= i < j <= size "indices must satisfy 1 ≤ i < j ≤ size"

    base = Matrix{ComplexF64}(I, size, size)

    base[i, i] = 0.0
    base[i, j] = 1.0
    base[j, i] = 1.0
    base[j, j] = 0.0

    return base
end

@doc Markdown.doc"""
    swap_block_symbolic(size::Int, position::Int)
    swap_block_symbolic(size::Int, indices::NTuple{2,Int})

Embed a symbolic 2×2 swap matrix into a `size × size` identity matrix.

This creates a swap with the standard structure `[[0, 1], [1, 0]]`
but using symbolic Basic types for further symbolic manipulation.

# Examples
```julia
# Adjacent indices
SWAP_sym = swap_block_symbolic(3, 1)

# Non-adjacent indices
SWAP_sym = swap_block_symbolic(5, (2, 5))
```
"""
function swap_block_symbolic(size::Int, position::Int)
    @assert 0 < position < size "position must satisfy 0 < position < size"
    @assert position + 1 <= size "position + 1 must be ≤ size"

    base = Matrix{Basic}(I, size, size)

    base[position, position] = Basic(0)
    base[position, position+1] = Basic(1)
    base[position+1, position] = Basic(1)
    base[position+1, position+1] = Basic(0)

    return base
end

function swap_block_symbolic(size::Int, indices::NTuple{2,Int})
    i, j = indices
    @assert 1 <= i < j <= size "indices must satisfy 1 ≤ i < j ≤ size"

    base = Matrix{Basic}(I, size, size)

    base[i, i] = Basic(0)
    base[i, j] = Basic(1)
    base[j, i] = Basic(1)
    base[j, j] = Basic(0)

    return base
end

"""
    su2_factorization(size::Int; skip_tail::Int = 0)

Build a random `size × size` unitary (SU(d)) as a product of SU(2) blocks; setting
`skip_tail` skips the last `skip_tail` layers of factors.
"""
function su2_factorization(size::Int; skip_tail::Int = 0)
    @assert size > skip_tail

    mat = Matrix{ComplexF64}(I, size, size)
    tmp = similar(mat)
    last_layer = size - (1 + skip_tail)
    @inbounds for layer in 1:last_layer
        @inbounds for position in layer:-1:1
            if position == 1
                angles = rand(Float64,3) .* π
                blk = su2_block(size, position, (angles[1], angles[2], angles[3]))
            else
                angles = rand(Float64,2) .* π
                blk = su2_block(size, position, (angles[1], angles[2], angles[1]))
            end
            mul!(tmp, mat, blk)
            mat, tmp = tmp, mat
        end
    end
    return mat
end

# Backward-compatible wrappers
simplefactorization(size::Int, quotient::Int) = su2_factorization(size; skip_tail = quotient)
simplefactorization(size::Int) = su2_factorization(size)

"""
    sud_from_angles(angles::Vector{Float64}, size::Int; quotient::Bool = false)

Construct an SU(d) matrix from a supplied list of Euler angles, optionally
zeroing the earliest layers when `quotient` is true.
"""
function sud_from_angles(angles::Vector{Float64}, size::Int; quotient::Bool = false)
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
                blk = su2_block(size, position, (ang1, ang2, ang3))
                mul!(tmp, blk, mat)
            else
                ang1, ang2 = angles_buf[idx], angles_buf[idx+1]
                idx += 2
                blk = su2_block(size, position, (ang1, ang2, ang1))
                mul!(tmp, blk, mat)
            end
            mat, tmp = tmp, mat
        end
    end

    return mat
end

const simple = sud_from_angles
# Backward compatibility for older callers
const bloquesun = su2_block
