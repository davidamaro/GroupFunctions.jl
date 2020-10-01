export bloquesun, simplefactorization, simple

function bloquesun(tamano::Int64, pos::Int64, angulos::Tuple{Float64,Float64,Float64})
    @assert pos < tamano && pos > 0

    base = zeros(Complex{Float64}, tamano, tamano)
    α,β,γ = angulos
    base[pos,pos] = exp(-im*(α+γ))*cos(β/2)
    base[pos,pos+1] = -exp(-im*(α-γ))*sin(β/2)
    base[pos+1,pos] = exp(im*(α-γ))*sin(β/2)
    base[pos+1,pos+1] = exp(im*(α+γ))*cos(β/2)
    for i in 1:tamano
        if i == pos || i == pos + 1
            continue
        end
        base[i,i] = 1.0
    end

    return base
end

function simplefactorization(size::Int64,cociente::Int64)
    @assert size > cociente

    mat::Array{Complex{Float64},2} = zeros(Complex{Float64}, size, size)
    @simd for i in 1:size
        @inbounds mat[i,i] = 1.0
    end
    @simd for tope in 1:size-(1 + cociente)
        @simd for max in tope:-1:1
            if max == 1
                angulos = rand(Float64,3)*π
                mat *= bloquesun(size, max, (angulos[1], angulos[2], angulos[3]))
            else
                angulos = rand(Float64,2)*π
                mat *= bloquesun(size, max, (angulos[1], angulos[2], angulos[1]))
            end
        end
    end
    mat
end

function simplefactorization(size::Int64)
    mat::Array{Complex{Float64},2} = zeros(Complex{Float64}, size, size)
    @simd for i in 1:size
        @inbounds mat[i,i] = 1.0
    end
    @simd for tope in 1:size-1
        @simd for max in tope:-1:1
            if max == 1
                angulos = rand(Float64,3)*π
                mat *= bloquesun(size, max, (angulos[1], angulos[2], angulos[3]))
            else
                angulos = rand(Float64,2)*π
                mat *= bloquesun(size, max, (angulos[1], angulos[2], angulos[1]))
            end
        end
    end
    mat
end

function simple(valores::Array{Float64,1}, n::Int64; quotient::Bool = false)
    angulos = valores |> copy
    len = angulos |> length

    @assert len ≈ n^2-1

    if quotient
        angulos[1:(3*(n-2) + (n-3)*(n-2))] .= 0
        #angulos[end-(2*(n-2) +3)+1:end] .= 0
    end

    mat::Array{Complex{Float64},2} = zeros(Complex{Float64}, n, n)
    @simd for i in 1:n
        @inbounds mat[i,i] = 1.0
    end

    #angulos
    @simd for tope in 1:n-1
        @simd for max in tope:-1:1
            if max == 1
                ang = angulos[1:3]#rand(Float64,3)*π
                #mat *= bloquesun(n, max, (ang[1], ang[2], ang[3]))
                mat = bloquesun(n, max, (ang[1], ang[2], ang[3])) * mat
                angulos = angulos[4:end]
            else
                ang = angulos[1:2]
                #mat *= bloquesun(n, max, (ang[1], ang[2], ang[1]))
                mat = bloquesun(n, max, (ang[1], ang[2], ang[1])) * mat
                angulos = angulos[3:end]
            end
        end
    end

    mat
end
