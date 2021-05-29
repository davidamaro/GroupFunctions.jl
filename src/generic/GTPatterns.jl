#export GTPattern
#export basis_states, obtener_diferencias_patron#, prematuro_pesos#, YoungTableau
import Base.isvalid


using Base.Iterators

##############################################################################
#
#   Partition type, AbstractVector interface
#
##############################################################################


@doc Markdown.doc"""

  GTPattern(Array of arrays, Array)

# Examples:

```julia
julia> GTPattern([[2,1,0],[2,1],[2]],[2])
```
"""
mutable struct GTPattern
    filas::Array{Array{Int64,1},1}
    ultima_fila::Array{Int64,1}
end
const Fila = Array{Int64,1}

function Base.show(io::IO, ::MIME"text/plain", G::GTPattern)
    list_of_rows = G.filas
    pattern = ""
    n = length(list_of_rows)
    mitad = n ÷ 2
    between = isodd(n)
    cont = 1
    i = 1
    length_first_row = length(list_of_rows[1])*2
    for row in list_of_rows
        pattern *= "│ "
        for number in row
            pattern *= string(number)
            pattern *= " "
        end
        pattern *= repeat(" ", length_first_row - 2*length(row))
        if i <= mitad
            pattern *= repeat(" ", cont - 1)
            pattern *= "╲\n"
            cont += 1
        elseif between
            pattern *= repeat(" ", cont - 1)
            pattern *= "╳\n"
            between = false
        else
            if !between && cont > mitad
                cont -= 1
            end
            pattern *= repeat(" ", cont - 1)
            pattern *= "╱\n"
            cont -= 1
        end
        i += 1
    end

    print(io, pattern)
end

function primera!(f::Fila, a::GTPattern)
    if length(a.ultima_fila) == 0
        a.ultima_fila = f
    end
    push!(a.filas, f)
    a
end

function determinar_siguientes(f::Fila)
    siguientes = UnitRange{Int64}[]
    for i in 1:length(f)-1
        push!(siguientes, f[i+1]:f[i] )
    end
    product(siguientes...)
end

function generar_siguiente_fila(tab::GTPattern)
    siguientes = determinar_siguientes(tab.ultima_fila)

    lista_patrones = GTPattern[]
    for sig in siguientes
        tmp = deepcopy(tab)
        tmp.ultima_fila = collect(sig)
        push!(tmp.filas, collect(sig))

        push!(lista_patrones, tmp)
    end
    lista_patrones
end

function generar_siguiente_fila(tabs::Array{GTPattern,1})
    lista_patrones = GTPattern[]
    for tab in tabs
        siguientes = determinar_siguientes(tab.ultima_fila)


        for sig in siguientes
            tmp = deepcopy(tab)
            tmp.ultima_fila = collect(sig)
            push!(tmp.filas, collect(sig))

            push!(lista_patrones, tmp)
        end
    end
    lista_patrones
end

function basis_states(λ::Fila)
    ejemplo = GTPattern([], [])
    primera!(λ, ejemplo)
    multitud_prueba = generar_siguiente_fila(ejemplo);
    for _ in 3:length(λ)
        multitud_prueba = generar_siguiente_fila(multitud_prueba)
    end
    multitud_prueba
end

##############################################################################
#
#   Codigo para la traduccion
#
##############################################################################

function obtener_diferencias_patron(tab::GTPattern,numerofila::Int64)
    filas = tab.filas
    if numerofila > length(filas)
        throw(BoundsError("te pasas"))
    end
    longitud = length(filas) + 1

    diferencias = Int64[0]
    max = 0
    for fil in reverse(filas[1:longitud - numerofila])
        diferencia = fil[numerofila] - max
        if diferencia > 0
            push!(diferencias, diferencia)
            max = fil[numerofila]
        else
            push!(diferencias, 0)
        end
    end
    contenido = Int64[]
    i = numerofila
    for punto in diferencias[2:end]
        for _ in 1:punto
            push!(contenido, i)
        end
        i += 1
    end
    contenido
end

function prematuro_pesos(tab::GTPattern)
    totales = [sum(x) for x in tab.filas]
    append!(totales,0)
    final = Int[]
    for i in (1:length(tab.filas)-1)
        push!(final, totales[i]-totales[i+1])
    end
    all(x -> x == 1, final)
end


@doc Markdown.doc"""
> Custom `isvalid` for GTPattern

  isvalid(x::GTPattern)

# Examples:

```julia
julia> t = GTPattern([[2,1,0],[2,1],[2]],[2])
julia> isvalid(t)
>true


julia> t = GTPattern([[2,1,0],[2,2],[2]],[2])
julia> isvalid(t)
>false
```
"""
function isvalid(x::GTPattern)
    rows = x.filas
    for mayor in 1:length(rows)-1
        arriba = rows[mayor]
        abajo = rows[mayor+1]
        
        for (i,px) in enumerate(abajo)
            if !(px <= arriba[i] && px >= arriba[i+1])
                
                return false
            end
        end
    end
    true
end

function siguientepatron!(x::GTPattern)
    if !isnothing(disminuible(x))
        fila, col = disminuible(x)
    else
        return nothing
    end

    x.filas[fila][col] -= 1

    for j in 1:col-1
        x.filas[fila][j] = x.filas[fila-1][j]
    end


    for fil in fila+1:length(x.filas)
        for co in 1:length(x.filas[fil])
            x.filas[fil][co] = x.filas[fil-1][co]
        end
    end
    x
end

function siguientepatron(x::GTPattern)
    fila, col = disminuible(x)
    rows = deepcopy(x.filas)

    rows[fila][col] -= 1

    for j in col+1:length(rows[fila])
        rows[fila][j] = rows[fila-1][j]
    end


    for fil in fila+1:length(rows)
        for co in 1:length(rows[fil])

            rows[fil][co] = rows[fil-1][co]
        end
    end
    GTPattern(rows,rows[end])
end

function disminuible(x::GTPattern)
    rows = deepcopy(x.filas)
    for j in length(rows):-1:2
        for i in eachindex(rows[j])
            if rows[j][i] == 0
                continue
            end
            rows[j][i] -= 1
            if isvalid(GTPattern(rows, rows[end]))
                return j,i
            else
                rows[j][i] += 1
            end
        end
    end
    return nothing
end

function gtinicial(irrep::Fila)
    filas = Array{Int64,1}[irrep]
    while length(filas[end]) > 1
        push!(filas, filas[end][1:end-1])
    end
    GTPattern(filas, filas[end])
end
