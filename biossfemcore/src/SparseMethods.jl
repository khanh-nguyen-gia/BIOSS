import SparseArrays: sparse, sparsevec

import Base: empty!

mutable struct SparseMatCOO{T<:Real}
    I :: Vector{Int}
    J :: Vector{Int}
    V :: Vector{T}
end

function SparseMatCOO()
    return SparseMatCOO{Float64}([], [], [])
end

function SparseMatCOO(Entries::Int)
    return SparseMatCOO{Float64}(zeros(Int,Entries),
    zeros(Int,Entries), zeros(Int,Entries))
end
# const SparseVectCOO = SparseMatCOO

function add!(A::SparseMatCOO, k::Int, I::Int, J::Int, V::Float64)
    A.I[k]=I
    A.J[k]=J
    A.V[k]=V

    return nothing
end

function add!(A::SparseMatCOO, I::Int, J::Int, V::Float64)
    push!(A.I,I)
    push!(A.J,J)
    push!(A.V,V)

    return nothing
end

function empty!(A::SparseMatCOO)
    empty!(A.I)
    empty!(A.J)
    empty!(A.V)
    return nothing
end

function append!(A::SparseMatCOO, B::SparseMatCOO)
    append!(A.I, B.I)
    append!(A.J, B.J)
    append!(A.V, B.V)
    return nothing
end

function isempty(A::SparseMatCOO)
    return isempty(A.I) && isempty(A.J) && isempty(A.V)
end

function add!(A::SparseMatCOO,index::Int, gdl1::AbstractVector{Int}, gdl2::AbstractVector{Int}, data)
    k = 1
    for j in gdl2
        for i in gdl1
            add!(A,index+k,i, j, data[k])
            k+=1
        end
    end
    return nothing
end

function add!(A::SparseMatCOO, gdl1::AbstractVector{Int}, gdl2::AbstractVector{Int}, data)
    k = 1
    for j in gdl2
        for i in gdl1
            add!(A, i, j, data[k])
            k+=1
        end
    end
    return nothing
end

SparseArrays.sparse(A::SparseMatCOO) = sparse(A.I, A.J, A.V)
Base.Matrix(A::SparseMatCOO) = Matrix(sparse(A))
Base.sizeof(A::SparseMatCOO) = sizeof(A.I) + sizeof(A.J) + sizeof(A.J)
