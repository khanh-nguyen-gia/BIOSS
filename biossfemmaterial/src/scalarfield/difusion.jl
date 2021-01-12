# Difussion of one species

mutable struct Difusion <: DifusionEquation
    DMat :: Array{Float64,2}  # Tensor de difusitividad
end

function Difusion()
    return Difusion(1.)
end

function Difusion(::Type{Dim}) where {Dim <: AbstractDim}
    DMat = zeros(getSymDim(Dim),getSymDim(Dim))
    for i in 1:getSymDim(Dim)
        DMat[i,i] = 1
    end
    return Difusion(DMat)
end

function Difusion(D::Float64,::Type{Dim}) where {Dim <: AbstractDim}
    DMat = zeros(getSymDim(Dim),getSymDim(Dim))
    for i in 1:getSymDim(Dim)
        DMat[i,i] = D
    end
    return Difusion(DMat)
end

mutable struct MultiDifusion{N} <: DifusionEquation
    # N : number of species
    DMat :: Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}}  # Tensor de difusitividad
end



"""
    function MultiDifusion(N,D,Dim)

    Constructor for MultiDifusion material.
        · N: Number of difusitive species
        · D: Vector containing difussitivity of each sustance
        · Dim: Number of spatial dimensions
"""
function MultiDifusion(N::Int, D::Vector{Float64},::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(Dim) # get number of spatial dimensions
    #DMat = zeros(Ndim*N,Ndim*N) # init Difusitive matrix, uncoupled
    DMat = spzeros(Ndim*N,Ndim*N)
    for i in 1:N:Ndim*N
        for j in 0:N-1
            DMat[i+j,i+j] = D[j+1]
        end
    end
    return MultiDifusion{N}(DMat)
end

"""
    function MultiDifusion(N,D,Dim)

    Constructor for MultiDifusion material.
        · N: Number of difusitive species
        · Dim: Number of spatial dimensions
"""
function MultiDifusion(N::Int,::Type{Dim}) where {Dim <: AbstractDim}
    return MultiDifusion(N,ones(N),Dim)
end

function MultiDifusion{N}(::Type{Dim}) where {N, Dim <: AbstractDim}
    getSpecies(MultiDifusion{N})
    return MultiDifusion(N,ones(N),Dim)
end

# return number of species of a MultiDifusion type
function getSpecies(::Type{MultiDifusion{N}}) where N
    return N
end

# return number of species of a MultiDifusion material
function getSpecies(Mat::MultiDifusion{N}) where N
    return getSpecies(typeof(Mat))
end
