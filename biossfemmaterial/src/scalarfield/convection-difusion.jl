# Convection-Difussion of one species

mutable struct ConvecDif <: ScalarMaterial
    DMat :: Array{Float64,2}  # Tensor de difusitividad
    Vel :: Vector{Float64}    # Fluid velocity

    rho :: Float64 #density
    Cv :: Float64 #heat constant
end


function ConvecDif(::Type{Dim}) where {Dim <: AbstractDim}
    DMat = zeros(getSymDim(Dim),getSymDim(Dim))
    Vel = zeros(getSymDim(Dim))

    for i in 1:getSymDim(Dim)
        DMat[i,i] = 1
        Vel[i] = 1
    end

    return ConvecDif(DMat, Vel, 1., 1.)
end

function ConvecDif(D::Float64,::Type{Dim}) where {Dim <: AbstractDim}
    DMat = zeros(getSymDim(Dim),getSymDim(Dim))
    Vel = zeros(getSymDim(Dim))
    for i in 1:getSymDim(Dim)
        DMat[i,i] = D
        Vel[i] = 1
    end
    return ConvecDif(DMat, Vel)
end


mutable struct MultiConvecDif{N} <: ScalarMaterial
    # N : number of species
    DMat :: Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}}  # Tensor de difusitividad
end



"""
    function MultiConvecDif(N,D,Dim)

    Constructor for MultiConvecDif material.
        · N: Number of difusitive species
        · D: Vector containing difussitivity of each sustance
        · Dim: Number of spatial dimensions

"""
function MultiConvecDif(N::Int, D::Vector{Float64},::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(Dim) # get number of spatial dimensions
    #DMat = zeros(Ndim*N,Ndim*N) # init Difusitive matrix, uncoupled
    DMat = spzeros(Ndim*N,Ndim*N)
    for i in 1:N:Ndim*N
        for j in 0:N-1
            DMat[i+j,i+j] = D[j+1]
        end
    end
    return MultiConvecDif{N}(DMat)
end

function MultiConvecDif(N::Int,::Type{Dim}) where {Dim <: AbstractDim}
    return MultiConvecDif(N,ones(N),Dim)
end

function MultiConvecDif{N}(::Type{Dim}) where {N, Dim <: AbstractDim}
    getSpecies(MultiConvecDif{N})
    return MultiConvecDif(N,ones(N),Dim)
end

# return number of species of a MultiConvecDif type
function getSpecies(::Type{MultiConvecDif{N}}) where N
    return N
end

# return number of species of a MultiConvecDif material
function getSpecies(Mat::MultiConvecDif{N}) where N
    return getSpecies(typeof(Mat))
end
