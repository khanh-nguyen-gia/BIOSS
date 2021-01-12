"""
    CellVegf

    Material for angiogenesis simulation, combining two species:
    endothelial cells and VEGF

"""
mutable struct CellVegf <: AngiogenesisMat
    DMat :: Array
    #Velocity :: Vector{Float64}
    k_Cell :: Float64
    d_Vegf :: Float64
end


function CellVegf(D::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = getMatDof[CellVegf]

    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    k_Cell = 0.5
    d_Vegf = 3
    dif_const = [k_Cell, d_Vegf] # motility constants of the different species


    for i in 1:Ndof
            for j in 1:Ndim
                k = (i-1)*Ndim+j
                DMat[k,k] = dif_const[i]
            end
    end

    # Velocity = ones(2*Ndim)

    return CellVegf(DMat, k_Cell, d_Vegf)
end

function CellVegf(D::Type{Dim}, k_Cell::Float64, d_Vegf:: Float64) where {Dim <: AbstractDim}
    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = getMatDof[CellVegf]

    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    dif_const = [k_Cell, d_Vegf] # motility constants of the different species


    for i in 1:Ndof
            for j in 1:Ndim
                k = (i-1)*Ndim+j
                DMat[k,k] = dif_const[i]
            end
    end

    # Velocity = ones(2*Ndim)

    return CellVegf(DMat, k_Cell, d_Vegf)
end
