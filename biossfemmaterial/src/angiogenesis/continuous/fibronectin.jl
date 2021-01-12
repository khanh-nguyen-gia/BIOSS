"""
    Fibronectin

    Material for Fibronectin, the main protein of the extracellular matrix.

"""
mutable struct Fibronectin <: AngiogenesisMat
    DMat :: Array{Float64,2}
    ECs :: Vector{Float64}


    k_Fibronectin :: Float64    #diffusion coeficient of fibronectin
    w_Cells :: Float64          # prodution rate coefficient by cells
    μ_Cells :: Float64          # uptake rate coefficient by
end


function Fibronectin(D::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = getMatDof[Fibronectin]

    ECs = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements
    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix
    k_Fibronectin = 1

    w_Cells = 0.05
    μ_Cells = 0.1



    for i in 1:Ndim
        DMat[i,i] = k_Fibronectin
    end


    return Fibronectin(DMat, ECs, k_Fibronectin, w_Cells, μ_Cells)
end





"""
    Fibronectin(D, k_Fibronectin, μ_Cells)

    Constructor for Fibronectin material.

        · D: Number of spatial dimensions
        · k_Fibronectin: Fibronectin diffusion constant
        · w_Cells : Fibronectin prodution rate coefficient by cells
        · μ_Cells: Cells' Fibronectin uptake coefficient
"""
function Fibronectin(D::Type{Dim}, k_Fibronectin::Float64, w_Cells::Float64,
                                μ_Cells::Float64) where {Dim <: AbstractDim}

    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = getMatDof[Fibronectin]

    ECs = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements

    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    for i in 1:Ndim
        DMat[i,i] = k_Fibronectin
    end


    return Fibronectin(DMat, Fibronectin, k_Fibronectin, w_Cells, μ_Cells)
end
