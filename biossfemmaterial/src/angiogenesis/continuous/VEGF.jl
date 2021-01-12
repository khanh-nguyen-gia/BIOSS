"""
    VEGF

    Material for VEGF, the growth factor involved in angiogenesis simulation

"""
mutable struct VEGF <: AngiogenesisMat
    DMat :: Array{Float64,2}
    ECs :: Vector{Float64}

    k_Vegf :: Float64
    l_Cells :: Float64  #consuption coefficient by cells
end


function VEGF(D::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = getMatDof[VEGF]

    ECs = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements
    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix
    k_Vegf = 1              #diffusion coeficient of VEGF
    l_Cells = 0.1           #uptake rate by endothelial cells


    for i in 1:Ndim
        DMat[i,i] = k_Vegf
    end


    return VEGF(DMat, ECs, k_Vegf, l_Cells)
end





"""
    VEGF(D, k_Vegf, l_Cells)

    Constructor for VEGF material.

        · D: Number of spatial dimensions
        · k_Vegf: VEGF diffusion constant
        · l_Cells: Cells' VEGF uptake coefficient
"""
function VEGF(D::Type{Dim}, k_Vegf::Float64, l_Cells::Float64) where {Dim <: AbstractDim}

    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = getMatDof[VEGF]

    ECs = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements

    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    for i in 1:Ndim
        DMat[i,i] = k_Vegf
    end


    return VEGF(DMat, VEGF, k_Vegf, l_Cells)
end
