"""
    CellVegf2

    Material for angiogenesis simulation, combining two species:
    endothelial cells and VEGF
    Similar to CellVegf, but this material focuses on E. cells given a VEGF distribution
    Contains:
    · DMat: Constitutive matrix of the diffusion equation of E. cells
    · VEGF: Vector containing the VEGF concentration in each node of the finite element
    · Fiber: Vector containing the fibronectin concentration in each node of the finite element
    · k_Cell: Endothelial cells diffusion constant
    · x_Vegf: Cells' chemotactic response coefficient (when VEGF concentration is ~ 0)
    · α_Vegf: Cells' saturation coefficient to chemotactic response
    · r_haptotaxis: Cells' haptotactic response coefficient

"""
mutable struct CellVegf2 <: AngiogenesisMat
    DMat :: Array{Float64,2}
    VEGF :: Vector{Float64}
    Fiber :: Vector{Float64}

    k_Cell :: Float64
    x_Vegf :: Float64
    α_Vegf :: Float64
    r_haptotaxis :: Float64
end


function CellVegf2(D::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = 1

    VEGF = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements
    Fiber = copy(VEGF)


    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    #Coefficients provided by Anderson & Chaplain, 1998
    k_Cell = 0.00035
    x_Vegf = 0.38
    α_Vegf = 0.6
    r_haptotaxis = 0.34    # Haptotaxis ON
    #r_haptotaxis = 0.       # Haptotaxis OFF


    for i in 1:Ndim
        DMat[i,i] = k_Cell
    end


    return CellVegf2(DMat, VEGF, Fiber, k_Cell, x_Vegf, α_Vegf, r_haptotaxis)
end




"""
    CellVegf2(D, k_Cell, x_Vegf)

    Constructor for CellVegf2 material.

        · D: Number of spatial dimensions
        · k_Cell: Endothelial cells diffusion constant
        · x_Vegf: Cells' chemotactic response coefficient (when VEGF concentration is ~ 0)
        · α_Vegf: Cells' saturation coefficient to chemotactic response
        · r_haptotaxis: Cells' haptotactic response coefficient
"""
function CellVegf2(D::Type{Dim}, k_Cell::Float64, x_Vegf::Float64,
        α_Vegf::Float64, r_haptotaxis::Float64) where {Dim <: AbstractDim}

    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = 1

    VEGF = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements
    Fiber = copy(VEGF)

    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    for i in 1:Ndim
        DMat[i,i] = k_Cell
    end


    return CellVegf2(DMat, VEGF, Fiber, k_Cell, x_Vegf, α_Vegf, r_haptotaxis)
end



"""
getChemotacticCoefficient

Function that calculates the chemotactic response of a CellVegf2 material χ(c)
in a specific point:

            χ(c)= χ/(1+α⋅c)

--------------------------------------------------------------------------
getChemotacticCoefficient(M, c)
        • M: CellVegf2 material containing the parameters involved
        • c: VEGF concentration vlue in the point of study
"""
function getChemotacticCoefficient(M::CellVegf2, c::Float64)
x = M.x_Vegf
α = M.α_Vegf

return x/(1+α*c)
end
