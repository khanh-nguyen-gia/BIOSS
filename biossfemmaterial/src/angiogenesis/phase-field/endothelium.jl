"""
    Endothelium

    Material for angiogenesis simulation using a phase -field model.
    This material simulates endothelial tissue using continuous equations
    Contains:
    · DMat: Constitutive matrix of the diffusion equation of E. cells
    · VEGF: Vector containing the VEGF concentration in each node of the finite element
    · M_Cell: Endothelial cells mobility (diffusion) constant
    · λ_Cell: Endothelial cells capilary wall width constant
    · α_Cell: Cells phenotype switch constant 1
    · β_Cell: Cells phenotype switch constant 2
    · f_act: Minimum TAF concentration required for Tip Cell activation


"""
mutable struct Endothelium <: AngiogenesisMat
    DMat :: Array{Float64,2}
    VEGF :: Vector{Float64}

    M_Cell :: Float64
    λ_Cell :: Float64
    α_Cell :: Float64
    β_Cell :: Float64
    f_act :: Float64
end


function Endothelium(D::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = 1

    VEGF = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements

    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    #Coefficients provided by Anderson & Chaplain, 1998
    M_Cell = 0.0577 #h⁻¹: Mobility coefficient with units #0.0875
    λ_Cell = 1.25 #μm: interface width with units
    α_Cell = 0.525
    β_Cell = 10. ^4
    #f_act = 0.001
    f_act = 0.3

    #Nondimensionalized values
    M_Cell *= (T0/3600)
    λ_Cell /= L0


    for i in 1:Ndim
        DMat[i,i] = M_Cell*(λ_Cell^2)
    end


    return Endothelium(DMat, VEGF, M_Cell, λ_Cell, α_Cell, β_Cell, f_act)
end
