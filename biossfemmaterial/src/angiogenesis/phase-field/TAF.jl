"""
    TAF

    Material for TAF, tumour angiogenic factor involved in angiogenesis simulation
    Particular material for Phase field models. Contains:

        · DMat :: Array{Float64,2} -> Diffusion matrix
        · ECs :: Vector{Float64} -> EC density in the nodes
        · HyCs :: Vector{Int64} -> Info on wether there are hypoxic cells in the nodes of the element (1 = yes, 0 = no)

        · k_TAF :: Float64 -> TAF diffusion coefficient

        · TAF_hyc :: Float64  -> TAF concentration inside hypoxic cells
        · Prod_hyc :: Float64 -> TAF production by hypoxic cells
        · Uu_EC :: Float64 -> TAF uptake by endothelial cells
        · Ud_EC :: Float64 -> TAF natural decay rate
        · ratio :: Float64 -> time ratio between TAF and EC evolution, time goes faster for TAF diffusion
"""
mutable struct TAF <: AngiogenesisMat
    DMat :: Array{Float64,2}
    ECs :: Vector{Float64}
    HyCs :: Vector{Int64}

    k_TAF :: Float64

    TAF_hyc :: Float64  #TAF concentration inside hypoxic cells
    Prod_hyc :: Float64 #TAF production by hypoxic cells
    Uu_EC :: Float64 #TAF uptake by endothelial cells
    Ud :: Float64 #TAF natural decay rate
    ratio:: Float64 #time ratio between TAF and EC evolution, time goes faster for TAF diffusion
end


function TAF(D::Type{Dim}) where {Dim <: AbstractDim}
    Ndim = getSymDim(D) # get number of spatial dimensions
    Ndof = getMatDof[TAF]

    ECs = zeros(2^Ndim) #It has the size of the number of nodes of the defined
                         #element, by default square/cubic elements
    HyCs = zeros(Int64, 2^Ndim)
    DMat = zeros(Ndof*Ndim, Ndof*Ndim) # init D matrix

    k_TAF = 360.5 #μm²/h # diffusion coeficient of TAF with units

    TAF_hyc = 1.
    P = 230.77 #h⁻¹  #TAF production rate by Hypoxic cells in 1/h

    Uu = 14.4 #100 #14.4 #h⁻¹  #TAF uptake rate by endothelial cells in 1/h
    Ud = 0.7 #0.231 #h⁻¹  #TAF natural decay rate in 1/h

    # dimensionless values
    k_TAF = k_TAF*(1/L0)^2*(T0/3600)
    P = P * (T0/3600)
    Uu = Uu * (T0/3600)
    Ud = Ud * (T0/3600)

    ratio = 4000. #unused

    for i in 1:Ndim
        DMat[i,i] = k_TAF
    end


    return TAF(DMat, ECs, HyCs, k_TAF, TAF_hyc, P, Uu, Ud, ratio)
end
