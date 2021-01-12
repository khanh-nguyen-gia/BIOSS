# Convection-Difussion of one species

mutable struct DifusionReac <: ScalarMaterial
    DMat :: Array{Float64,2}  # Tensor de difusitividad
    Vel :: Vector{Float64}    # Fluid velocity

    chi :: Float64 # reaction coef

    rho :: Float64 #density
    Cv :: Float64 #heat constant
end


function DifusionReac(::Type{Dim}) where {Dim <: AbstractDim}
    DMat = zeros(getSymDim(Dim),getSymDim(Dim))
    Vel = zeros(getSymDim(Dim))

    for i in 1:getSymDim(Dim)
        DMat[i,i] = 1
        Vel[i] = 1
    end

    return DifusionReac(DMat, Vel, 1., 1., 1.)
end
