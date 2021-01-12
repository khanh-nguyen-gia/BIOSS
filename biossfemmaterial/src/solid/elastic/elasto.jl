using StaticArrays


mutable struct Elastic1D <: SolidMaterial
    DMat :: Array{Float64,2}   # E * A
end

function Elastic1D(::Type{Dim}) where {Dim <: AbstractDim}
    EA = 10

    D = zeros(1,1) .+ EA

    return Elastic1D(D)
end


mutable struct Elastic2D <: SolidMaterial
    DMat :: Array{Float64,2}  # Tensor elasticidad
end

function Elastic2D(::Type{Dim}) where {Dim <: AbstractDim}
    E = 10
    vu = .2
    D = zeros(3,3)

    D[1:2,1:2] = [1. vu; vu 1.]
    D[3,3] = (1. - vu)/2.
    D = D.* (E/(1. -vu*vu))
    return Elastic2D(D)
end

mutable struct Elastic3D <: SolidMaterial
    DMat ::  Union{Array{Float64,2},MArray{Tuple{6,6}}}  # Tensor elasticidad
end

function Elastic3D(::Type{Dim}) where {Dim <: AbstractDim}
    E = 10
    vu = .2
    D = @MMatrix zeros(6,6)

    D[1:3,1:3] = [(1. -vu) vu vu;
                vu (1. -vu) vu;
                vu vu (1. -vu)]
    D[4,4] =1. -2vu
    D[5,5] =1. -2vu
    D[6,6] =1. -2vu
    D = D.* (E/((1. +vu)*(1. -2*vu)))
    return Elastic3D(D)
end


mutable struct Elastic <: SolidMaterial
    DMat :: Union{Array{Float64,2},MArray{Tuple{6,6}}}  # Tensor elasticidad
end


function Elastic(::Type{_1D})
    E = 10
    vu = .2
    D = zeros(3,3)

    D[1:2,1:2] = [1. vu; vu 1.]
    D[3,3] = (1. - vu)/2.
    D = D.* (E/(1. -vu*vu))
    return Elastic(D)
end


function Elastic(::Type{_2D})
    E = 10
    vu = .2
    D = zeros(3,3)

    D[1:2,1:2] = [1. vu; vu 1.]
    D[3,3] = (1. - vu)/2.
    D = D.* (E/(1. -vu*vu))
    return Elastic(D)
end
