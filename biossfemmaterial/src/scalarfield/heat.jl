mutable struct Heat <: DifusionEquation
    DMat :: Array{Float64,2}  # Tensor de conductividad
    
    rho :: Float64 #density
    Cv :: Float64 #heat constant
end

function Heat()
    return Heat(1.,1.,1.)
end

function Heat(::Type{Dim}) where {Dim <: AbstractDim}
    K = zeros(getSymDim(Dim),getSymDim(Dim))
    for i in 1:getSymDim(Dim)
        K[i,i] = 1.
    end
    return Heat(K,1.,1.)
end

function Heat(D::Float64,::Type{Dim}) where {Dim <: AbstractDim}
    K = zeros(getSymDim(Dim),getSymDim(Dim))
    for i in 1:getSymDim(Dim)
        K[i,i] = D
    end
    return Heat(K,1.,1.)
end

function Heat(D::Float64,rho::Float64, Cv:: Float64, ::Type{Dim}) where {Dim <: AbstractDim}
    K = zeros(getSymDim(Dim),getSymDim(Dim))
    for i in 1:getSymDim(Dim)
        K[i,i] = D
    end
    return Heat(K,rho,Cv)
end
