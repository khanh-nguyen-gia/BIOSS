mutable struct TestNonElastic1D <: SolidMaterial
    DMat :: Array{Float64,2}   # E * A
end

function TestNonElastic1D(::Type{Dim}) where {Dim <: AbstractDim}
    EA = 10

    D = zeros(1,1) .+ EA

    return TestNonElastic1D(D)
end
