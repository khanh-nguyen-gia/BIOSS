"""
    function MassModel(el, Point)

    Return the Garlekin formulation of the mass
    a given physical equation.

    Input:

        el :: Element                            -> Element
        Point :: Union{Vector{Float64},Float64}  -> Integration point to calculate the model

    Output:

        M_i :: Array{Float64,2}  -> Local mass evaluated at the integration point
"""
function MassModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B<: AbstractBasis, M<:AbstractMaterial}

    Jac = ElementJacobian(el,Point)
    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    # Calc B matrix
    Nmat = getNmat(el,Point)
    # calculate local stiffnes matrix of the Element
    M_i = Nmat'* Nmat*J

    return M_i
end

# fix temporal para elementos 1D :TODO: FIX
function MassModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: _1D,
    B<: AbstractBasis, M<:AbstractMaterial}

    Jac = ElementJacobian(el,Point)
    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    # Calc B matrix
    Nmat = getNmat(el,Point...)
    # calculate local stiffnes matrix of the Element
    M_i = Nmat'* Nmat*J

    return M_i
end

# MassModel for heat equation
function MassModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B<: AbstractBasis, M<:Union{ConvecDif,Heat}}

    Jac = ElementJacobian(el,Point)
    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    # Calc B matrix
    Nmat = getNmat(el,Point)
    # calculate local stiffnes matrix of the Element
    M_i = Nmat'* Nmat*J* el.Mat.rho*el.Mat.Cv

    return M_i
end
