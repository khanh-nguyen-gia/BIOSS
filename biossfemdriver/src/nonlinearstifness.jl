# NOTE: This file is provisional, it will be upgraded in future versions

"""
    LocalStiffness(el)

    Calculate stiffness matrix of a given
    element

"""
function LocalStiffness(el::AbstractElement{D,B,M}, u0) where {D <: AbstractDim, B<: AbstractBasis, M<:TestNonElastic1D}

    N = get_basis_gdl(el)
    K = @MMatrix zeros(N*el.ndof,N*el.ndof)

    DMat = get_Dmat(el)

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points
    Npoints = size(IntPoints,1)

    # calculate element K, loop trough integration points

    for i in 1:Npoints
        # Calculate element jacobian matrix
        Jac = ElementJacobian(el,IntPoints[i,:])
        # Calculate element jacobian
        J = det(Jac)
        # Inverse of the jacobian, to calculate B matrix
        J_1 = inv(Jac)
        # Calc B matrix
        Bmat = getBmat(el,J_1,IntPoints[i,:])
        # display(Bmat)
        # display(DMat)
        Dm = DMat*(1. .-Bmat*u0) ## Material subroutine
        # Evaluate K at point of integration
        K += Bmat'*Dm*Bmat*J*Weights[i]
    end

    return K
end

T1= Vector{Union{Missing, Float64}}

function LocalForces(el::Element{D,B,M},
                   u::Union{T1,Vector{Float64}},
                   u0::Vector{Float64}) where {D <: AbstractDim, B<: AbstractBasis, M<:TestNonElastic1D}

    # get IntPoints

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N*el.ndof)

    DMat = get_Dmat(el)

    for i in 1:N*el.ndof
        for j in 1:N*el.ndof
            if u[j] !== missing
                keab = 0
                for m in 1:Npoints
                    # Calculate element jacobian matrix
                    Jac = ElementJacobian(el,IntPoints[m,:])
                    # Calculate element jacobian
                    J = det(Jac)
                    # Inverse of the jacobian, to calculate B matrix
                    J_1 = inv(Jac)
                    # Calc B matrix
                    Bmat = getBmat(el,J_1,IntPoints[m,:])
                    # Evaluate K at point of integration
                    Dm = DMat*(1 -CalcB(el,u0,IntPoints[i,:])[1])
                    keab += Bmat[:,i]'*Dm* Bmat[:,j]*J*Weights[m]
                end

                f[i] += -u[j]*keab
            end
        end
    end

    return f
end
