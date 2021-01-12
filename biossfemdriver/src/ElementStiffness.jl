
# TODO: OPTIMIZE this calculation, the matrix is SYM

"""
    LocalStiffness(el)

    Calculate the local stiffness matrix of a given
    element

    Input:

        el :: AbstractElement{D,B,M}  -> Element to calculate stiffness matrix

    Output:

        K :: Array{Float64}         -> Element stiffness matrix

"""
function LocalStiffness(el::AbstractElement{D,B,M}) where {D <: AbstractDim, B<: AbstractBasis, M<:AbstractMaterial}

    N = get_basis_gdl(el)
    K = @MMatrix zeros(N*el.ndof,N*el.ndof)

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points
    Npoints = size(IntPoints,1)

    # calculate element K, loop trough integration points

    for i in 1:Npoints
        # Calculate element physical Ecuation at a given point
        Ki = StiffnessModel(el,IntPoints[i,:])
        # Evaluate K at point of integration
        K += Ki*Weights[i]
    end

    return K
end

# va a ser más eficiente calcular la K y la f a la vez
# en ocasiones

T1= Vector{Union{Missing, Float64}}

# TODO: OPTIMIZE this calculation

"""
    LocalForces(el,u)

    Calculate the local forcer of a given
    element, caused by a given displacement u

    Input:

        el :: AbstractElement{D,B,M}   -> Element to calculate stiffness matrix
        u :: Union{T1,Vector{Float64}} -> Displacement

    Output:

        f :: Vector{Float64} -> Resulting forces on the element

"""
function LocalForces(el::Element{D,B,M},
                   u::Union{T1,Vector{Float64}}) where {D <: AbstractDim,
                   B <: AbstractBasis, M <:AbstractMaterial}

    # get IntPoints

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N*el.ndof)

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
                    keab += Bmat[:,i]'*el.Mat.DMat* Bmat[:,j]*J*Weights[m]
                end

                f[i] += -u[j]*keab
            end
        end
    end

    return f
end

function LocalForces(el::Element{D,B,ConvecDif},
                   u::Union{T1,Vector{Float64}}) where {D<: AbstractDim,B<: AbstractBasis}

    # get IntPoints

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N*el.ndof)

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
                    Nmat = getNmat(el,IntPoints[i,:])
                    # Evaluate K at point of integration
                    #(el.Mat.Vel'*Bmat)'*Nmat')
                    keab += (Bmat[:,i]'*el.Mat.DMat* Bmat[:,j] + Bmat[:,j]'el.Mat.Vel*Nmat[i]')*J*Weights[m]
                    #keab += (Bmat[:,i]'*DMat* Bmat[:,j] + (el.Mat.Vel'*Bmat[:,j]*Nmat[i]'))*J*Weights[m]
                end

                f[i] += -u[j]*keab
            end
        end
    end

    return f
end

function calculo_f2(el::Element{D,B,M},
                   u::Array{Float64,1}) where {D <: AbstractDim, B<: AbstractBasis, M<:AbstractMaterial}

    # get IntPoints

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N*el.ndof)

    DMat = get_Dmat(el)

    K = spzeros(N*el.ndof,N*el.ndof)

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
        # Evaluate K at point of integration
        K += Bmat'*DMat*Bmat*J*Weights[i]
    end

    f = - K * u

    return f
end

export calculo_f2

# fuerzas superficie, validas para elementos 2D
# con parametrización regular

function calculo_fsup(el::Element{_2D,B,M},
                   h::Array{Float64,1}) where {B<: AbstractBasis, M<:AbstractMaterial}

    # get IntPoints

    (IntPoints, Weights) = get_BD_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N)

    for m in 1:Npoints

        J = norm(el.Coords[2,:]-el.Coords[3,:])
        Nmat = getNmat(el,[1,IntPoints[m]])

        f += Nmat.*h*Weights[m]*J
    end

    return f
end

# fuerzas en el volumen

function calculo_fvol(el::Element{_2D,B,M},
                   l::Array{Float64,1}) where {B<: AbstractBasis, M<:AbstractMaterial}

    # get IntPoints

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N)


    for m in 1:Npoints

        # Calculate element jacobian matrix
        Jac = ElementJacobian(el,IntPoints[m,:])
        # Calculate element jacobian
        J = det(Jac)
        Nmat = getNmat(el,[IntPoints[m,:]])

        f += Nmat.*l*Weights[m]*J
    end

    return f
end

#forces value due to reaction terms in  TAF equation
function calculo_fTAF(el::Element{D,B,M}) where {D <: AbstractDim, B<: AbstractBasis, M<:TAF}

    # get IntPoints

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N)


    for m in 1:Npoints

        # Calculate element jacobian matrix
        Jac = ElementJacobian(el,IntPoints[m,:])
        # Calculate element jacobian
        J = det(Jac)
        Nmat = getNmat(el,IntPoints[m,:])

        hycs = Nmat*el.Mat.HyCs #Concentracion de celulas hipoxicas en el punto
        TAF_hyc = el.Mat.TAF_hyc
        if hycs[1] >= 0.5
            P = el.Mat.Prod_hyc
        else
            P = 0.
        end

        f += Nmat'*P*TAF_hyc*Weights[m]*J
    end

    return f
end
