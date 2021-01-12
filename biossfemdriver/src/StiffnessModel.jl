"""
    function StiffnessModel(el, Point)

    Return the Garlekin formulation of the stiffness
    a given physical equation

    Input:

        el :: Element                            -> Element
        Point :: Union{Vector{Float64},Float64}  -> Integration point to calculate the model

    Output:

        K_i :: Array{Float64,2}  -> Local stiffnes evaluated at the integration point
"""
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim, B<: AbstractBasis, M<:AbstractMaterial}
    println("No se ha creado ningún modelo de este tipo: ")
    println("Material: $M")
    return nothing
end

LinearTypes = Union{DifusionEquation, Elastic1D, Elastic2D, Elastic3D}
"""
    function StiffnessModel(el{D,B,M}, Point) where {M <: LinearTypes}

    Return the Garlekin formulation of the stiffness
    a given physical equation
    StiffnessModel for linear Laplace equation:

        (1)       -Δu = 0

    has associated the following weak formulation of the
    local stiffness

        (2)     K_el = ∫ B^T D B dx
"""
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B<: AbstractBasis, M<:LinearTypes}

    Jac = ElementJacobian(el,Point)
    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el,J_1,Point)
    # calculate local stiffnes matrix of the Element
    K_i = Bmat'*el.Mat.DMat*Bmat*J

    return K_i
end

# StiffnessModel for a scalar convection-difusion equation
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim, B<: AbstractBasis, M<:ConvecDif}

    Jac = ElementJacobian(el,Point)

    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el, J_1, Point)
    Nmat = getNmat(el, Point)

    K_i = (Bmat'*el.Mat.DMat*Bmat + (el.Mat.Vel'*Bmat)'*Nmat')*J

    return K_i
end

# stiffnes model for difussion reaction
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: _1D,
    B <: AbstractBasis, M <: DifusionReac}

    Jac = ElementJacobian(el,Point)

    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el, J_1, Point)
    Nmat = getNmat(el, Point...)

    λ = el.Mat.chi

    K_i = (Bmat'*el.Mat.DMat*Bmat - λ*Nmat'*Nmat)*J

    return K_i
end

# StiffnessModel for cellvegf equation
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B<: AbstractBasis, M <: CellVegf2}

    Jac = ElementJacobian(el,Point)

    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el, J_1, Point)
    Nmat = getNmat(el, Point)

    VEGF = el.Mat.VEGF #vector con la concentracion de VEGF en los nodos
    c = dot(Nmat, VEGF)      #VEGF concentration in the point analyzed
    x_Vegf = getChemotacticCoefficient(el.Mat, c) #obtain chemotactic response coefficient

    Fiber =el.Mat.Fiber #vector containing fibronectin concentration in each node
    r_hapt = el.Mat.r_haptotaxis #vhaptotactic response coefficient


    K_i = (Bmat'*el.Mat.DMat*Bmat - Bmat'*Bmat*(x_Vegf*VEGF+r_hapt*Fiber)*Nmat)*J

    return K_i
end

# StiffnessModel for cellvegf equation
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B <: AbstractBasis, M <: VEGF}

    Jac = ElementJacobian(el,Point)

    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el, J_1, Point)
    Nmat = getNmat(el, Point)

    rho = Nmat*el.Mat.ECs #concentracion de ECs en el punto

    #λ = el.Mat.l_Cells
    #K_i = (Bmat'*el.Mat.DMat*Bmat - λ*Nmat'*rho*Nmat)*J

    #Modelo de comportamiento como problema difusivo de momento, modelo de reaccion va aparte
    K_i = Bmat'*el.Mat.DMat*Bmat*J
    return K_i
end

# StiffnessModel for fibronectin equation
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B <: AbstractBasis, M <: Fibronectin}

    Jac = ElementJacobian(el,Point)

    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el, J_1, Point)
    Nmat = getNmat(el, Point)

    #=
    rho = Nmat*el.Mat.ECs #concentracion de ECs en el punto

    μ = el.Mat.μ_Cells
    K_i = (Bmat'*el.Mat.DMat*Bmat - μ*Nmat'*rho*Nmat)*J
    =#

    #Modelo de comportamiento como problema difusivo de momento, modelo de reaccion va aparte
    #=SI SE UTILIZA MODELO DE REACCION DIFUSION HAY QUE INCLUIR TERMINO FORZANTE
    CORRESPONDIENTE A PRODUCCION POR ECs=#

    K_i = Bmat'*el.Mat.DMat*Bmat*J
    return K_i
end


function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B<: AbstractBasis, M<:Endothelium}

    Jac = ElementJacobian(el,Point)
    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el,J_1,Point)
    # calculate local stiffnes matrix of the Element
    K_i = Bmat'*el.Mat.DMat*Bmat*J

    return K_i
end

# StiffnessModel for TAF equation
function StiffnessModel(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B <: AbstractBasis, M <: TAF}

    Jac = ElementJacobian(el,Point)

    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el, J_1, Point)
    Nmat = getNmat(el, Point)

    rho = dot(Nmat,el.Mat.ECs) #concentracion de ECs en el punto
    hycs = dot(Nmat,el.Mat.HyCs) #Concentracion de celulas hipoxicas en el punto

    if hycs >= 0.5
        P = el.Mat.Prod_hyc
    else
        P = 0.
    end

    if rho >= 0.
        U = el.Mat.Uu_EC*rho
    else
        U = -el.Mat.Ud*rho
    end

    #Modelo de comportamiento como problema difusivo - reactivo
    #K_i = (Bmat'*el.Mat.DMat*Bmat + Nmat'*(P+U)*Nmat)*J
    K_i = (Bmat'*el.Mat.DMat*Bmat)*J
    return K_i
end
