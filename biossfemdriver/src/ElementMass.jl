# TODO: OPTIMIZE this calculation

"""
    function LocalMass(el)

    Calculate the local mass matrix of a given
    element

    Input:

        el :: AbstractElement{D,B,M}  -> Element to calculate mass matrix

    Output:

        M :: Array{Float64}         -> Element stiffness matrix
"""
function LocalMass(el::AbstractElement{D,B,M}) where {D <: AbstractDim,
    B<: AbstractBasis, M<:AbstractMaterial}

    N = get_basis_gdl(el)
    Mmat = zeros(N*el.ndof,N*el.ndof)
    #M = spzeros(N*el.ndof,N*el.ndof)

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points
    Npoints = size(IntPoints,1)

    # calculate element K, loop trough integration points
    for i in 1:Npoints
        # Calculate element jacobian matrix

        M_i = MassModel(el,IntPoints[i,:])

        # Evaluate K at point of integration

        Mmat += M_i*Weights[i]
    end

    return Mmat
end
