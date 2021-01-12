#=
    This file contains routines for the elements

=#


#TODO: Reordenar este fichero


function get_el_type(el::AbstractElement{D,B,M}) where {D, B,M}
    return B
end

function getElementDim(el::AbstractElement{D,B,M}) where {D, B,M}
    return D
end

function get_el_id(el::AbstractElement{D,B,M}) where {D,B,M}
    return el.Id
end

function get_el_nodes(el::AbstractElement{D,B,M}) where {D,B,M}
    return el.Nodes
end

function get_basis_props(el::AbstractElement{D,B,M}) where {D,B,M}
    return get_basis_props(B)
end

function get_basis_gdl(el::AbstractElement{D,B,M}) where {D,B,M}
    return get_basis_gdl(B)
end

function getShape(el::AbstractElement{D,B,M}) where {D,B,M}
    return getShape(D,B)
end

function getShapeGrad(el::AbstractElement{D,B,M}) where {D,B,M}
    return getShapeGrad(D,B)
end

function get_int_points(el::AbstractElement{D,B,M})  where {D,B,M}
    return get_int_points(B,el.IntP)
end

function get_BD_int_points(el::AbstractElement{D,B,M})  where {D,B,M}
    return get_BD_int_points(B,el.IntP)
end


"""

    dict = group_by_type(x)

    Function that returns a dictionary
    containing the different elements,
    grouped by type.
"""
function group_by_type(x::Array{AbstractElement})

    # Entender la funcion MAP

    return println("Todavía no se ha implantado")
end



"""
    function SpatialNodesToIndices(Nodes::Vector{Int}, ::Type{M})

    Return the asociated degrees of freedom (DOF) for each node
    in the array Nodes

    Input:

        Nodes :: Vector{Int}    -> Array of the nodes
        Type{M}                 -> Material


    Output:

        Indices :: Vector{Int}  -> Array of the DOF Indices
"""
function SpatialNodesToIndices(Nodes::Vector{Int},
    ::Type{M}) where{M <: AbstractMaterial}

    Ndof = getMatDof[M] # get material DOF
    NofNodes = length(Nodes)

    Indices = zeros(Int,Ndof*NofNodes)

    cont = 0

    for i in 1:NofNodes
        for m in 1:Ndof
            cont += 1
            Indices[cont] = Ndof*(Nodes[i]-1) + m
        end
    end

    return Indices
end

function SpatialNodesToIndices(Nodes::Vector{Int},
    ::Type{M}, values::Union{Vector{Union{Missing, Float64}},Vector{Float64}}) where{M <: AbstractMaterial}

    Ndof = getMatDof[M] # get material DOF
    NofNodes = length(Nodes)

    Indices = zeros(Int,Ndof*NofNodes)

    cont = 0

    for i in 1:NofNodes
        for m in 1:Ndof
            cont += 1
            if values[m] !== missing
                Indices[cont] = Ndof*(Nodes[i]-1) + m
            end
        end
    end

    return Indices
end

function SpatialNodesToIndices(El::Element{D,B,M}) where{D <: AbstractDim, B <: AbstractBasis,M <: AbstractMaterial}
    return SpatialNodesToIndices(El.Nodes, M)
end
