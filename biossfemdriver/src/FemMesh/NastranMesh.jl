"""
    function NastreanMesh(::Type{D}, :: Type{B}, :: Type{M},
                        Content::NastranContent)

    This function is used to create a mesh from. It returns a
    NastranMesh structure, used to create a model.

    Input:

        D :: Type                 -> Dimension of the mesh
        B :: Type                 -> Basis of the elements shape functions
        M :: Type                 -> Material of the elements
        Content :: NastranContent -> Content from a .bdf file

    Output:

        FemMesh :: NastranMesh -> Struct that contains the mesh
"""
function NastranMesh(::Type{D}, :: Type{B}, :: Type{M},
    Content::NastranContent)  where {D <:AbstractDim, B <: AbstractBasis, M <: AbstractMaterial}

    # Get the number of Nodes
    NofNodes = length(Content.GridNodes)

    # Get the number of degrees od freedom
    DOF = NofNodes * getMatDof[M]

    # Create the elements
    Elements = CreateNastranElements(D, B, M, Content)

    return NastranMesh(Elements, NofNodes,
    DOF, Content.GridNodes, Content.GridCoords, D())
end


function CreateNastranElements(::Type{D}, :: Type{B}, :: Type{M},
    Content::NastranContent)  where {D <:AbstractDim, B <: AbstractBasis, M <: AbstractMaterial}

    Els::Array{AbstractElement,1} = []

    for i in 1:length(Content.ID)

        element_i = Element(D,B,M,Content.ID[i], Content.ElementNodes[i])

        cont = 1
        for node_i in element_i.Nodes
            index = findall(x -> x == node_i, Content.GridNodes)
            element_i.Coords[cont,:] = Content.GridCoords[:,index]
            cont += 1
        end
        push!(Els,element_i)
    end

    return Els
end
