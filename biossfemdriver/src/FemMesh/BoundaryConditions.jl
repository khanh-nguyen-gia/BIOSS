abstract type AbstractBC end


struct BCRight <: AbstractBC end

struct BCLeft <: AbstractBC end

struct BCDirichlet <: AbstractBC end

struct BCForces <: AbstractBC end
struct BCNeumann <: AbstractBC end

# TODO: OPTIMIZE this

# Declaration of types

T1 = Vector{Vector{Union{Missing,Float64}}}
ValueType = Union{Vector{Float64},Vector{Union{Missing,Float64}}}


"""
    struct BoundaryConditions{K<: AbstractBC}

    Struct that contains the information of
    a given boundary condition

    Nodes :: Vector{Int}

        Nodes that are affected by a given condition

    Indices :: Vector{Int}

        Indices (DOF) that are affected by a given condition

    Elements :: Vector{Int}

        Element that are affected by a given condition


    Cond :: Union{T1,Vector{Vector{Float64}}}

        Condition that affects each element

    Kind :: K

        Kind of the boundary condition:

            BCDirichlet
            BCNeumann
            BCForces
"""
mutable struct BoundaryConditions{K<: AbstractBC}

    Nodes :: Vector{Int}

    Indices :: Vector{Int}

    Elements :: Vector{Int}

    Cond :: Union{T1,Vector{Vector{Float64}}}

    Kind :: K
end

"""
    function FindElements(Els, Nodes)

    Find all the elements that have a node in Nodes

    Input:

        Els :: Vector{AbstractElement}  -> Elements to find
        Nodes :: Vector{Int}            -> Nodes to look for

    Output:

        ElementsIndices
"""
function FindElements(Els::Vector{AbstractElement},
    Nodes::Vector{Int})

    ElementsIndices::Vector{Int} = []

    cont = 0

    ElementNodes::Vector{Vector{Int}} = []

    for iel in 1:length(Els)
        first = true
        for i in Nodes
            if i in Els[iel].Nodes
                if first
                    push!(ElementsIndices,iel)
                    cont += 1
                    push!(ElementNodes, [i])
                    first = false
                else
                    push!(ElementNodes[cont], i)
                end

            end
        end
    end

    return ElementsIndices, ElementNodes
end


"""
    function GetCondition(Elements, Value)

    Calculate the condition associated to each element
    in Elements, given by Value

    Input:

        FemMesh :: AbstractMesh
        ElementsIndex :: Vector{Int}
        Value :: Vector{Union{Missing,Float64}}

    Output:

        Condition

"""
function GetCondition(FemMesh::AbstractMesh,
    ElementsIndices::Vector{Int},
    ElementNodes::Vector{Vector{Int}}, Value::ValueType)

    # create the array that stores the boundary conditions
    # for each element
    Condition::T1 = []

    # iterate through all the elements defined in Elementsindex
    for el_index in 1:length(ElementsIndices)

        # Get the actual element
        element_i = FemMesh.Elements[ElementsIndices[el_index]]
        # init the condition for element_i
        condition_i = zeros(Union{Missing,Float64},
        length(element_i.Nodes)*element_i.ndof)

        condition_i .= missing

        # iterate all nodes indices

        for node in ElementNodes[el_index]
            internal_element_node = findall2(x -> x == node, element_i.Nodes)

            if length(internal_element_node) > 0

                internal_element_index = SpatialNodesToIndices(internal_element_node,typeof(element_i.Mat))

                for m in 1:length(Value)
                    condition_i[internal_element_index[m]] = Value[m]
                end
            end
        end

        push!(Condition,condition_i)
    end

    return Condition
end

function CreateBoundary(FemMesh::AbstractMesh{D}, ::Type{M},
    Nodes::Vector{Int}, Value::ValueType, ::Type{BC}) where {D <: AbstractDim,
    M <: AbstractMaterial, BC <: AbstractBC}

    # get the indices of the FemModel Degrees of freedom
    Indices = SpatialNodesToIndices(Nodes, M, Value)

    # get the indices of the elements that have a node in nodes
    ElementsIndices, ElementNodes = FindElements(FemMesh.Elements, Nodes)

    Condition = GetCondition(FemMesh, ElementsIndices, ElementNodes, Value)

    Kind = BC()

    return BoundaryConditions(Nodes, Indices, ElementsIndices, Condition, Kind)
end


function CreateBoundary(FemMesh::AbstractMesh{D}, ::Type{M},
    Value::ValueType, ::Type{BCLeft}, ::Type{BC}) where {D <: _2D,
    M <: AbstractMaterial, BC <: AbstractBC}

    Nodes = FemMesh.Nodes[1,:]

    return CreateBoundary(FemMesh, M, Nodes, Value, BC)
end

function CreateBoundary(FemMesh::AbstractMesh{D}, ::Type{M},
    Value::ValueType, ::Type{BCRight}, ::Type{BC}) where {D <: _2D,
    M <: AbstractMaterial, BC <: AbstractBC}

    Nodes = FemMesh.Nodes[end,:]

    return CreateBoundary(FemMesh, M, Nodes, Value, BC)
end

function CreateBoundary(FemMesh::AbstractMesh{D}, ::Type{M},
    Value::ValueType, ::Type{BCLeft}, ::Type{BC}) where {D <: _3D,
    M <: AbstractMaterial, BC <: AbstractBC}

    Nodes = vec(FemMesh.Nodes[1,:,:])

    return CreateBoundary(FemMesh, M, Nodes, Value, BC)
end

function CreateBoundary(FemMesh::AbstractMesh{D}, ::Type{M},
    Value::ValueType, ::Type{BCRight}, ::Type{BC}) where {D <: _3D,
    M <: AbstractMaterial, BC <: AbstractBC}

    Nodes = vec(FemMesh.Nodes[end,:,:])

    return CreateBoundary(FemMesh, M, Nodes, Value, BC)
end

#creamos condiciones de contorno para cuando tengamos el valor de F en todo el campo
function CreateBoundary(Mesh::AbstractMesh{D}, Fext::Array) where {D <: AbstractDim}

   Nodes = vec(Mesh.Nodes)

   Indices = SpatialNodesToIndices(Nodes,typeof(Mesh.Elements[1].Mat))
   Els = collect(reshape(1:length(Mesh.Elements),length(Mesh.Elements)))
   Forces::Vector{Vector{Float64}} = []

   for el in Els
       push!(Forces, Fext[Mesh.Elements[el].Nodes])
   end

   return BoundaryConditions(Nodes, Indices, Els, Forces, BCForces())
end
