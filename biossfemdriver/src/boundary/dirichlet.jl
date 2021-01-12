abstract type AbstractBC end

struct BCRight <: AbstractBC end

struct BCLeft <: AbstractBC end

struct BCDirichlet <: AbstractBC

end

# TODO: Ordenar bien este fichero

function CreateBoundary(Mesh::RegularMesh{D}, ::Type{M}, value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    Nodes::Vector{Int}, ::Type{BC}) where {D <: AbstractDim,
    M <: AbstractMaterial, BC <: AbstractBC}

    Indices = SpatialNodesToIndices(Nodes,M,value)

    @debug "Encontrando elementos"

    Els, ElIndices = FindElements(Mesh.Elements,Nodes)

    return CreateBoundary(Mesh, Nodes, Indices, ElIndices, Els, BC,value)
end

function CreateBoundary(Mesh::RegularMesh{D}, ::Type{M}, value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    Nodes::Int, ::Type{BC}) where {D <: AbstractDim,
    M <: AbstractMaterial, BC <: AbstractBC}

    return CreateBoundary(Mesh,M,value,[Nodes],BC)
end

function CreateBoundary(Mesh::RegularMesh{_2D}, ::Type{M}, value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    ::Type{BCLeft},::Type{BCDirichlet}) where {M <: AbstractMaterial}

    ElementTable = collect(reshape(1:length(Mesh.Elements), Mesh.DimNodes...))
    Nodes = Mesh.Nodes[1,:]
    Indices = SpatialNodesToIndices(Nodes,M)

    return CreateBoundary(Mesh, Nodes, Indices, ElementTable[1,:], BCDirichlet, value)
end

function CreateBoundary(Mesh::RegularMesh{_2D}, ::Type{M}, value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    ::Type{BCRight},::Type{BCDirichlet})  where {M <: AbstractMaterial}

    ElementTable = collect(reshape(1:length(Mesh.Elements), Mesh.DimNodes...))
    Nodes = Mesh.Nodes[end,:]
    Indices = SpatialNodesToIndices(Nodes,M)

    return CreateBoundary(Mesh, Nodes, Indices, ElementTable[end,:] , BCDirichlet, value)
end

function CreateBoundary(Mesh::RegularMesh{_3D}, ::Type{M}, value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    ::Type{BCLeft},::Type{BCDirichlet}) where {M <: AbstractMaterial}

    ElementTable = collect(reshape(1:length(Mesh.Elements), Mesh.DimNodes...))

    Nodes = reshape(Mesh.Nodes[1,:,:],length(Mesh.Nodes[1,:,:]))
    Indices = SpatialNodesToIndices(Nodes,M)

    Els = reshape(ElementTable[1,:,:],length(ElementTable[1,:,:]))

    return CreateBoundary(Mesh,Nodes, Indices,Els, BCDirichlet,value)
end

function CreateBoundary(Mesh::RegularMesh{_3D}, ::Type{M}, value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    ::Type{BCRight}, ::Type{BCDirichlet}) where {M <: AbstractMaterial}

    ElementTable = collect(reshape(1:length(Mesh.Elements), Mesh.DimNodes...))

    Nodes = reshape(Mesh.Nodes[end,:,:],length(Mesh.Nodes[end,:,:]))
    Indices = SpatialNodesToIndices(Nodes,M)

    Els = reshape(ElementTable[end,:,:],length(ElementTable[end,:,:]))

    return CreateBoundary(Mesh,Nodes, Indices, Els, BCDirichlet,value)
end

function getForces(Mesh::RegularMesh{D}, Nodes::Vector{Int},
    Indices::Vector{Int}, Elements::Vector{Int},
    value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}}) where {D,M}

    Forces::Vector{Vector{Float64}} = []
    # cambiar por el.GDL para fijar más parámetros

    for el_i in Elements
        element = Mesh.Elements[el_i]
        elForces = zeros(Float64,length(element.Nodes)*element.ndof)
        for i in Nodes
            internalNode = findall(x-> x == i,element.Nodes)
            internalIndex = SpatialNodesToIndices(internalNode,typeof(element.Mat))

            # println(internalIndex)
            # println(value)
            # TODO: FIX this!!
            if internalIndex != Int64[]
                for m in 1:length(value)
                    if value[m] != nothing
                        elForces[internalIndex[m]] += value[m]
                    end
                end
            end
        end
        push!(Forces,elForces)
    end

    return Forces
end

function getForces(Mesh::RegularMesh{D}, Nodes::Vector{Int},
    Indices::Vector{Int},ElIndices::Vector{Vector{Int}}, Elements::Vector{Int},
    value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    ::BC) where {D,M,BC}

    Forces::Vector{Vector{Float64}} = []
    # cambiar por el.GDL para fijar más parámetros

    Nels = length(Elements)

    for el_i in 1:Nels
        element = Mesh.Elements[el_i]
        elForces = zeros(Float64,length(element.Nodes)*element.ndof)
        for i in ElIndices[el_i]
            internalNode = findall(x-> x == i,element.Nodes)
            internalIndex = SpatialNodesToIndices(internalNode,typeof(element.Mat))
            # println(internalIndex)
            # println(value)
            # TODO: FIX this!!
            if internalIndex != Int64[]

                for m in 1:length(value)
                    elForces[internalIndex[m]] = 1.
                    if value[m] != nothing
                        elForces[internalIndex[m]] += value[m]
                    end
                end
            end
        end
        push!(Forces,elForces)
    end

    return Forces
end

function getForces(Mesh::RegularMesh{D}, Nodes::Vector{Int},
    Indices::Vector{Int},ElIndices::Vector{Vector{Int}}, Elements::Vector{Int},
    value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}},
    ::Type{BCDirichlet}) where {D,M}

    Forces::T1 = []
    # cambiar por el.GDL para fijar más parámetros

    Nels = length(Elements)

    for el_i in 1:Nels
        # solucion provisional ...
        element = Mesh.Elements[Elements[el_i]]

        elForces = zeros(Union{Missing,Float64},length(element.Nodes)*element.ndof)
        elForces .= missing
        for i in ElIndices[el_i]
            internalNode = findall(x-> x == i,element.Nodes)
            internalIndex = SpatialNodesToIndices(internalNode,typeof(element.Mat))
            # println(internalIndex)
            # println(value)
            # TODO: FIX this!!
            if internalIndex != Int64[]
                for m in 1:length(value)
                    elForces[internalIndex[m]] = 1.
                    if value[m] != nothing
                        elForces[internalIndex[m]] = value[m]
                    end
                end
            end
        end
        push!(Forces,elForces)
    end

    return Forces
end

function CreateBoundary(Mesh::RegularMesh{D}, Nodes::Vector{Int},
    Indices::Vector{Int}, Els::Vector{Int}, ::Type{BC},
    value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}}) where {D <: AbstractDim, BC <: AbstractBC}
    Forces = getForces(Mesh, Nodes, Indices, Els, value)
    return BoundaryConditions(Nodes,Indices, Els, Forces, BC())
end

function CreateBoundary(Mesh::RegularMesh{D}, Nodes::Vector{Int},
    Indices::Vector{Int}, ElIndices::Vector{Vector{Int}}, Els::Vector{Int}, ::Type{BC},
    value::Union{Vector{Union{Nothing, Float64}},Vector{Float64}}) where {D <: AbstractDim, BC <: AbstractBC}
    Forces = getForces(Mesh, Nodes, Indices,ElIndices, Els, value, BC)
    return BoundaryConditions(Nodes,Indices,Els, Forces, BC())
end

# creamos condiciones de contorno para cuando tengamos el valor de F en todo el campo
function CreateBoundary(Mesh::RegularMesh{D}, Fext::Array) where {D <: AbstractDim}
    Nodes = vec(Mesh.Nodes)

    Indices = SpatialNodesToIndices(Nodes,typeof(Mesh.Elements[1].Mat))
    Els = collect(reshape(1:length(Mesh.Elements),length(Mesh.Elements)))
    Forces::Vector{Vector{Float64}} = []

    for el in Els
        push!(Forces, Fext[Mesh.Elements[el].Nodes])
    end

    return BoundaryConditions(Nodes, Indices, Els, Forces, BCForces())
end
