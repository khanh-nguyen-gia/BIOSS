struct BCForces <: AbstractBC

end


function CreateBoundary(Mesh::RegularMesh{_2D}, ::Type{M}, value::Vector{Float64},
    ::Type{BCLeft},::Type{BCForces}) where {M <: AbstractMaterial}

    ElementTable = collect(reshape(1:length(Mesh.Elements), Mesh.DimNodes...))
    Nodes = Mesh.Nodes[1,:]
    Indices = SpatialNodesToIndices(Nodes,M)

    return CreateBoundary(Mesh, Nodes, Indices, ElementTable[1,:], BCDirichlet, value)
end

function CreateBoundary(Mesh::RegularMesh{_2D}, ::Type{M}, value::Vector{Float64},
    ::Type{BCRight},::Type{BCForces}) where {M <: AbstractMaterial}

    ElementTable = collect(reshape(1:length(Mesh.Elements), Mesh.DimNodes...))
    Nodes = Mesh.Nodes[end,:]
    Indices = SpatialNodesToIndices(Nodes,M)

    return CreateBoundary(Mesh, Nodes, Indices, ElementTable[1,:], BCDirichlet, value)
end

# function CreateBoundary(Mesh::RegularMesh{D}, Nodes::Vector{Int},
#     Indices::Vector{Int}, Els::Vector{Int}, ::Type{BCForces},
#     value::Vector{Float64}) where {D <: AbstractDim}
#
#     Forces = getForces(Mesh, Nodes, Indices, Els, value)
#     return BoundaryConditions(Nodes,Indices, Els, Forces, BC())
# end
