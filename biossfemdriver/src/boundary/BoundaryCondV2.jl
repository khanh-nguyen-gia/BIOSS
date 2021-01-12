function CreateFemBC(FemMesh, ::Type{M},
    value, Nodes, ::Type{BC}) where {M <: AbstractMaterial, BC <: AbstractBC}

    Els_ID::Vector{Int} = []
    Els::Vector{Int} = []

    Elements = FemMesh.Elements

    for iel in 1:length(Elements)
        for i in Nodes

            # si i < Id del elemento, salimos (en malla regular)
            if i in Elements[iel].Nodes
                push!(Els_ID,Elements[iel].Id)
                push!(Els,iel)
                break
            end
        end
    end

    ElIndices::Vector{Vector{Int}} = []

    cont = 0

    for i_el_id in Els
        cont += 1
        push!(ElIndices,[])
        for i in Nodes
            if i in Elements[i_el_id].Nodes
                push!(ElIndices[cont],i)
            end
        end
    end

    Indices = SpatialNodesToIndices(Nodes,M)

    Forces = getForces(FemMesh, Nodes, Indices, ElIndices, Els, value, BC)

    return BoundaryConditions(Nodes, Indices, Els, Forces, BC())

end
