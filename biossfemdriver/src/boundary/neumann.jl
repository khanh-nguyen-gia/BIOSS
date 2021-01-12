struct BCNeumann <: AbstractBC

end


function CreateBoundary(Mesh::RegularMesh{_2D},value::Float64,side::String,::Type{BCNeumann})


    # array that contais the element ID and
    # the boundary condition value
    C::Array{Tuple{Int,Array{Float64,1}},1} = []

    multi = getMulti(typeof(Mesh.Elements[1].Basis))

    if side == "Left"

        Nodes = zeros(Int,multi*Mesh.DimNodes[2]+1)
        # nodes that contais the bc

        for j in 1:multi*Mesh.DimNodes[2]
            Nodes[j] = (j-1)*(multi*Mesh.DimNodes[1]+1)+1
        end

        for j in 1:Mesh.DimNodes[2]
            # array con el id del elemento y con el valor de la cc
            push!(C, ((j-1)*Mesh.DimNodes[1]+1, setElBoundaryL(Mesh.Elements[(j-1)*Mesh.DimNodes[1]+1],value/2)))
        end

        Nodes[multi*Mesh.DimNodes[2]+1] = (multi*Mesh.DimNodes[2])*(multi*Mesh.DimNodes[1]+1)+1

        Bc = BoundaryConditions(Nodes, C, BCNeumann())

    elseif side  == "Right"

        Nodes = zeros(Int,multi*Mesh.DimNodes[2]+1)

        for j in 1:multi*Mesh.DimNodes[2]
            Nodes[j] = (j-1)*(multi*Mesh.DimNodes[1]+1)+1+multi*Mesh.DimNodes[1]
        end

        for j in 1:Mesh.DimNodes[2]
            push!(C,((j-1)*(Mesh.DimNodes[1])+Mesh.DimNodes[1],setElBoundaryR(Mesh.Elements[(j-1)*Mesh.DimNodes[1]+Mesh.DimNodes[1]],value/2)))
        end

        Nodes[multi*Mesh.DimNodes[2]+1] = (multi*Mesh.DimNodes[2])*(multi*Mesh.DimNodes[1]+1)+1+multi*Mesh.DimNodes[1]

        Bc = BoundaryConditions(Nodes, C, BCNeumann())

    else
        println("Introduce un valor válido")
        return None

    end

    return Bc
end
