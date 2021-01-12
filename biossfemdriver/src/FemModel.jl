"""
    struct FemModel

    Struct that contains the Finite element model

    name :: String
        Name of the problem

    dim :: String
        Spatial dimensions of the problem

    DOF :: Int
        number of degrees of fredom

    Elements :: Array{AbstractElement,1}
        Array that contains the elements

    RestrictedNodes :: Vector{Int}
        Restricted nodes

    FreeNodes :: Vector{Int}
        Free Nodes

    Bc :: Vector{BoundaryConditions}
        Array of all the boundary conditions

    K :: Union{SparseMatCOO{Float64},SparseMatrixCSC{Float64,Int},Array{Float64,2}} # Stiffness matrix
        Stiffness Matrix

    M :: Union{SparseMatCOO{Float64},SparseMatrixCSC{Float64,Int},Array{Float64,2}} # Mass matrix
        Mass Matrix
"""
mutable struct FemModel

    name :: String     # name of the problem
    dim :: String

    DOF :: Int #number of degrees of fredom

    # elements
    Elements :: Array{AbstractElement,1}

    # Restricted nodes
    RestrictedNodes :: Vector{Int}
    # free nodes
    FreeNodes :: Vector{Int}

    # boundary Conditions
    Bc :: Array{BoundaryConditions,1}

    K :: Union{SparseMatCOO{Float64},SparseMatrixCSC{Float64,Int},Array{Float64,2}} # Stiffness matrix
    M :: Union{SparseMatCOO{Float64},SparseMatrixCSC{Float64,Int},Array{Float64,2}} # Mass matrix
end

"""
    function FemModel(name, dim, DOF, Elements,Bc;sparse=true)

    Constructor for the FemModel struct

    Input:

        name :: String                      -> Name of the problem
        dim :: String                       -> Spatial dimensions of the problem
        DOF :: Int                          -> number of degrees of fredom
        Elements :: Vector{AbstractElement} -> Array that contains the elements
        Bc :: Vector{BoundaryConditions}    -> Array of all the boundary conditions

    Optional Input:

        sparse :: Bool  -> Flag that indicates the use of sparse arrays (default true)
"""
function FemModel(name::String, dim::String, DOF::Int,
                 Elements :: Vector{AbstractElement},
                 Bc::Union{Vector{BoundaryConditions{Kind}},Vector{BoundaryConditions}};
                 sparse=true) where {Kind <: AbstractBC}


    FreeNodes = collect(1:DOF)

    # Check the restricted nodes

    RestrictedNodes = []

    for B in Bc
        if B.Kind == BCDirichlet()
            RestrictedNodes = vcat(RestrictedNodes, B.Indices)
        end
    end

    # Remove the restricted nodes from the FreeNodes

    for i in RestrictedNodes
        @debug println("Filtrando: ", i)
        filter!(x->x≠i, FreeNodes)
    end

    if sparse
        K = SparseMatCOO()
        M = SparseMatCOO()
    else
        K = zeros(Float64,DOF,DOF)
        M = zeros(Float64,DOF,DOF)
    end

    return FemModel(name,dim,DOF,Elements,RestrictedNodes,FreeNodes,Bc,K,M)
end


"""
    function FemModel(name, dim, DOF, Elements,Bc;sparse=true)

    Constructor for the FemModel struct

    Input:

        name :: String                      -> Name of the problem
        dim :: String                       -> Spatial dimensions of the problem
        FemMesh :: AbstractMesh             -> Fem Mesh of the problem
        Bc :: Vector{BoundaryConditions}    -> Array of all the boundary conditions

    Optional Input:

        sparse :: Bool  -> Flag that indicates the use of sparse arrays (default true)
"""
function FemModel(name::String, dim :: String,
    mesh::AbstractMesh{D}, Bc::Union{Vector{BoundaryConditions{Kind}},Vector{BoundaryConditions}};
    sparse=true) where {D <: AbstractDim, Kind <: AbstractBC}
    return FemModel(name, dim, mesh.DOF, mesh.Elements, Bc,sparse=sparse)
end

"""
    function Forces(model)

    Calculate the resulting forces in a model

    Input:

        model :: FemModel -> Finite element model

    Output:

        f :: Vector{Float64} -> Resulting forces
"""
function Forces(model::FemModel)
    f = spzeros(Float64, model.DOF)
    for B in model.Bc

        if B.Kind == BCDirichlet()
            for i = 1:length(B.Elements)

                local_nodes = SpatialNodesToIndices(model.Elements[B.Elements[i]])

                f[local_nodes] += LocalForces(model.Elements[B.Elements[i]], B.Cond[i])
            end
        elseif B.Kind == BCNeumann()
            for i = 1:length(B.Elements)
                local_nodes = SpatialNodesToIndices(model.Elements[B.Elements[i]])
                f[local_nodes] += calculo_fsup(
                    model.Elements[B.Elements[i]],
                    B.Cond[i])
            end
        elseif B.Kind == BCForces()

            for i = 1:length(B.Elements)

                local_nodes = SpatialNodesToIndices(model.Elements[B.Elements[i]])
                for index in 1:length(local_nodes)
                    if B.Cond[i][index] !== missing
                        f[local_nodes[index]] += B.Cond[i][index]/get_basis_gdl(model.Elements[B.Elements[i]])
                    end
                end

            end
        end
    end

    if typeof(model.Elements[1].Mat) == TAF
        f_reac = calculo_fTAF(model)
        f += f_reac

        #=
        for el in model.Elements
            local_nodes = SpatialNodesToIndices(el)
            f[local_nodes] += calculo_fTAF(el)
        end
        =#
    end
    return Array(f)
end

"""
    function uval(model)

    Calculate the restricted displacements of a Fem Model

    Input:

        model :: FemModel -> Finite element model

    Output:

        u :: Vector{Float64} -> Displacements
"""
function uval(p::FemModel)
    u = zeros(Float64, p.DOF)

    for B in p.Bc
        if B.Kind == BCDirichlet()
            for i = 1:length(B.Elements)
                local_nodes = SpatialNodesToIndices(p.Elements[B.Elements[i]])
                cont = 0
                for m in local_nodes
                    cont += 1
                    if B.Cond[i][cont] !== missing
                        u[m] = B.Cond[i][cont]
                    end
                end
            end
        end
    end
    return u
end
