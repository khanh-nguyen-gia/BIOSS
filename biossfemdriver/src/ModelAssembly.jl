# TODO: Pensar en una integracion más eficiente,
# con varios procesos

# https://docs.julialang.org/en/v1/manual/parallel-computing/

"""
    Assembly!(p)

    Assembly the stiffness matrix for a given model
"""
function Assembly!(p::FemModel)
    for i in 1:length(p.Elements)
        Assembly!(p.K,p.Elements[i])
    end

    if p.Elements[1].Mat == TAF
        K_reac = Stiffnes_TAF(p)
        p.K += K_reac
    end

    return p
end

# para si sparse matrix
function Assembly!(K::SparseMatCOO{Float64},el::AbstractElement)
    nodes = SpatialNodesToIndices(el)
    add!(K,nodes,nodes,LocalStiffness(el))
end

# para no sparse
function Assembly!(K::Array{Float64,2},el::AbstractElement)
    nodes = SpatialNodesToIndices(el)
    K[nodes,nodes] += LocalStiffness(el)
end

# for regular mesh
function fastRegularAssembly!(p::FemModel)
    K = LocalStiffness(p.Elements[1])
    for i in 1:length(p.Elements)
        nodes = SpatialNodesToIndices(p.Elements[i])
        add!(p.K,nodes,nodes,K)
    end
end

export fastRegularAssembly!

# # for regular mesh
# function fastElementAssembly!(p::FemModel) # for a given element type
#     (IntPoints, Weights) = get_int_points(p.Elements[1])
#     Npoints = size(IntPoints,1)
#
#     bmatsize = getBmat(p.Elements[1],[0. 0.])
#     Bmats = zeros(Npoints,bmatsize...)
#
#     for m in 1:Npoints
#         Bmats[m,:,:] = getBmat(el,J_1,IntPoints[i,:])
#     end
#
#     for i in 1:length(p.Elements)
#         nodes = SpatialNodesToIndices(p.Elements[i])
#         add!(p.K,nodes,nodes,LocalStiffness(p.Elements[i],Bmats))
#     end
# end
#
# export fastElementAssembly!


function AssemblyMass!(p::FemModel)
    for i in 1:length(p.Elements)
        AssemblyMass!(p.M,p.Elements[i])
    end
end

# para si sparse matrix
function AssemblyMass!(M::SparseMatCOO{Float64}, el::AbstractElement)
    nodes = SpatialNodesToIndices(el)
    add!(M,nodes,nodes,LocalMass(el))
end


"""
    getKll(K,DOF)

    Return the reduced stiffness matrix at the given degrees of freedom
"""
function getKll(A::Array{Float64,2},DOF::Array{Int})
    return A[DOF,DOF]
end

function getKll(A::SparseMatCOO,DOF::Array{Int})
    return sparse(A)[DOF,DOF]
end

function getKll(A::SparseMatrixCSC{Float64,Int},DOF::Array{Int})
    return A[DOF,DOF]
end


"""
    Dismantle!(p)

    Clears the stiffness matrix of a given model
"""
function Dismantle!(p::FemModel)
    if  typeof(p.K) == SparseMatCOO{Float64}
        p.K = SparseMatCOO()
    else
        p.K = zeros(size(p.K)[1], size(p.K)[2])
    end
end

"""
    DismantleMass!(p)

    Clears the Mass matrix of a given model
"""
function DismantleMass!(p::FemModel)
    if  typeof(p.M) == SparseMatCOO{Float64}
        p.M = SparseMatCOO()
    else
        p.M = zeros(size(p.M)[1], size(p.M)[2])
    end
end

"""
    Stiffness_TAF(model::FemModel)

    Contribution to stiffness matrix of reactive terms in TAF behaviour equations
"""
function Stiffness_TAF(model::FemModel)
    K_reac = zeros(Float64, model.DOF, model.DOF)

    nodes_HyCs = BIOSSFemDriver.nodes_HyCs(model)
    nodes_ECs = BIOSSFemDriver.nodes_ECs(model)
    nodes_no_HyCs = setdiff(1:model.DOF, nodes_HyCs)
    nodes_no_ECs = setdiff(1:model.DOF, nodes_ECs)

    d_ECs = BIOSSFemDriver.extract_EC_distribution(model)

    P = model.Elements[1].Mat.Prod_hyc
    Uu = model.Elements[1].Mat.Uu_EC
    Ud = model.Elements[1].Mat.Ud

    K_reac[nodes_HyCs, nodes_HyCs] .+= P
    K_reac[nodes_ECs, nodes_ECs] += Uu*d_ECs[nodes_ECs]
    K_reac[nodes_no_ECs, nodes_no_ECs] += -Ud*d_ECs[nodes_no_ECs]

    return K_reac
end

function calculo_fTAF(model::FemModel)
    f_reac = spzeros(Float64, model.DOF)

    nodes_HyCs = BIOSSFemDriver.nodes_HyCs(model)
    nodes_no_HyCs = setdiff(1:model.DOF, nodes_HyCs)

    TAF_hyc = model.Elements[1].Mat.TAF_hyc
    P = model.Elements[1].Mat.Prod_hyc

    f_reac[nodes_HyCs] .+= P*TAF_hyc

    return f_reac
end
