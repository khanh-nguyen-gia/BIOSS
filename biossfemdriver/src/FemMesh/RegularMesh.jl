
"""
    function RegularMesh(::Type{D}, :: Type{B}, :: Type{M},
                         DimNodes :: Array{Int,1},
                         EndPoint :: Array{Float64,1})

    This function is used to create a parametric mesh. It returns a
    RegularMesh structure, used to create a model.

    Input:

        D :: Type               -> Dimension of the mesh
        B :: Type               -> Basis of the elements shape functions
        M :: Type               -> Material of the elements
        DimNodes :: Vector{Int} -> Number of elements in each dimension
        EndPoint :: Vector{Int} -> Coordinates of the mesh endpoint

    Output:

        FemMesh :: RegularMesh -> Struct that contains the mesh
"""
function RegularMesh(::Type{D}, :: Type{B}, :: Type{M},
                    DimNodes :: Array{Int,1},
                    EndPoint :: Array{Float64,1} ) where {D <:AbstractDim, B <: AbstractBasis, M <: AbstractMaterial}

    # coef de nodos por cada base, para saber cuantos se cogen
    multi = getMulti(B)
    # total number of nodes
    NofNodes = prod(multi.*DimNodes .+ 1)
    DOF = NofNodes * getMatDof[M]
    # number of spatial dimensions
    Ndim = getSymDim(D)
    # increment in distance for each element
    dx = EndPoint ./ DimNodes

    #array containing coordinates
    Coordinates = zeros(Float64, Ndim, DimNodes.*multi.+1...)
    #array containing nodes
    Nodes = collect(reshape(1:NofNodes, DimNodes.*multi.+1...))

    Threads.@threads for I in CartesianIndices(Nodes)
        Coordinates[:,I] = (I.I .-1).*dx
    end

    Dim = D()

    Els = CreateRegularElements(D, B, M, DimNodes, Nodes, Coordinates)
    #Els::Array{AbstractElement,1}= []
    return RegularMesh(DimNodes, EndPoint, Els, NofNodes, DOF, Nodes, Coordinates, Dim)
end

# TODO: Add the correct input types

"""
    function CreateRegularElements(::Type{D},::Type{B},::Type{M}, DimNodes, Nodes, Coordinates)

    Function used to create the elements of a regular mesh

    Input:

        D :: Type                     -> Spatial dimensions of the mesh
        B :: Type                     -> Basis of the elements shape functions
        M :: Type                     -> Material of the elements
        DimNodes :: Vector{Int}       -> Number of elements in each dimension
        EndPoint :: Vector{Float64}   -> Coordinates of the mesh endpoint
        Coordinates :: Array{Float64} -> Coordinates of all the nodes

    Output:

        Els :: Array{AbstractElement,1} -> Array containing all the elements
"""
function CreateRegularElements(::Type{D},::Type{B},::Type{M}, DimNodes, Nodes,
    Coordinates) where {D <: AbstractDim, B <: AbstractBasis, M<:AbstractMaterial}

    Els::Array{AbstractElement,1} = []


    NofElements = prod(DimNodes)
    ElementTable = collect(reshape(1:NofElements, DimNodes...))

    for I in CartesianIndices(ElementTable)
        elNodes = getElementNodes(I.I,DimNodes)
        element_i = Element(D,B,M,ElementTable[I],elNodes)
        for k in 1:length(elNodes)
            index = CartesianIndices(Nodes)[elNodes[k]]
            element_i.Coords[k,:] = Coordinates[:,CartesianIndices(Nodes)[elNodes[k]]]
        end
        push!(Els, element_i)
    end
    return Els
end

"""
    function getElement(coords, Mesh)

    Get the element number inside which the point with the given coordinates is located
"""
function getElement(coords::Vector{Float64}, Mesh::RegularMesh)
    ElementDims = Mesh.EndPoint ./ Mesh.DimNodes

    #Code fragment for finding NaN bug
    if coords .÷ ElementDims == NaN
        prinln("NaN bug")
        println("Coords: $coords")
        println("ElementDims: $ElementDims, coords .÷ ElementDims = $(coords .÷ ElementDims)")
    end


    Lowerindices = floor.(Int64, coords .÷ ElementDims)
    fact = ones(Int64, size(Mesh.DimNodes))
    for i in 1:size(Lowerindices)[1]
        for j in 1:i-1
            fact[i] *= Mesh.DimNodes[j]
        end
    end
    index = Lowerindices'*fact+1
    return index
end

"""
    getElementNodes(I, DimNodes)

    Get the nodes of each element at a given index I

    Input:

        I :: CartesianIndex      -> Index of the element
        DimNodes :: Vector{Int}  -> Number of elements in each dimension
"""
function getElementNodes(I,DimNodes)
    Ndim = length(DimNodes)
    dx = zeros(Int,Ndim)
    for i in 1:Ndim
        dx[i] = prod(DimNodes[1:i-1].+1)
    end
    Corner = getCorner(Ndim)
    elNodes = zeros(Int,2^Ndim)

    for i in 1:2^Ndim
        elNodes[i] = (Corner[i] + collect(I).-1 )'*dx + 1
    end

    return elNodes
end


# TODO: Explicar esto mejor, buscar la hoja

"""
    getCorner(N::Int)

    Get all the corners of the element of a given dimension

    Input:

        NDims :: Int -> Number of spatial dimensions
"""
function getCorner(NDims::Int)
    if NDims > 2
        PrevCorner = getCorner(NDims-1)
        return vcat(push!.(deepcopy(PrevCorner),0),push!.(deepcopy(PrevCorner),1))
    elseif NDims == 2
        return [[0, 0], [1, 0], [1, 1], [0, 1]]
    elseif NDims == 1
        return [[0], [1]]
    end
end

# Used to create Quad9 meshes, TODO: OPTIMIZE this
function getMulti(El::Element{D,B,M}) where{D <: AbstractDim, B <: AbstractBasis,M <: AbstractMaterial}
    return getMulti(B)
end

function getMulti(::Type{B}) where {B <: AbstractBasis}
    return 1
end

function getMulti(::Type{Quad4})
    return 1
end

function getMulti(::Type{Brick8})
    return 1
end

function getMulti(::Type{Quad9})
    return 2
end

function getMulti(::Type{Lin3})
    return 2
end


"""
get_node_coordinates(node::Int64, Mesh::RegularMesh)

Function that returns an array with the corresponding of a specific node in a
given mesh using the "Coordinates" matrix of that mesh.

Examples
≡≡≡≡≡≡≡≡≡≡≡

julia> get_node_coordinates(1, FemMesh) #Where FemMesh is a 1D RegularMesh
1-element Array{Float64,1}:
 0.0

 julia> get_node_coordinates(1, FemMesh) #Where FemMesh is a 2D RegularMesh
 2-element Array{Float64,1}:
  0.0
  0.0

"""
function get_node_coordinates(node::Int64, Mesh::RegularMesh)
    indices = CartesianIndices(Mesh.Nodes)[node]
    coordinates = Mesh.Coordinates[:,indices]

    return coordinates
end

"""
get_local_coordinates(globalcoords, el)

Returns a vector with the local coordinates in the given element

Only usable for Quad4, Brick8 elements without angular deformation
"""
function get_local_coordinates(globalcoords::Vector{Float64}, el::AbstractElement)

    localcoords = (globalcoords-el.Coords[1,:]) ./ (el.Coords[end-1,:]-el.Coords[1,:])*2 .- 1

     #Comprobation local coords are OK:
     #=Nmat = getNmat(el, localcoords)
     if globalcoords - (Nmat*el.Coords)' != zeros(size(localcoords))
         println("Error in the local coordinates of the thenter of TipCell #$(TC.ID)")
     end
     =#

    return localcoords
end

"""
    get_value

    Returns the value of the Scalar field in the given point
"""
function get_value(localcoords::Vector{Float64}, el::AbstractElement, d_u::Vector{Float64})
    # Calculation VEGF gradient
    u_nodes = d_u[el.Nodes]
    Nmat = getNmat(el, localcoords)
    val = Nmat*u_nodes
    return val
end

"""
    get_gradient

    Returns the gradient of the Scalar field in the given point
"""
function get_gradient(localcoords::Vector{Float64}, el::AbstractElement, d_u::Vector{Float64})
    # Calculation VEGF gradient
    u_nodes = d_u[el.Nodes]

    Jac = ElementJacobian(el, localcoords)
    J_1 = inv(Jac)
    Bmat = getBmat(el, J_1, localcoords)
    grad = Bmat*u_nodes
    return grad
end

"""
    checkifinbounds(globalcoords, Mesh)

    Checks if the given point is inside the rectangular Mesh
    Returns a Bool type:
        · true: if the point is strictly inside the Mesh
        · false: if the point is outside or in the boundaries
"""
function checkifinbounds(globalcoords::Vector{Float64}, Mesh::RegularMesh)
    inbounds = true
    for i in 1:length(Mesh.EndPoint)
        if globalcoords[i]>=Mesh.EndPoint[i]
            inbounds = false
        elseif globalcoords[i] <= 0.
            inbounds = false
        end
    end
    return inbounds
end

"""
function distances_Matrix(Mesh::RegularMesh)

    Returns a NofNodes*NofNodes Symmetric matrix containing the euclidean distances
    between the nodes in the given Mesh

"""
function distances_Matrix(Mesh::RegularMesh)
    Nodes = Mesh.Nodes
    n = Mesh.NofNodes

    # Distances matrix stored in a 2 dymensional Array
    Distances = zeros(n,n)
    for i in Nodes
        coords_i = get_node_coordinates(i, Mesh)
        for j in 1:i
            coords_j = get_node_coordinates(j, Mesh)
            distance_ij = norm(coords_i-coords_j)
            #Distances[i,j] = distance_ij
            Distances[j,i] = distance_ij
        end
     end

    #=
    #Distances stored in a 2*N dymensional Array
    dim= size(Mesh.Nodes) #dimensions of nodes matrix
    Dim = vcat(collect(dim),collect(dim))   #dimensions of Distances array
    Distances = zeros(Tuple(Dim))
    for i in Nodes
        coords_i = get_node_coordinates(i, Mesh)
        indices_i = CartesianIndices(Nodes)[i]
        for j in Nodes
            coords_j = get_node_coordinates(j, Mesh)
            indices_j = CartesianIndices(Nodes)[j]

            distance_ij = norm(coords_i-coords_j)
            Distances[indices_i, indices_j] = distance_ij
        end
     end
    =#

    return Symmetric(Distances)
end

"""
    adjacent4(Mesh, distances_Matrix)

    Returns a dictionary containing the 4-adjacent nodes to each one of the Mesh.
    That is to say, the ones adjacent to each node in the 4 (6 in 3D) main directions.

    Example:

    M = distances_Matrix(Mesh)
    adj4 = adjacent4(Mesh,M)
    typeof(adj4) = Dict{Int,Array{Int64, 1}}

    For a 2D mesh with 3 nodes per row/column:

        adj4[1] => [2, 4]
        adj4[2] => [1, 3, 5]
        adj4[5] => [3, 4, 6, 8]

    If the mesh were 3D, with 3 nodes in each direction:

        adj4[1] => [2, 4, 10]
        adj4[14] => [5, 11, 13, 15, 17, 23]

"""
function adjacent4(Mesh::RegularMesh, distances_Matrix::Symmetric{Float64,Array{Float64,2}})
    Ndim = length(Mesh.DimNodes)
    d_4 = Mesh.EndPoint ./ Mesh.DimNodes
    Adjacent_nodes = Dict{Int,Vector{Int64}}()
    for i in 1:length(Mesh.Nodes)
        vector=empty([1])
        for j in 1:length(Mesh.Nodes)
            for k in 1:length(d_4)
                if (distances_Matrix[i,j] <= d_4[k]+200*eps())&(distances_Matrix[i,j] >= d_4[k]-200*eps())
                    if !(j in vector) push!(vector, j) end
                end
            end
        end
        Adjacent_nodes[i] = vector
    end
    #Comprobation of the number of boundary nodes
    numberofcornernodes = 0
    for i in 1:Mesh.NofNodes
        if length(Adjacent_nodes[i])==Ndim
            numberofcornernodes += 1
        end
    end
    real_NofCorners = 2^Ndim
    if real_NofCorners != numberofcornernodes
        println("Error in corner detection. $numberofcornernodes corners found when $real_NofCorners were expected")
    end

    numberofaristnodes = 0
    for i in 1:Mesh.NofNodes
        if length(Adjacent_nodes[i])<1+Ndim
            numberofaristnodes += 1
        end
    end
    real_NofAristNodes = 2^(Ndim-1)*dot(ones(Int64,length(Mesh.DimNodes)), Mesh.DimNodes .- 1) + real_NofCorners
    if real_NofAristNodes != numberofaristnodes
        println("Error in arist detection. $numberofaristnodes arist nodes found when $real_NofAristNodes were expected")
    end

    return Adjacent_nodes
end
