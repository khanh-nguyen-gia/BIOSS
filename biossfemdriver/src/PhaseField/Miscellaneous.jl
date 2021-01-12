#Includes a varied set of functions used in Angiogenesis simulation

"""
    get_gradient(TC::TipCell, Mesh::RegularMesh, d_u::Vector{Float64})

    method for calculating the value of the gradient of d_u in the center of the TipCell
"""
function get_gradient(TC::TipCell, Mesh::RegularMesh, d_u::Vector{Float64})

    # Cell position in local coordinates
    el_ID = TC.Element
    el = Mesh.Elements[el_ID]
    localcoords = get_local_coordinates(TC.coords, el)
    return get_gradient(localcoords, el, d_u)
end

"""
    get_value(TC::TipCell, Mesh::RegularMesh, d_u::Vector{Float64})

    method for calculating the value of d_u in the center of the TipCell
"""
function get_value(TC::TipCell, Mesh::RegularMesh, d_u::Vector{Float64})

    # Cell position in local coordinates
    el_ID = TC.Element
    el = Mesh.Elements[el_ID]
    localcoords = get_local_coordinates(TC.coords, el)
    return get_value(localcoords, el, d_u)
end

"""
    checkifinsight(direction, TC_orientation, theta)

    Answers if direction is in the area of sight around TC_orientation delimited by theta
"""
function checkifinsight(direction::Vector{Float64}, TC_orientation::Vector{Float64}, theta::Float64)
    deviation_cos = dot(TC_orientation,direction)/(norm(TC_orientation)*norm(direction))
    if deviation_cos >= cos(theta) return true
    else return false
    end
end

"""
    distance_nodetoTC(Node, TC_ID, List, distances_Matrix)

    Returns the distance between a node and a Tip Cell of the given List
"""
function distance_nodetoTC(Node::Int64, TC_ID::Int64, List::TipCellList, distances_Matrix::Symmetric{Float64,Array{Float64,2}})
    index = findfirst(x -> x==TC_ID, List.IDs)
    alldistances = [List.distances[index,j]+distances_Matrix[Node, List.Nodes[index,j]] for j in 1:length(List.Nodes[index,:])]
    distance_min = minimum(alldistances)

    return distance_min
end

"""
    Find_irrigated_HypoxCells!(HyC_vec, VesselInfo)

    Deactivates those Hypoxic cells that are now irrigated from the Hypoxic cells vector
"""
function Find_irrigated_HypoxCells!(HyC_vec::Vector{HypoxCell}, vessels::VesselInfo, List::TipCellList)
    R = HyC_vec[1].radius
    active_nodes = vessels.active_nodes
    d_vessels = vessels.distribution
    deactivation_IDs = empty([0])
    for Cell in HyC_vec
        if Cell.active == true
            ID = Cell.ID
            deactivation = false
            for i in Cell.Near_Nodes
                if active_nodes[i] == 1 #when an active vessel is nearby
                #if d_vessels[i] > -1 #when a growing vessel is nearby
                    deactivation = true
                end
            end

            if deactivation == true
                println("Hypox Cell #$ID now has an active vessel in its proximity")
                println("Procceeding its deactivation.")
                push!(deactivation_IDs, ID)
            end
        end
    end

    HypoxCellDeactivation!(deactivation_IDs, HyC_vec)

    number_active_HyCs = 0
    for Cell in HyC_vec
        if Cell.active == true
            number_active_HyCs += 1
        end
    end

    if number_active_HyCs == 0 & !Base.isempty(List.IDs)
        println("All the initial hypoxic cells are now irrigated")
        println("procceeding to deactivate all tip cells")
        TipCellDeactivation!(List.IDs, List)
    end
    return HyC_vec, List
end
