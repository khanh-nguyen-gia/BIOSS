
import BIOSSFemMaterial: TipCell, HypoxCell

"""
    function TipCell(ID, birthdate, Genealogy, coords, Mesh, d_EC)
        Constructor for TipCell material given
            · Birthdate and genealogical history of the TipCell
            · The initial position of the TipCell
            · The Endothelial cell distribution: the new cell will be orientated against
              the EC gradient in its center point
"""
function TipCell(ID::Int64, date::Float64, TC_Genealogy::Vector{Int64},
    TC_coords::Vector{Float64}, Mesh::RegularMesh, d_EC::Vector{Float64})

    TC_Element_ID = getElement(TC_coords, Mesh)
    TC_Element = Mesh.Elements[TC_Element_ID]
    TC_Nodes = TC_Element.Nodes
    TC_distances = [norm(Mesh.Elements[TC_Element_ID].Coords[i, :]-TC_coords) for i in 1:size(TC_Nodes)[1]]

    localcoords = get_local_coordinates(TC_coords, TC_Element)
    grad_EC = get_gradient(localcoords, TC_Element, d_EC)
    TC_orientation = -grad_EC/norm(grad_EC)

    return TipCell(ID, date, TC_Genealogy, TC_Element_ID, TC_Nodes, TC_coords, TC_orientation, TC_distances)
end

function TipCell(ID::Int64, date::Float64, Genealogy::Vector{Int64}, node:: Int64, Mesh::RegularMesh, d_EC::Vector{Float64})
    coords = get_node_coordinates(node, Mesh)
    return TipCell(ID,date, Genealogy, coords, Mesh, d_EC)
end

"""
    TipCell(ID, date, node, Mesh, d_EC, d_vessels, List)
        Constructor for TipCell material given:
            · ID
            · BirthDate
            · Initial node of the TipCell
            · Mesh: RegularMesh, to locate the TipCell
            · d_EC: to orientate the new TipCell against its gradient
            · d_vessels: Distribution of vessels in the given mesh, in order to know the mother vessel of the newborn
            · List: TipCellList containing the Genealogical history of the vessels in d_vessels, in order to obtain the new TC's
"""
function TipCell(ID::Int64, date::Float64, node:: Int64, Mesh::RegularMesh, d_EC::Vector{Float64}, d_vessels::Vector{Int64}, List::TipCellList)
    coords = get_node_coordinates(node, Mesh)
    mother_vessel = d_vessels[node]
        if mother_vessel == 0
            Genealogy = [0]
        else
            Genealogy = copy(List.Genealogy[mother_vessel])
            push!(Genealogy, mother_vessel)
        end


    return TipCell(ID,date, Genealogy, coords, Mesh, d_EC)
end

"""
    HypoxCell(ID, coords, Mesh)
    Constructor for HypoxCell material, when only its ID and position are known
"""
function HypoxCell(ID::Int64, coords::Vector{Float64}, Mesh::RegularMesh)
    nodes = empty([0])
    near_nodes = empty([0])
    HyC = HypoxCell(ID, coords, nodes, near_nodes)
    for node in Mesh.Nodes
        node_coords = get_node_coordinates(node, Mesh)
        distance = norm(node_coords-coords)
        if distance <= 6*HyC.radius #3 in case deactivation happens when any EC is nearby
            push!(HyC.Near_Nodes, node)
            if distance <= HyC.radius
                push!(HyC.Nodes, node)
            end
        end
    end
    return HyC
end

"""
    VesselInfo

    Data structure containing all the information about the vessels in angiogenesis:
    · distribution: In which nodes there is a vessel and which one it is
    · active_nodes: Wether blood circulation is active in a vessel or not
    · nodes_BDs: the time(date) in which the Tip cell that reated the vessel was in that specific location
    · notch_viability :: the information on which nodes cannot be branching points because of the notch pathway (a tip cell has been activated nearby)
"""
mutable struct VesselInfo <: AngiogenesisMat
    distribution::Vector{Int64}
    active_nodes::Vector{Int64}
    nodes_BDs::Vector{Float64}
    notch_viability :: Vector{Bool}
end



include("ContinuousBehaviour.jl")

include("Miscellaneous.jl")

include("TipCellBehaviour.jl")
