"""
    HypoxCell

    Agent material for Hypoxic cells in angiogenesis simulation
        Contains:
    · ID: ID number for the Hypoxic Cell
    · coords: Given as a vectorof FLoat64 elements
    · Nodes: containing the IDs of the nodes inside the Hypoxic Cell
    · active: True if active, false if inactive

        Also the following parameters
    · radius: Size of the hypoxic cell
    · TAF_hyc: TAF concentration inside the hypoxic cell
"""
mutable struct HypoxCell <: AngiogenesisMat
    ID :: Int64

    coords :: Vector{Float64}
    Nodes :: Vector{Int64}
    Near_Nodes :: Vector{Int64}
    active :: Bool

    #Parameters
    radius :: Float64
    TAF_hyc :: Float64
end


"""
HypoxCell(ID, coords, Nodes)

Constructor for HypoxCell agent Material when nodes and position are known

"""
function HypoxCell(ID::Int64, coords::Vector{Float64}, Nodes::Vector{Int64}, Near_Nodes::Vector{Int64})
    radius = 5.
    radius /= L0

    TAF_hyc = 1.
    active = true

    return HypoxCell(ID, coords, Nodes, Near_Nodes, active, radius, TAF_hyc)
end

"""
    HypoxCellDeactivation!

    Deactivates the indicated HypoxicCells
"""
function HypoxCellDeactivation!(ID::Int64, HyC_Vector::Vector{HypoxCell})
    println("Deactivating Hypoxic Cell #$ID")
    for i in HyC_Vector
        if i.ID == ID
            if i.active == true
                i.active = false
                println("Hypoxic Cell #$ID succesfully deactivated")
            else
                println("Hypoxic Cell #$ID had already been deactivated")
            end
        end
    end
    return HyC_Vector
end
function HypoxCellDeactivation!(IDs::Vector{Int64}, HyC_Vector::Vector{HypoxCell})
    for id in IDs
        HypoxCellDeactivation!(id, HyC_Vector)
    end
    return HyC_Vector
end
