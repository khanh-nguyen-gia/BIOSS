"""
    TipCell

    Agent material for Endothelial tip cells in angiogenesis simulation
        Contains:
    · ID: ID number for the Tip Cell
    · BirthDate: time of creation
    · Genealogy: vector containing IDs of mother, grandmother... all the way to initial
    · Element: Id of the node where it is located at that moment
    · Nodes: containing the IDs of the nodes surrounding the tip cell center (nodes in Element)
    · coords: Given as a vectorof FLoat64 elements
    · orientation: Orientation of the tip Cell
    · distances: containing the distances from the tip cell center to the nodes of the element
    · velocity: Given as a vector of FLoat64 elements


        Also the following parameters
    · radius: Size of the tip cell
    · rho: Endothelial cell density in the area occupied by a tip cell
    · theta: Tip cell "width of sight" (max angle between velocity and orientation)
    · x_Vegf: Tip cells' chemotactic response coefficient (when VEGF concentration ~ 0)
    · G_max: Norm of VEGF gradient for maximum tip cell speed
    · G_min: Norm of minimum VEGF gradient for tip cell activation
    · r_Notch: Notch Signaling Pathway travel distance for Tip Cell fate supression close to active TCs

"""
mutable struct TipCell <: AngiogenesisMat
    ID :: Int64
    BirthDate :: Float64
    Genealogy :: Vector{Int64}

    Element :: Int64
    Nodes :: Vector{Int64}
    coords :: Vector{Float64}
    orientation :: Vector{Float64}
    distances :: Vector{Float64}
    velocity :: Vector{Float64}

    #Parameters
    radius :: Float64
    rho :: Float64
    theta :: Float64
    x_Vegf :: Float64
    G_max :: Float64
    G_min :: Float64
    r_Notch :: Float64
end


"""
TipCell(ID, BirthDate, Element, Nodes, coords, orientation, distances)

Constructor for TipCell agent Material when birth, node and position are known

"""


function TipCell(ID::Int64, BirthDate::Float64, Genealogy::Vector{Int64}, Element::Int64, Nodes::Vector{Int64},
     coords::Vector{Float64}, orientation::Vector{Float64}, distances::Vector{Float64})

    Ndim = size(coords)[1] # get number of spatial dimensions

    radius = 5. # in μm
    radius /= L0
    r_Notch = 8*radius

    rho = 1.
    theta = deg2rad(60.)

    #Coefficients provided by Anderson & Chaplain, 1998
    max_speed = 21 #tip cell speed in μm/h
    G_max = 0.000000000001 #gradient for max speed in μm⁻¹

    max_speed = max_speed/L0*(T0/3600) #dimensionless parameter
    G_max *= L0

    G_min = G_max/2

    x_Vegf = max_speed/G_max



    velocity = x_Vegf*orientation/norm(orientation)


    return TipCell(ID, BirthDate, Genealogy, Element, Nodes, coords, orientation, distances, velocity, radius, rho, theta, x_Vegf,G_max, G_min, r_Notch)
end

function TipCell(ID::Int64, Element::Int64, Nodes::Vector{Int64},
     coords::Vector{Float64}, orientation::Vector{Float64}, distances::Vector{Float64})

     BirthDate = 0.
     Genealogy = [0]

     return TipCell(ID, BirthDate, Genealogy, Element, Nodes, coords, orientation, distances)
 end

"""
Constructor for a generic TipCell containing no information on its position, velocity
"""

function TipCell(ID::Int64, BirthDate::Float64, Genealogy::Vector{Int64}, D::Type{Dim}) where {Dim<:AbstractDim}

    Ndim = getSymDim(D) # get number of spatial dimensions

    Element = 1
    coords = zeros(Ndim)
    orientation = zeros(Ndim)
    velocity = zeros(Ndim)
    Nodes = zeros(Int64, 2^Ndim)
    distances = zeros(2^Ndim)

    return TipCell(ID, BirthDate, Genealogy, Element, Nodes, coords, orientation, distances)
end

function TipCell(ID::Int64, D::Type{Dim}) where {Dim<:AbstractDim}
    BirthDate = 0.
    Genealogy = [0]
    return TipCell(ID, BirthDate, Genealogy, D)
end

"""
getChemotacticCoefficient(TC::TipCell, c::Float64)

Calculates the chemotactic response of a Tip cell
    • TC: TipCell material containing the parameters involved
    • c: VEGF concentration value in the point of study
"""
function getChemotacticCoefficient(TC::TipCell, c::Float64)
    x = TC.x_Vegf
    α = TC.G_max

    return x/(1+α*c)
end


"""
    TipCellList

    Agent material for the List of Endothelial tip cells in angiogenesis simulation
        Contains:
    · IDs: vector containing ID number of each Tip Cell
    · BirthDates: Dictionary contining BirthDates of its TipCells tied to their ID
    · Genealogy: Dictionary containing vectors with the IDs of mother, grandmother... for each vessel

    · Nodes: Vector containind Id of the node where it is located at that moment
    · coords:  List of coordinates of the tip cells (Vector of FLoat64 vectors)
    · orientations: List of orientations of the tip Cells
    · velocities: List of velocities of the tip cells (Vector of FLoat64 vectors)

        Also the following parameters
    · radius: Size of the tip cells
    · rho: Endothelial cell density in the area occupied by a tip cell
    · theta: Tip cell "width of sight" (max angle between velocity and orientation)
    · x_Vegf: Tip cells' chemotactic response coefficient (when VEGF concentration ~ 0)
    · G_max: Norm of VEGF gradient for maximum tip cell speed
    · G_min: Norm of minimum VEGF gradient for tip cell activation
    · r_Notch: Notch Signaling Pathway travel distance for Tip Cell fate supression close to active TCs



"""
mutable struct TipCellList <: AngiogenesisMat
    IDs :: Vector{Int64}
    BirthDates :: Dict{Int, Float64}
    Genealogy :: Dict{Int, Vector{Int64}}
    Elements :: Vector{Int64}
    Nodes :: Array{Int64,2}
    coords :: Vector{Array{Float64,1}}
    orientations :: Vector{Array{Float64,1}}
    distances :: Array{Float64,2}
    velocities :: Vector{Array{Float64,1}}

    #Parameters
    radius :: Float64
    rho :: Float64
    theta :: Float64
    x_Vegf :: Float64
    G_max :: Float64
    G_min :: Float64
    r_Notch :: Float64
end


"""
function TipCellList(List::Vector{TipCell})

Constructor for a TipCellList material given a list of TipCells

"""
function TipCellList(List::Vector{TipCell})
    IDs = [i.ID for i in List]
    Elements = [i.Element for i in List]
    coords = [i. coords for i in List]
    orientations = [i.orientation for i in List]

    BirthDates = Dict{Int64,Float64}()
    Genealogy = Dict{Int, Vector{Int64}}()
    Nodes = Array{Int64}(undef, 0, size(List[1].Nodes)[1])
    distances = Array{Float64}(undef, 0, size(List[1].Nodes)[1])
    for i in 1:size(List)[1]
        BirthDates[IDs[i]] = List[i].BirthDate
        Genealogy[IDs[i]] = copy(List[i].Genealogy)
        Nodes = vcat(Nodes, List[i].Nodes')
        distances = vcat(distances, List[i].distances')
    end

    velocities = [i.velocity for i in List]
    radius = List[1].radius
    rho= List[1].rho
    theta = List[1].theta
    x_Vegf = List[1].x_Vegf
    G_max = List[1].G_max
    G_min = List[1].G_min
    r_Notch = List[1].r_Notch

    return TipCellList(IDs, BirthDates, Genealogy, Elements, Nodes, coords, orientations, distances, velocities, radius, rho, theta, x_Vegf, G_max, G_min, r_Notch)
end

"""
function TipCellList(D::Type{Dim}) where Dim<: AbstractDim

Generic TipCellList constructor of List with no elements
Contains only the TipCell parameters, obtained from a generic TipCell

"""
function TipCellList(D::Type{Dim}) where Dim<: AbstractDim
    Ndim = getSymDim(D)
    IDs = empty([1])
    BirthDates = Dict{Int64,Float64}()
    Genealogy = Dict{Int, Vector{Int64}}()
    Elements = empty([1])
    Nodes = Array{Int64}(undef, 0, 2^Ndim)
    coords = empty([[0.]])
    orientations = empty([[0.]])
    distances = Array{Float64}(undef, 0, 2^Ndim)
    velocity = empty([[0.]])

    #Tip-cell parameters are obtained from a generic TipCell
    generic_TC = TipCell(0,D)

    radius = generic_TC.radius
    rho = generic_TC.rho
    theta = generic_TC.theta
    x_Vegf = generic_TC.x_Vegf
    G_max = generic_TC.G_max
    G_min = generic_TC.G_min
    r_Notch = generic_TC.r_Notch

    return TipCellList(IDs, BirthDates, Genealogy, Elements, Nodes, coords, orientations, distances, velocity, radius, rho, theta, x_Vegf,G_max, G_min, r_Notch)
end

function TipCellList(D::Type{Dim}, TC_vector::Vector{TipCell}) where Dim<: AbstractDim
    List = TipCellList(D)
    expandTipCellList!(List, TC_vector)
    return List
end



"""
expandTipCellList!(List, List2) -> List

Adds the tip cells in List2 to List

Has 2 methods:

expandTipCellList!(List:: TipCellList, TC_vector::Vector{TipCell})

Argument TC_vector is a vector containing the endothelial TipCells that are to be added to List

"""
function expandTipCellList!(List:: TipCellList, TC_vector::Vector{TipCell})

    List2= TipCellList(TC_vector)
    return expandTipCellList!(List, List2)
end

"""
expandTipCellList!(List:: TipCellList, List2:: TipCellList)

Both arguments are TipCellList, storages the data from the second one into the first one
"""
function expandTipCellList!(List:: TipCellList, List2:: TipCellList)

    if List.IDs != setdiff(List.IDs, List2.IDs)
        println("At least one of the tip cells that are to be added has the same ID as some of the main list.")
        println("This will produce conflict and errors going forward.")
    end
    List.IDs = vcat(List.IDs, List2.IDs)
    for ID in List2.IDs
        List.BirthDates[ID] = List2.BirthDates[ID]
        List.Genealogy[ID] = copy(List2.Genealogy[ID])
    end
    List.Elements = vcat(List.Elements, List2.Elements)
    List.Nodes = vcat(List.Nodes, List2.Nodes)
    List.distances = vcat(List.distances, List2.distances)
    List.coords = vcat(List.coords, List2.coords)
    List.orientations = vcat(List.orientations, List2.orientations)
    List.velocities = vcat(List.velocities, List2.velocities)

    return List
end

"""
    TipCellDeactivation!

    Returns the TipCell List without the TipCell(s) whose ID(s) are introduced
"""
function TipCellDeactivation!(ID::Int64, List::TipCellList)
    println("")
    println("Deacivating TipCell #$ID.")

    indices = findall(x -> x!=ID, List.IDs) #indices of the REMAINING TCs
    if length(indices)==length(List.IDs)

        println("TipCell #$ID has not been found in the List")
        if Base.isempty(List.IDs)
            println("...Because the List is empty")
        else
            println("TipCells in the List: $(List.IDs)")
        end

    else

        List.IDs = List.IDs[indices]
        List.Elements = List.Elements[indices]
        List.Nodes = List.Nodes[indices, :]
        List.coords = List.coords[indices]
        List.orientations = List.orientations[indices]
        List.distances = List.distances[indices, :]
        List.velocities = List.velocities[indices]

        if Base.isempty(List.IDs)
            println("All the tip cells in the List have been deactivated")
        else
            println("Remaining TipCells: $(List.IDs)")
        end

    end
    return List
end

function TipCellDeactivation!(IDs::Vector{Int64}, List::TipCellList)

    println("")
    println("Deactivating TipCells: #$(IDs[:]).")
    Remaining_IDs = setdiff(List.IDs, IDs)

    if Remaining_IDs == List.IDs
        println("None of the TipCells have been found in the List")
        if Base.isempty(List.IDs)
            println("...Because the List is empty")
        else
            println("TipCells in the List: $(List.IDs)")
        end
    else
        if length(Remaining_IDs)!=(length(List.IDs)-length(IDs)) #Not all TCs have been found
            found_IDs = setdiff(List.IDs, Remaining_IDs)
            nonfound_IDs = setdiff(IDs, found_IDs)
            println("The following TipCells have not been found in the List:
            #$nonfound_IDs")
        end

        indices = [findfirst(x -> y==x, List.IDs) for y in Remaining_IDs] #indices of Remaining IDs
        List.IDs = List.IDs[indices]
        List.Elements = List.Elements[indices]
        List.Nodes = List.Nodes[indices, :]
        List.coords = List.coords[indices]
        List.orientations = List.orientations[indices]
        List.distances = List.distances[indices, :]
        List.velocities = List.velocities[indices]

        if Base.isempty(List.IDs)
            println("All the tip cells in the List have been deactivated")
        else
            println("Remaining TipCells: #$(List.IDs)")
        end
    end

    return List
end

"""
getChemotacticCoefficient(TC_List::TipCellList, c_vector::Vector{Float64})

Calculates the chemotactic response of all the Tip cells of the list
    • TC_List: TipCell material containing the parameters involved
    • c_vector: vector containing VEGF concentration value in each TipCell
"""
function getChemotacticCoefficient(TC_List::TipCellList, c_vector::Vector{Float64})
    if size(TC_List.IDs)!=size(c_vector)
        println("Error associating VEGF concentration to TipCell list")
    end
    x = TC_List.x_Vegf
    α = TC_List.G_max

    return [x/(1+α*c) for c in c_vector]
end


"""
    mothervessel
    Returns the "mother" vessel of the indicated TipCell, i.e., the vessel from which it branched
"""
function mothervessel(TC::TipCell)
    return TC.Genealogy[end]
end

function mothervessel(TC_ID::Int64, List::TipCellList)
    return List.Genealogy[TC_ID][end]
end
