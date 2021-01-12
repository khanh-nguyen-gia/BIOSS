mutable struct NastranContent
    ID :: Vector{Int}
    ElementNodes :: Vector{Vector{Int}}
    GridNodes :: Vector{Int}
    GridCoords :: Array{Float64,2}
end


"""

    function ImportNastranMesh(path)

    Import a .bdf mesh created in Nastran


    Input:

        inputfile :: String
"""
function ImportNastranMesh(dims, inputfile::String)


    file = open(inputfile,"r")

    text = readlines(file)

    escribe = false
    grid = true

    totalNodes = zeros(Int,0)
    Id = zeros(Int,0)
    Nodes = []

    if dims == 2
        totalCoords = zeros(Float64,2)
    end

    for line in text

        if occursin("BULK",line)
            escribe = true
            println("Comenzando a leer elementos")
            println("")
        elseif occursin("Material",line)
            escribe = false
        end

        if escribe

            if occursin("CQUAD4",line)

                el_id, el_nodes = parseElement(line)
                push!(Id, el_id)
                push!(Nodes, el_nodes)
            end

            if occursin("GRID",line)
                grid = true
                node, coords = parseGrid(dims,line)
                totalCoords = cat(dims = 2, totalCoords, coords)
                totalNodes = vcat(totalNodes,node)
            end

        end
    end

    close(file)

    NasContent = NastranContent(Id, Nodes, totalNodes, totalCoords[:,2:end])

    for i in 1:length(NasContent.ElementNodes)
        NasContent.ElementNodes[i] .-= minimum(NasContent.GridNodes) - 1
    end
    NasContent.GridNodes .-= minimum(NasContent.GridNodes) - 1

    NasContent.ID .-= minimum(NasContent.ID) - 1

    return NasContent

end



function parseElement(line::String)
    Inc = 8
    startNodes = false
    id = parse(Int,line[1*8:2*8])
    nodes = []

    NumOfDiv = Int(floor(length(line)/Inc))


    for i in 3:NumOfDiv-1
        fragment = line[i*Inc+1:i*Inc+7]
        node = parse(Int,fragment)
        push!(nodes, node)
    end
    fragment = line[NumOfDiv*Inc+1:end]
    node = parse(Int,fragment)
    push!(nodes,node)

    return id, nodes
end

function parseGrid(Dims::Int ,line::String)

    Inc = 8
    startCoords = false

    Coords = zeros(Dims)
    Coord_index = 1

    node = parse(Int,line[1*8:2*8])

    for i in 3:Int(floor(length(line)/Inc))-1
        fragment = line[i*Inc+1:i*Inc+7]

        if occursin(".", fragment)
            startCoords = true
        end
        if startCoords
            Coords[Coord_index] = parse(Float64, fragment)
            Coord_index += 1
        end
    end
    return node, Coords
end
