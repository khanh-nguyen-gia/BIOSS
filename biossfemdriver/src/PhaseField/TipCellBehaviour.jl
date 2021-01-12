
"""
    Update_TipCell_chemotaxis!

    Returns the given TipCell, Vector of TipCells or TipCellList
    with their velocities updated according to chemotactic response
"""
function Update_TipCell_chemotaxis!(TC::TipCell, Mesh::RegularMesh, d_VEGF::Vector{Float64})
    if (date-LTC.BirthDate)> 2*TC.radius/(TC.x_Vegf*TC.G_min)
        VEGF_grad = get_gradient(TC, Mesh, d_VEGF)
        G_max = TC.G_max
        #Update TipCell Velocity
        G = norm(VEGF_grad)

        # METODO UNA VELOCIDAD PARA INMADURA Y OTRA PARA MADURA
        TC.velocity = TC.x_Vegf*G_max*VEGF_grad/G

        #= METODO VELOCIDAD PROPORCIONAL AL GRADIENTE
        if G < G_max
            TC.velocity = TC.x_Vegf*VEGF_grad
        else
            TC.velocity = TC.x_Vegf*G_max*VEGF_grad/G
        end
        =#
    else #underage

        #METODO UNA VELOCIDAD PARA INMADURA Y OTRA PARA MADURA
        TC.velocity = TC.velocity/norm(TC.velocity)*G_min*x_Vegf

        # METODO VELOCIDAD PROPORCIONAL AL GRADIENTE
        #TC.velocity = TC.velocity/norm(TC.velocity)*G*x_Vegf
    end
    return TC
end

function Update_TipCell_chemotaxis!(TC_vector::Vector{TipCell}, Mesh::RegularMesh, d_VEGF::Vector{Float64}, date::Float64)
    for TC in TC_vector
        Update_TipCell_chemotaxis!(TC, Mesh, d_VEGF,date)
    end
    return TC_vector
end

function Update_TipCell_chemotaxis!(List_TCs::TipCellList, Mesh::RegularMesh, d_VEGF::Vector{Float64}, date::Float64)

    deactivation_IDs = empty([1])
    G_max = List_TCs.G_max
    G_min = List_TCs.G_min
    theta = List_TCs.theta
    for i in 1:length(List_TCs.IDs)
        TC_ID = List_TCs.IDs[i]

        el_ID = List_TCs.Elements[i]
        el = Mesh.Elements[el_ID]
        globalcoords = List_TCs.coords[i]
        localcoords = get_local_coordinates(globalcoords,el)

        VEGF_val = get_value(localcoords, el, d_VEGF)

        if VEGF_val[1] >= 0.001 #the tip cell is receiving VEGF stimmulus
            VEGF_grad = get_gradient(localcoords, el, d_VEGF)
            G = norm(VEGF_grad)
            TC_orientation =List_TCs.orientations[i]

            if (date-List_TCs.BirthDates[TC_ID])> 2*List_TCs.radius/(List_TCs.x_Vegf*G_min)

                #To avoid dividing by 0
                if G > 200*eps()
                    VEGF_grad_dir = VEGF_grad/G
                else
                    VEGF_grad_dir = TC_orientation
                end


                #=Pequeña correcion para evitar dividir por cero en el caso en que el gradiente sea opuesto a la orientación, pequeño giro horario =#
                if norm(VEGF_grad_dir+TC_orientation)<= 20*eps()
                    VEGF_grad_dir[2] += 0.05*VEGF_grad_dir[1]
                    VEGF_grad_dir[1] -= 0.05*VEGF_grad_dir[2]
                    if (length(VEGF_grad_dir)>2) #evitar caso 3D, vertical
                            VEGF_grad_dir[2] += 0.05*VEGF_grad_dir[3]
                    end
                    VEGF_grad_dir /= G
                end

                if checkifinsight(VEGF_grad_dir, TC_orientation, theta)
                    vel_dir = VEGF_grad_dir #vel_dir is the unitary vector in the direction of the new velocity vector
                else
                    #deviation dir is the unitary vector which is perpendicular to TC_dir and intersects with VEGF_grad_dir
                    deviation_dir = VEGF_grad_dir/dot(VEGF_grad_dir, TC_orientation)-TC_orientation
                    deviation_dir /= norm(deviation_dir)

                    vel_dir = TC_orientation + (1.0-10*eps())*deviation_dir*tan(theta)
                    vel_dir /= norm(vel_dir)
                end

                #### METODO UNA VELOCIDAD PARA INMADURA Y OTRA PARA MADURA #########
                List_TCs.velocities[i] = List_TCs.x_Vegf*G_max*vel_dir

            else #underage TipCell
                vel_dir = TC_orientation

                #### METODO UNA VELOCIDAD PARA INMADURA Y OTRA PARA MADURA #########
                List_TCs.velocities[i] = List_TCs.x_Vegf*G_max*vel_dir
            end

            #### METODO VELOCIDAD PROPORCIONAL AL GRADIENTE #########
            #=
            if G <= G_max
                List_TCs.velocities[i] = List_TCs.x_Vegf*G*vel_dir
            else
                List_TCs.velocities[i] = List_TCs.x_Vegf*G_max*vel_dir
            end
            =#


        else #the tip cell does not sense any VEGF
            ID = List_TCs.IDs[i]
            println("Tip Cell #$ID no longer senses VEGF presence")
            println("Procceeding its deactivation.")
            push!(deactivation_IDs, ID)
        end
    end

    if !Base.isempty(deactivation_IDs)
        TipCellDeactivation!(deactivation_IDs, List_TCs)
    end

    return List_TCs
end

function Update_TipCell_chemotaxis!(List_TCs::TipCellList, Mesh::RegularMesh, d_VEGF::Vector{Float64}, date::Float64, orientation::Bool)
    if orientation
        Update_TipCell_chemotaxis!(List_TCs, Mesh, d_VEGF, date)
    else
        deactivation_IDs = empty([1])
        G_max = List_TCs.G_max
        G_min = List_TCs.G_min
        x_Vegf = List_TCs.x_Vegf
        for i in 1:length(List_TCs.IDs)
            TC_ID = List_TCs.IDs[i]
            el_ID = List_TCs.Elements[i]
            el = Mesh.Elements[el_ID]
            globalcoords = List_TCs.coords[i]
            localcoords = get_local_coordinates(globalcoords,el)

            VEGF_val = get_value(localcoords, el, d_VEGF)
            if VEGF_val[1] >= 0.001 #the tip cell is receiving VEGF stimmulus
                VEGF_grad = get_gradient(localcoords, el, d_VEGF)
                G = norm(VEGF_grad)
                    if G  > 200*eps()

                        #### METODO UNA VELOCIDAD PARA INMADURA Y OTRA PARA MADURA #########
                        if (date-List_TCs.BirthDates[TC_ID])> 2*List_TCs.radius/(List_TCs.x_Vegf*G_min)

                            List_TCs.velocities[i] = List_TCs.x_Vegf*(G_max/G)*VEGF_grad

                        else #underage
                            dir = List_TCs.velocities[i]/norm(List_TCs.velocities[i])
                            List_TCs.velocities[i] *= G_min*x_Vegf
                        end

                        #### METODO VELOCIDAD PROPORCIONAL AL GRADIENTE #########
                        #=
                        if (date-List_TCs.BirthDates[TC_ID])> 2*List_TCs.radius/(List_TCs.x_Vegf*G_min)
                            if G <= G_max
                                List_TCs.velocities[i] = List_TCs.x_Vegf*VEGF_grad
                            else
                                List_TCs.velocities[i] = List_TCs.x_Vegf*(G_max/G)*VEGF_grad
                            end
                        else #underage
                            dir = List_TCs.velocities[i]/norm(List_TCs.velocities[i])
                            if G <= G_max
                                List_TCs.velocities[i] *= G*x_Vegf
                            else
                                List_TCs.velocities[i] *= G_max*x_Vegf
                            end
                        end
                        =#
                    end
            else #the tip cell does not sense any VEGF
                    ID = List.IDs[i]
                    println("Tip Cell #$ID no longer senses VEGF presence")
                    println("Procceeding its deactivation.")
                    push!(deactivation_IDs, ID)
            end
        end

        if !Base.isempty(deactivation_IDs)
            TipCellDeactivation!(deactivation_IDs, List_TCs)
        end
    end

    return List_TCs
end


"""
    Update_TipCell_anastomosis!(List, Mesh, Distances,d_EC, d_vessels)

    Detects which TipCells are close enough to other vessels to detect them and
    move towards them in order to initiate anastomosis. In this case, it overwrites
    the velocity and direction of the tip cell towards the other vessel

"""
function Update_TipCell_anastomosis!(List::TipCellList, Mesh::RegularMesh, Distances::Symmetric{Float64,Array{Float64,2}},
    date::Float64, d_EC::Vector{Float64}, d_vessels::Vector{Int64})

    #TODO: Implement alternative method without TC orientation
    G_max = List.G_max
    G_min = List.G_min
    theta = List.theta
    anast_dist = 4*List.radius
    anastomosis_fixed = empty([0]) #List of TipCells that have begun anastomosis
    anastomosis_directions = Dict{Int,Vector{Float64}}() #Dictionary containing those directions
    for TC in 1:length(List.IDs)
        TC_ID = List.IDs[TC]
        if (date-List.BirthDates[TC_ID])> 2*List.radius/(List.x_Vegf*G_min)
            globalcoords = List.coords[TC]
            TC_orientation = List.orientations[TC]
            detected = false #will become true if a vessel is detected

            detected_nodes = filter(j -> distance_nodetoTC(j, TC_ID, List, Distances)<anast_dist, Vector(1:Mesh.NofNodes))
            #Filter those nodes associated with own vessel and initial conditions
            #TODO: change initial conditions for mother vessel once TC activation is implemented
            filter!(j -> (d_vessels[j]>0)&(d_vessels[j]!=TC_ID), detected_nodes)
            #Last comprobation that EC density is close to 1
            filter!(j -> (d_EC[j]>0.7), detected_nodes) #TODO change this to -0.9 once in final model
            #Only keep those in sight
            filter!(j -> checkifinsight((get_node_coordinates(j,Mesh)-globalcoords),TC_orientation,theta), detected_nodes)

            detected_distances = [distance_nodetoTC(j, TC_ID, List, Distances) for j in detected_nodes]

            if !(Base.isempty(detected_nodes)) # That TipCell has found another vessel
                push!(anastomosis_fixed, TC_ID)

                objective_node = detected_nodes[argmin(detected_distances)] #The one closest to the TC
                println("TipCell #$(TC_ID) is going towards vessel #$(d_vessels[objective_node])")
                direction = (get_node_coordinates(objective_node, Mesh)-globalcoords) #non unitary
                anastomosis_directions[TC_ID] = direction/norm(direction)
                #List.orientations[TC] = anastomosis_directions[TC_ID]
                List.velocities[TC] = norm(List.velocities[TC])*anastomosis_directions[TC_ID]
            end
        end
    end
    return List
end

"""
    TipCell_displacement!

    Time integration for TipCell displacement according to velocities and update of the list
"""
function TipCell_displacement!(List::TipCellList, Mesh::RegularMesh, dt::Float64)

    deactivation_IDs = empty([1])
    turning_radius = 2*List.radius
    #TODO: make this radius a property of TipCells, meve this parameter to BIOSSFemMaterial


    for i in 1:length(List.IDs)
        coords1 = List.coords[i] + dt*List.velocities[i]
        distance = norm(coords1-List.coords[i])
        orientation1 = (coords1-List.coords[i])/distance
        if !checkifinsight(orientation1, List.orientations[i], List.theta)
            println("Warning: TipCell #$(List.IDs[i]) is trying to move outside its area of sight")
        end

        theta_turn = 2*atan(0.5*dt*norm(List.velocities[i])/turning_radius) # max angle of turn for a tip cell, trigonometry

        if !checkifinsight(List.velocities[i], List.orientations[i], theta_turn)
            #in order to get softer turns
            vel_dir = List.velocities[i]/norm(List.velocities[i])
            #deviation dir is the unitary vector which is perpendicular to TC_dir and intersects with vel_dir
            deviation_dir = vel_dir/dot(vel_dir, List.orientations[i])-List.orientations[i]
            deviation_dir /= norm(deviation_dir)

            vel_dir1 = List.orientations[i] + deviation_dir*tan(theta_turn)
            vel_dir1 /= norm(vel_dir1)

            coords1 = List.coords[i] + dt*vel_dir1*norm(List.velocities[i])
            orientation1 = (coords1-List.coords[i])/norm(coords1-List.coords[i])

        end

        if checkifinbounds(coords1, Mesh)
            List.coords[i] = coords1
            if distance >= 2000*eps()
                List.orientations[i] = orientation1
            else
                println("TipCell # $(List.IDs[i]) has stopped moving")
                println("Procceeding to its deactivation.")
                push!(deactivation_IDs, List.IDs[i])
            end
            List.Elements[i] = getElement(List.coords[i], Mesh)
            List.Nodes[i, :] = Mesh.Elements[List.Elements[i]].Nodes
            List.distances[i,:] = [norm(Mesh.Elements[List.Elements[i]].Coords[j, :]-List.coords[i]) for j in 1:size(List.Nodes)[2]]
        else

            println("Tip Cell #$(List.IDs[i]) has reached the boundaries of the Mesh.")
            println("Procceeding to its deactivation.")
            push!(deactivation_IDs, List.IDs[i])
        end
    end

    if !Base.isempty(deactivation_IDs)
        TipCellDeactivation!(deactivation_IDs, List)
    end

    return List
end

"""
update_rho_TipCells!(List_TCs::TipCellList, Mesh::RegularMesh, Distances::Array{Float64,2}, d0_EC::Array{Float64,1})

Updates endothelial cell density in the nodes inside the Tip Cells of the list

"""
function update_rho_TipCells!(List_TCs::TipCellList, Mesh::RegularMesh,
    Distances::Symmetric{Float64,Array{Float64,2}}, d0_EC::Array{Float64,1})

    radius = List_TCs.radius
    for ID in List_TCs.IDs
        for j in Mesh.Nodes
            r = distance_nodetoTC(j, ID, List_TCs, Distances)
            if r<=radius
                rho = rho_distrib_TC(r/radius)
                if d0_EC[j]<rho     d0_EC[j]=rho       end
            end
        end
    end
    #=
    for i in 1:length(List_TCs.Nodes)
        node_i = List_TCs.Nodes[i]
        if (List_TCs.distances[i]<=List_TCs.radius)&(d0_EC[node_i]<List_TCs.rho)
             d0_EC[node_i]=List_TCs.rho
        end
        for j in Mesh.Nodes
            if ( (Distances[node_i,j]+List_TCs.distances[i])<=List_TCs.radius )&(d0_EC[j]<List_TCs.rho)
                d0_EC[j]=List_TCs.rho
            end
        end
    end=#
    return d0_EC
end

#TODO: change deactivation to also occur when reaching a vessel 0, but be conditioned by cell age
function update_rho_TipCells!(List_TCs::TipCellList, Mesh::RegularMesh,
    Distances::Symmetric{Float64,Array{Float64,2}}, d0_EC::Array{Float64,1},
    vessels::VesselInfo, date::Float64)

    d_vessels = vessels.distribution
    nodes_BDs = vessels.nodes_BDs
    radius = List_TCs.radius
    merging_IDs = empty([1])
    contact_nodes = empty([1])
    for TC in 1:length(List_TCs.IDs)
        reached = false
        reachedvessel = -1
        contact_node = 0
        TC_ID = List_TCs.IDs[TC]
        mother_ID = mothervessel(TC_ID, List_TCs)

        for i in 1:length(List_TCs.Nodes[TC,:])
            node_i = List_TCs.Nodes[TC,i]
            r = List_TCs.distances[TC,i]
            if r <= radius
                rho =  rho_distrib_TC(r/radius)
                if d0_EC[node_i]<rho
                    d0_EC[node_i]=rho
                end
                if (d_vessels[node_i]>-1) & (d_vessels[node_i]!=mother_ID) & (d_vessels[node_i]!=TC_ID) & (nodes_BDs[node_i]!=Inf64)
                    reached = true
                    reachedvessel = d_vessels[node_i]
                    contact_node = node_i
                elseif d_vessels[node_i] == -1
                     d_vessels[node_i] = TC_ID
                     nodes_BDs[node_i] = date
                end
            end
        end

        for j in setdiff(Mesh.Nodes, List_TCs.Nodes[TC,:])  #the other nodes in the mesh
            r = distance_nodetoTC(j, TC_ID, List_TCs, Distances)
            if r <= radius
                rho = rho_distrib_TC(r/radius)
                if d0_EC[j]<rho
                    d0_EC[j]=rho
                end
                if (d_vessels[j]>-1)&(d_vessels[j]!=TC_ID)&(d_vessels[j]!=mother_ID)
                    #in contact but still not close enough

                    #reached = true
                    #reachedvessel = d_vessels[j]
                elseif d_vessels[j] == -1
                    d_vessels[j] = TC_ID
                    nodes_BDs[j] = date
                end
            end
        end

        if reached
            println("TipCell #$TC_ID reached vessel #$reachedvessel.")
            println("Procceeding to its deactivation")
            push!(merging_IDs, TC_ID)
            push!(contact_nodes, contact_node)
        end
    end
    if !Base.isempty(merging_IDs)
        Merge!(merging_IDs, contact_nodes, List_TCs, vessels)
    end

    return d0_EC, List_TCs, vessels
end


"""
function branching_nodes_PF(TCList, Mesh, Distances, d0_EC, d0_VEGF)

    Returns a vector containing the nodes of the given Mesh that meet the criteria to activate a TipCell

"""
function branching_nodes_PF(TCList::TipCellList, Mesh::RegularMesh, Distances::Symmetric{Float64,Array{Float64,2}},
    Adjacent_nodes::Dict{Int64, Vector{Int64}}, vessels::VesselInfo, d0_VEGF::Vector{Float64}, date::Float64)

    d_vessels = vessels.distribution
    notch_viability =vessels.notch_viability
    PossibleNodes = reshape(copy(Mesh.Nodes), length(Mesh.Nodes))
    Ndim = length(Mesh.DimNodes)
    minAge = 4*TCList.radius/(TCList.x_Vegf*TCList.G_min)
    minDistance = TCList.r_Notch

    #filter!(j -> length(Adjacent_nodes[j])>(2^(Ndim-1)+1) , PossibleNodes) #Those in the limits of the mesh give error
    filter!(j -> d_vessels[j] >= 0, PossibleNodes)
    filter!(j -> notch_viability[j], PossibleNodes)
    nodesinsprouts = filter(j -> (d_vessels[j]>0), PossibleNodes)
    underage = filter(j -> (date-TCList.BirthDates[d_vessels[j]])<minAge, nodesinsprouts)
    setdiff!(PossibleNodes, underage)
    filter!(j -> d0_VEGF[j] > 0.001, PossibleNodes)

    #Eliminate nodes within 8 radii from the closest active TC (Notch Pathway)
    for i in TCList.IDs
        filter!(j -> distance_nodetoTC(j,i, TCList, Distances) > minDistance, PossibleNodes)
    end

    #Finally, eliminate those nodes who are surrounded in a 65% by endothelium
    nodes_too_surrounded = empty([1])
    for i in PossibleNodes
        surrounding_nodes = filter(j -> Distances[i,j]<2.5*TCList.radius, Mesh.Nodes)
        surrounding_nodes_vessels = filter(j->d_vessels[j]>-1, surrounding_nodes)
        if length(surrounding_nodes_vessels)/length(surrounding_nodes)>= 0.55
            push!(nodes_too_surrounded, i)
        end
    end
    setdiff!(PossibleNodes, nodes_too_surrounded)

    branching_nodes = empty([1])
    while !Base.isempty(PossibleNodes)
        #Randomly select one of the list using a Uniform distribution
        node_i = PossibleNodes[ceil(Int64, rand(Uniform(0+eps(), length(PossibleNodes))))]
        push!(branching_nodes, node_i)
        filter!(j -> Distances[node_i, j] > minDistance, PossibleNodes)
    end
    return branching_nodes
end

"""
    Branching!

    Evaluates the nodes in which new TipCells activate and adds them to the List
"""
function Branching!(List::TipCellList, date::Float64, nodes:: Vector{Int64},
    Mesh::RegularMesh, Distances::Symmetric{Float64,Array{Float64,2}},
    d_EC::Vector{Float64}, d_TAF::Vector{Float64}, vessels::VesselInfo)

    G_max = List.G_max #Gradient norm for max speed
    G_min =List.G_min #Minimum gradient for Tip Cell activation
    d_vessels = vessels.distribution
    notch_viability =vessels.notch_viability
    minDistance = List.r_Notch
    AllNodes = reshape(Mesh.Nodes, length(Mesh.Nodes))

    if Base.isempty(List.BirthDates)
        next_ID = 1 # ID for the next tip cell to create
    else
        next_ID = maximum(List.BirthDates)[1] + 1
    end

    #Create an empty Vector{TipCell}
    TC_vector = [TipCell(0,_2D)]
    filter!(e -> e.ID != 0, TC_vector)
    final_nodes = empty([0])

    for node in nodes
        ID = next_ID
        TC_i = TipCell(ID, date, node, Mesh, d_EC, d_vessels, List)
        Grad_TAF = get_gradient(TC_i, Mesh, d_TAF)
        #if norm(Grad_TAF)>= G_min
            push!(TC_vector, TC_i)
            push!(final_nodes, node)
            next_ID += 1
            Notch_nodes = filter(j -> Distances[node, j] <= minDistance/2, AllNodes)
            notch_viability[Notch_nodes] .= false
        #end
    end

    if !Base.isempty(TC_vector)
        println("Branching in process. New TipCells: # $([i.ID for i in TC_vector]) in nodes #$final_nodes")
        expandTipCellList!(List, TC_vector)
    end


    return List, vessels
end

"""
    Merge!

    Function that simulates the biological process that occur when a TipCell
    merges with another vessel in the anastomosis process:
    · Activates blood circulation in the nodes where it is then viable
    · Deactivates the TipCell phenotype in the cells that are now connected to another vessel
"""
function Merge!(Merging_IDs::Vector{Int64}, contact_nodes::Vector{Int64},
    List::TipCellList, vessels::VesselInfo)

    d_vessels = vessels.distribution
    nodes_BDs = vessels.nodes_BDs
    active_nodes = vessels.active_nodes
    nodes = collect(1:length(d_vessels))


    for i in 1:length(Merging_IDs)
        TC_ID = Merging_IDs[i]
        Genealogy_TC = List.Genealogy[TC_ID]

        #activate nodes in TC's vessel
        nodes2activate = filter(j -> d_vessels[j]==TC_ID, nodes)

        #activate nodes in TipCell ancestors
        sprout_ID = TC_ID
        for j in length(Genealogy_TC):-1:2
            mother_ID = Genealogy_TC[j]
            sprout_BD = List.BirthDates[sprout_ID]
            activate_in_mother = filter(k -> (d_vessels[k]==mother_ID)&(nodes_BDs[k]<=sprout_BD), nodes)
            nodes2activate = vcat(nodes2activate, activate_in_mother)
            sprout_ID = mother_ID
        end

        contact_node = contact_nodes[i]
        reached_ID = d_vessels[contact_node]
        reached_BD = nodes_BDs[contact_node] #BD of the reached NODE
        if reached_ID == 0
            reached_Genealogy = [0]
        else
            reached_Genealogy = List.Genealogy[reached_ID]
        end

        #activate nodes in the reached vessel
        activate_in_reached = filter(k -> (d_vessels[k]==reached_ID)&(nodes_BDs[k]<=reached_BD), nodes)
        nodes2activate = vcat(nodes2activate, activate_in_reached)

        #activate nodes in the reached vessel's ancestors
        sprout_ID = reached_ID
        for j in length(reached_Genealogy):-1:2
            mother_ID = reached_Genealogy[j]
            sprout_BD = List.BirthDates[sprout_ID]
            activate_in_mother = filter(k -> (d_vessels[k]==mother_ID)&(nodes_BDs[k]<=sprout_BD), nodes)
            nodes2activate = vcat(nodes2activate, activate_in_mother)
            sprout_ID = mother_ID
        end

        active_nodes[nodes2activate] .= 1
    end

    TipCellDeactivation!(Merging_IDs, List)

    return List, vessels
end
