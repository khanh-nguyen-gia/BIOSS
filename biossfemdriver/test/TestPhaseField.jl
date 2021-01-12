@time begin
    println("Compilando paquetes")

    using Revise

    using BIOSSFemCore
    using BIOSSFemMaterial
    using BIOSSFemElements
    using BIOSSFemDriver

    #using Plots
    using LinearAlgebra
    using SparseArrays

    include("developing_functions.jl")
    dir= "Desktop/Phase-Field_nuevo"
end

@time begin
    println("Creacion de mallas y condiciones de contorno")

    init_getNmat(_2D,Quad4, Endothelium) # inicio de la funcion de forma
    init_getNmat(_2D,Quad4, TAF)
    init_getNmat(_2D, Quad4, Fibronectin)

    init_getBmat(_2D,Quad4, Endothelium) # inicio de la funcion cinematica
    init_getBmat(_2D,Quad4, TAF)
    init_getBmat(_2D,Quad4, Fibronectin)

    initJacobian(_2D,Quad4)

    N = 180 # elementos en x
    M = 180 # elementos en y

    x̂ = 300. # x_max en μms adimensionalizada con L0 = 1.25μm
    ŷ = 300. # y_max en μms adimensionalizada con L0 = 1.25μm

    Mesh = RegularMesh(_2D, Quad4, Endothelium, [N, M], [x̂, ŷ]) # Malla
    MeshTAF = RegularMesh(_2D, Quad4, TAF, [N, M], [x̂, ŷ]) # Malla

    distances = distances_Matrix(Mesh) #Matrix containing distances between nodes of the mesh
    Adjacent_nodes = adjacent4(Mesh, distances)


    BcEC1 = CreateBoundary(Mesh, Endothelium, [1.], BCLeft, BCDirichlet)
    BcEC2 = CreateBoundary(Mesh, Endothelium, [1.], BCRight, BCDirichlet)
    BcEC = [BcEC1, BcEC2]
    modelEC = FemModel("Endothelial cells", "2D", Mesh, BcEC)

    BcTAF1 = CreateBoundary(MeshTAF, TAF, [0.], BCLeft, BCForces)
    BcTAF = [BcTAF1]
    modelTAF = FemModel("TAF", "2D", MeshTAF, BcTAF)
end
    ## #################################################3
@time begin
    println("Comienzo de la simulacion")

    d0_EC, d0_TAF, HyC_vec, modelTAF = initial_cond_PhaseField(Mesh, modelTAF)

    vessels = init_VesselInfo(d0_EC)

    d_vessels = vessels.distribution
    active_nodes = vessels.active_nodes
    nodes_BDs = vessels.nodes_BDs


    #se guarda la condicion inicial en Paraview
    println("Guardando el instante t = 0")
    SaveToParaview(Mesh,[d0_EC], ["Cells"], "cell0", dir; replace=true,append=true)
    SaveToParaview(Mesh,[d0_TAF], ["TAF"], "TAF0",dir; replace=true,append=true)
    #SaveToParaview(Mesh,[d0_Fiber], ["Fibronectin"], "Fiber0",dir; replace=true,append=true)


    #=  Graficar valor de μ en funcion de ρ para un valor de concentracion de TAF dado
    Capillar_material= Endothelium(_2D)
    ρ = collect(-2.:0.01:2.)
    f = fill!(copy(ρ), 0.001)
    μ = BIOSSFemDriver.ChemicalContribution(Capillar_material, ρ, f)

    plot(ρ, μ)
    =#


    #Calculo de la derivada temporal de la densidad de celulas endoteliales
    Dismantle!(modelEC)
    Assembly!(modelEC)
    f_EC = Forces(modelEC)

    AssemblyMass!(modelEC)

    # Escala temporal adimensionalizada con t0 = 5460s = 91 min
    t= 0.
    tmax = 22.0
    dt=0.05
    dt1 = 0.001 #for TAF evolution in the initial time step
    dt2 = 0.01 #for TAF evolution


    List = TipCellList(_2D)
    orientation=true
    #modelTAF, d0_TAF = TAF_evolution(modelTAF, d0_TAF, d0_EC, HyC_vec, 2.35, dt2)
    for i in 1:Int64(tmax÷dt)
        global modelEC, modelTAF, List, distances
        global d0_EC, d0_TAF, HyC_vec, d_vessels, active_nodes, vessels, dt, dt1, dt2, t

        t += dt
        println("Operando el instante t = $t")
        Find_irrigated_HypoxCells!(HyC_vec, vessels, List)

        if i == 1
            modelTAF, d0_TAF = Update_TAF_evolution!(modelTAF, d0_TAF, d0_EC, HyC_vec, dt, dt1)
        else
            modelTAF, d0_TAF = Update_TAF_evolution!(modelTAF, d0_TAF, d0_EC, HyC_vec, dt, dt2)
        end

        Update_TipCell_chemotaxis!(List, Mesh, d0_TAF, t, orientation)
        Update_TipCell_anastomosis!(List, Mesh, distances, t, d0_EC, d_vessels)
        TipCell_displacement!(List, Mesh, dt)
        update_rho_TipCells!(List, Mesh, distances, d0_EC, vessels, t)
        Update_VesselInfo!(d0_EC, vessels, Adjacent_nodes)
        BranchingNodes = branching_nodes_PF(List, Mesh, distances, Adjacent_nodes, vessels, d0_TAF, t)
        Branching!(List, t, BranchingNodes, Mesh, distances, d0_EC, d0_TAF, vessels)
        update_rho_TipCells!(List, Mesh, distances, d0_EC, vessels, t)
        Update_VesselInfo!(d0_EC, vessels, Adjacent_nodes)
        v0_EC = time_derivative_Endothelium(modelEC, f_EC, d0_EC, d0_TAF)
        d0_EC = d0_EC + v0_EC*dt

        SaveToParaview(Mesh,[d0_TAF], ["TAF"], "TAF$i", dir; replace=true,append=true)
        SaveToParaview(Mesh,[d0_EC], ["Cells"], "cell$i", dir; replace=true,append=true)
        SaveToParaview(Mesh,[d_vessels], ["Vessels"], "vessel$i", dir; replace=true,append=true)
        SaveToParaview(Mesh,[nodes_BDs], ["Birthdates"], "Birthdates$i", dir; replace=true,append=true)
        SaveToParaview(Mesh,[active_nodes], ["Active nodes"], "active$i", dir; replace=true,append=true)
    end
end






#=
#Comprobacion de qué nodos ha detectado como aristas
d_comprobacion = zeros(length(d0_EC))
nodos = collect(1:length(d0_EC))
Ndim = 2
filter!(j -> length(Adjacent_nodes[j])<=(2^(Ndim-1)+1) , nodos)
d_comprobacion[nodos] .= 1.
SaveToParaview(Mesh,[d_comprobacion], ["comprovacion"], "comprobacion", dir; replace=true,append=true)
=#

#Filtrado de los nodos que no cumplen condiciones para aparecer TipCells en ellos

#PossibleNodes = filter_nodes_PF(List, Mesh, distances, d0_EC, d0_TAF)
