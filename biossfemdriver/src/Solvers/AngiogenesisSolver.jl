#solver for angiogenesis problem, with VEGF is simulated by an steady state solution for a diffusion equation

function SolveProgramAngiogenesis(MeshVEGF, modelVEGF ::FemModel, MeshEC, modelEC::FemModel,
    dt:: Float64, TotalTime::Float64, dir::String; replace=true, append=true)

    if replace && ispath(homedir()*dirSep*dir)
        rm(homedir()*dirSep*dir, recursive=true)
    end

    d0_VEGF = SolveProgram(modelVEGF) #condicion inicial VEGF difundido por el tejido
    d0_EC = uval(modelEC)

    AssemblyMass!(modelVEGF)
    AssemblyMass!(modelEC)

    f_EC = Forces(modelEC)

    Nsteps = round(Int, TotalTime/dt)

    #d1_VEGF = d0_VEGF
    for i in 1:Nsteps

         for Element_i in modelEC.Elements
             Element_i.Mat.VEGF = d0_VEGF[Element_i.Nodes] # concentracion de VEGF en los nodos
         end

        # ActualizaVEGF(modeloVEGF,modeloCelulas)

        # for Element_i in model.Elements
        #     Element_i.Mat.Vel = grad(d0_VEGF) # gradiente de vegf
        # end
        Dismantle!(modelEC)
        Assembly!(modelEC)

        if i == 1
            global v0_EC
            v0_EC = sparse(modelEC.M) \ (f_EC - sparse(modelEC.K)*d0_EC)
        end

        d1_EC, v1_EC = dFormIntegration(modelEC, dt, f_EC, d0_EC, v0_EC, alpha = 1)
        v0_EC = v1_EC; d0_EC = d1_EC

        if mod(i,1) == 0
            println("Guardando el paso $i")
            SaveToParaview(MeshEC,[d1_EC], ["Cells"], "cell$i",
            dir; replace=replace,append=append)
        end
    end
end



function SolveProgramAngiogenesis2(MeshVEGF, modelVEGF ::FemModel, MeshEC, modelEC::FemModel,
    dt:: Float64, TotalTime::Float64, dir::String; replace=true, append=true)

    if replace && ispath(homedir()*dirSep*dir)
        rm(homedir()*dirSep*dir, recursive=true)
    end

    #To simulate initial diffusion of VEGF, EC concentration in every element must be 0
    d0_EC=zeros(MeshVEGF.DOF)
    for element_i in modelVEGF.Elements
        element_i.Mat.ECs=d0_EC[element_i.Nodes]
    end


    d0_VEGF = SolveProgram(modelVEGF) #condicion inicial VEGF difundido por el tejido
    d0_EC = uval(modelEC)


    #se guarda la condicion inicial en Paraview
    println("Guardando el paso 0")
    SaveToParaview(MeshEC,[d0_EC], ["Cells"], "cell0",
    dir; replace=replace,append=append)
    SaveToParaview(MeshVEGF,[d0_VEGF], ["VEGF"], "VEGF0",
    dir; replace=replace,append=append)


    DismantleMass!(modelVEGF)
    AssemblyMass!(modelVEGF)
    DismantleMass!(modelEC)
    AssemblyMass!(modelEC)

    f_EC = Forces(modelEC)
    f_VEGF = Forces(modelVEGF)

    Nsteps = round(Int, TotalTime/dt)

    #d1_VEGF = d0_VEGF
    for i in 1:Nsteps

        #Updating VEGF distribution in EC model
         for Element_i in modelEC.Elements
             Element_i.Mat.VEGF = d0_VEGF[Element_i.Nodes] # concentracion de VEGF en los nodos
         end

        Dismantle!(modelEC)
        Assembly!(modelEC)


        v0_EC = sparse(modelEC.M) \ (f_EC - sparse(modelEC.K)*d0_EC)

        d1_EC, v1_EC = dFormIntegration(modelEC, dt, f_EC, d0_EC, v0_EC, alpha = 1)
        d0_EC = d1_EC


        #Updating ECs distribution in VEGF model
        modelVEGF, d0_VEGF = updateVEGFmodel!(modelVEGF,d0_VEGF,d0_EC)

        #=
        for element_i in modelVEGF.Elements
            element_i.Mat.ECs=d0_EC[element_i.Nodes]

        end
        Dismantle!(modelVEGF)
        Assembly!(modelVEGF)
        =#

        v0_VEGF = sparse(modelVEGF.M) \ (f_VEGF - sparse(modelVEGF.K)*d0_VEGF)

        d1_VEGF, v1_VEGF = dFormIntegration(modelVEGF, dt, f_VEGF, d0_VEGF, v0_VEGF, alpha = 1)
        d0_VEGF = d1_VEGF



        #if mod(i,1) == 0
            println("Guardando el paso $i")
            SaveToParaview(MeshEC,[d1_EC], ["Cells"], "cell$i",
            dir; replace=replace,append=append)
            SaveToParaview(MeshVEGF,[d0_VEGF], ["VEGF"], "VEGF$i",
            dir; replace=replace,append=append)
        #end
    end
end


"""
SolverAngiogenesisAC
Solver of the angiogenic problem according to Anderson & Chaplain (1998) Method
"""
function SolverAngiogenesisAC(modelVEGF ::FemModel, modelFiber::FemModel,
    MeshEC, modelEC::FemModel, dt:: Float64, TotalTime::Float64,
    dir::String; replace=true, append=true)

    if replace && ispath(homedir()*dirSep*dir)
        rm(homedir()*dirSep*dir, recursive=true)
    end

    #To simulate initial diffusion of VEGF, EC concentration in every element must be 0
    d0_0=zeros(MeshEC.NofNodes)
    for element_i in modelVEGF.Elements
        element_i.Mat.ECs=d0_0[element_i.Nodes]
    end
    for element_i in modelFiber.Elements
        element_i.Mat.ECs=d0_0[element_i.Nodes]
    end

    #=
    d0_VEGF = SolveProgram(modelVEGF) #condicion inicial VEGF difundido por el tejido
    d0_Fiber = SolveProgram(modelFiber)
    d0_EC = uval(modelEC) #TODO cambiarlo para que sea de acuerdo a unas condiciones iniciales

    onlyPositiveValues!([d0_EC, d0_VEGF, d0_Fiber])
    =#

    d0_EC, d0_VEGF, d0_Fiber = initial_cond_AndersonChaplain(MeshEC)

    #se guarda la condicion inicial en Paraview
    println("Guardando el instante t = 0")
    SaveToParaview(MeshEC,[d0_EC], ["Cells"], "cell0",
    dir; replace=replace,append=append)
    SaveToParaview(MeshEC,[d0_VEGF], ["VEGF"], "VEGF0",
    dir; replace=replace,append=append)
    SaveToParaview(MeshEC,[d0_Fiber], ["Fibronectin"], "Fiber0",
    dir; replace=replace,append=append)


    DismantleMass!(modelEC)
    AssemblyMass!(modelEC)

    f_EC = Forces(modelEC)

    Nsteps = round(Int, TotalTime/dt)

    #d1_VEGF = d0_VEGF
    for i in 1:Nsteps

        #Updating VEGF and fibronectin distribution in EC model
        updateECmodel!(modelEC, d0_VEGF, d0_Fiber)

        Dismantle!(modelEC)
        Assembly!(modelEC)

        #Equation for EC distribution
        v0_EC = sparse(modelEC.M) \ (f_EC - sparse(modelEC.K)*d0_EC)

        d1_EC, v1_EC = dFormIntegration(modelEC, dt, f_EC, d0_EC, v0_EC, alpha = 1)
        d0_EC = d1_EC
        onlyPositiveValues!(d0_EC)

        # Calculates the result of the equations for the other species
        v0_VEGF, v0_Fiber = reactionEquationsAnderson(modelVEGF, modelFiber, d0_EC, d0_VEGF, d0_Fiber)

        # Integration for VEGF and Fibronectin substances
        d0_VEGF += v0_VEGF*dt
        d0_Fiber += v0_Fiber*dt


        onlyPositiveValues!([d0_VEGF, d0_Fiber])

        t = i*dt
        if mod(t,1) == 0
            println("Guardando el instante t = $t")
            t = round(Int64, t)
            SaveToParaview(MeshEC,[d1_EC], ["Cells"], "cell$t",
            dir; replace=replace,append=append)
            SaveToParaview(MeshEC,[d0_VEGF], ["VEGF"], "VEGF$t",
            dir; replace=replace,append=append)
            SaveToParaview(MeshEC,[d0_Fiber], ["Fibronectin"], "Fibronectin$i",
            dir; replace=replace,append=append)
        end
    end
end










"""
function updateVEGFmodel(modelVEGF, d0_VEGF, d0_EC)

Updates the VEGF model with the endothelial cell distribution in each element
Updates VEGF distribution according to the reaction equation

    •modelVEGF: FEM model made with VEGF material, containing reaction coefficients
    •d0_VEGF: VEGF distribution in the previous time step, intended to be updated following reaction equation:
    •d0_EC: Endothelial cell distribution
Returns:
    •modelVEGF: with updated EC distribution
    •d0_VEGF: updated VEGF distribution

"""

function updateVEGFmodel!(modelVEGF::FemModel, d0_VEGF::Vector{Float64}, d0_EC::Vector{Float64})

    #Model of VEGF consumption by ECs
    λ = modelVEGF.Elements[1].Mat.l_Cells  #VEGF consumption coefficient by ECs. Must be the same across the mesh
    d0_VEGF .-= λ*d0_EC.*d0_VEGF
    for VEGF_val in d0_VEGF
        if VEGF_val < 0.
            VEGF_val = 0.
        end
    end

    #Updating EC values in VEGF model
    for element_i in modelVEGF.Elements
        element_i.Mat.ECs = d0_EC[element_i.Nodes]
    end

    return modelVEGF, d0_VEGF
end



"""
function reactionEquationsAnderson

    Models Reaction equations for VEGF and Fibronectin substances according to model by
    Anderson & Chaplain (1998). These are:

    ∂c/∂t= -λ*n*c
    ∂f/∂t = ω*n - μ*n*f

    Where:  • c is the VEGF concentration value,
            • n is the EC density
            • f is the fibronectin concentration
            • λ is the VEGF uptake coefficient by ECs
            • ω is the production rate of fibronectin by ECs
            • μ is the uptake rate of fib. by ECs

    Returns:
        • v0_VEGF: vector containing ∂c/∂t in each node
        • v0_Fiber: vector containing ∂f/∂t in each node
"""

function reactionEquationsAnderson(modelVEGF::FemModel, modelFiber::FemModel,
    d0_EC::Vector{Float64}, d0_VEGF::Vector{Float64}, d0_Fiber::Vector{Float64})

    #=Model of VEGF uptake by ECs. Parameters involved must be the same
    across the whole mesh =#
    λ = modelVEGF.Elements[1].Mat.l_Cells  #VEGF uptake coefficient by ECs

    #Equation:
    v0_VEGF = -λ*d0_EC.*d0_VEGF

    #= Model of Fibronectin behaviour. Parameters involved must be the same
    across the whole mesh =#
    ω = modelFiber.Elements[1].Mat.w_Cells # production rate of fibronectin by ECs
    μ = modelFiber.Elements[1].Mat.μ_Cells # uptake rate of fib. by ECs

    #Equation:
    v0_Fiber = ω*d0_EC - μ*d0_EC.*d0_Fiber

    return v0_VEGF, v0_Fiber

 end


 """
 function updateECmodel(modelEC, d0_VEGF, d0_Fiber)

 Updates the EC model with the VEGF and fiber distribution in each element

     •modelEC: FEM model made with CellVegf2 material
     •d0_VEGF: VEGF distribution
     •d0_Fiber: Fibronectin distribution
 Returns:
     •modelEC: with updated VEGF and Fibronectin distribution

 """

 function updateECmodel!(modelEC::FemModel, d0_VEGF::Vector{Float64}, d0_Fiber::Vector{Float64})
      #Updating EC values in VEGF model
     for element_i in modelEC.Elements
         element_i.Mat.VEGF = d0_VEGF[element_i.Nodes]
         element_i.Mat.Fiber = d0_Fiber[element_i.Nodes]
     end

     return modelEC
 end


"""
function onlyPositiveValues!

    Substitutes by 0. those values in a Vector or Vector of Vectors that are  negative
"""

function onlyPositiveValues!(v::Vector{Float64})
    for i in 1:size(v)[1]
        if v[i] < 0.
            v[i] = 0.
        end
    end
    return v
end

function onlyPositiveValues!(A::Array{Array{Float64,1},1})
    for v in A
        onlyPositiveValues!(v)
    end
    return A
end



function EC_distrib_AndersonChaplain(x::Float64, y::Float64)
    f = exp(-x^2/0.001) *(sin(3*pi*y))^2
    return f
end

function VEGF_distrib_AndersonChaplain(x::Float64, y::Float64)
    f = exp(-(1-x)^2/0.45)
    return f
end

function Fiber_distrib_AndersonChaplain(x::Float64, y::Float64)
    f = 0.75*exp(-x^2/0.45)
    return f
end

"""
initial_cond_AndersonChaplain(Mesh::RegularMesh)

Returns the initial conditions for Anderson&Chaplain's simulation of angiogenesis

Returns d0_EC, d0_VEGF, d0_Fiber, where:
    · d0_EC is the initial distribution of Endothelial Cells
    · d0_VEGF is the initial distribution of VEGF
    · d0_Fiber is the initial distribution of fibronectin

"""
function initial_cond_AndersonChaplain(MeshEC::RegularMesh)
    n = MeshEC.NofNodes # number of nodes
    d0_EC = zeros(n)
    d0_VEGF = zeros(n)
    d0_Fiber = zeros(n)
    for node in 1:n
        coords = get_node_coordinates(node, MeshEC)
        d0_EC[node] = EC_distrib_AndersonChaplain(coords[1], coords[2])
        d0_VEGF[node] = VEGF_distrib_AndersonChaplain(coords[1], coords[2])
        d0_Fiber[node] = Fiber_distrib_AndersonChaplain(coords[1], coords[2])
    end

    return d0_EC, d0_VEGF, d0_Fiber
end
