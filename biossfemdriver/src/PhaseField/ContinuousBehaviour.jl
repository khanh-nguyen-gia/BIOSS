# Collection of functions used for simulation of angiogenesis using a Phase-Field model

function EC_distrib_PhaseField(x::Float64, y::Float64)
    if x <= 8
        f = 1.
    elseif x >= 300-8
        f = 1.
    else
        f = -1.
    end
    return f
end

function VEGF_distrib_PhaseField(x::Float64, y::Float64)
    r = sqrt((x-1)^2+(y-0.5)^2)
    v = (sqrt(5)-0.1)/(sqrt(5)-1.)
    if r <= 0.1
        f = 1.
    else
        f = (v-r)^2/(v-0.1)^2
    end
    #f = exp(-(1-x)^2/0.45)
    return f
end

function HyC_distrib_PhaseField(MeshEC::RegularMesh, modelTAF::FemModel)
    x0 = MeshEC.EndPoint[1]/2
    y0 = MeshEC.EndPoint[2]/2

    #Hypoxic cell distribution:
    HyC_number = 20 #number of HyCs

    #Follows a Normal Distribution regarding distance to the middle of the mesh
    #95% of these cells will be within R = 1/4 D from the middle
    D = minimum(MeshEC.EndPoint) #Shortest dimension of the mesh
    R = 0.15*D
    Zref= 1.967 #97.5% of the times z<Zref; z = (r-μ)/σ
    μ = 0.
    σ = (R-μ)/Zref
    r_distr = Normal(μ, σ)
    r_distr = truncated(r_distr, 0.0, Inf) #Only positive values

    #Distribution of θ in polar coordinates to give full positions
    angle_distr = Uniform(0., 2*pi)

    r = rand(r_distr, HyC_number)
    θ = rand(angle_distr, HyC_number)

    if length(MeshEC.DimNodes) == 3 #3D mesh

        z0 = MeshEC.EndPoint[3]/2
        ϕ = rand(Uniform(-pi/2, pi/2), HyC_number)
        HyC_vec = [HypoxCell(i, [ x0+r[i]*cos(ϕ[i])*cos(θ[i]), y0+r[i]*cos(ϕ[i])*sin(θ[i]), z0+r[i]*sin(ϕ[i]) ], MeshEC) for i in 1:HyC_number]

    else
        HyC_vec = [HypoxCell(i, [ x0+r[i]*cos(θ[i]), y0+r[i]*sin(θ[i]) ], MeshEC) for i in 1:HyC_number]
    end

    return HyC_vec
end

function Fiber_distrib_PhaseField(x::Float64, y::Float64)
    f = 0.75*exp(-x^2/0.45)
    return f
end

"""
initial_cond_PhaseField(Mesh::RegularMesh)

Returns the initial conditions for a simulation of angiogenesis using a Phase-Field Model

Returns d0_EC, d0_VEGF, d0_Fiber, where:
    · d0_EC is the initial distribution of Endothelial Cells
    · d0_VEGF is the initial distribution of VEGF
    · d0_Fiber is the initial distribution of fibronectin

"""
function initial_cond_PhaseField(MeshEC::RegularMesh, modelTAF::FemModel)

    n = MeshEC.NofNodes # number of nodes

    d0_EC = zeros(n)
    for node in 1:n
        coords = get_node_coordinates(node, MeshEC)
        d0_EC[node] = EC_distrib_PhaseField(coords[1], coords[2])
    end

    HyC_vec = HyC_distrib_PhaseField(MeshEC, modelTAF)


    AssemblyMass!(modelTAF)
    UpdateTAFmodel!(modelTAF, d0_EC, HyC_vec)
    d0_VEGF = zeros(n)#SolveProgram(modelTAF)

    #In order to reduce initial errors propagation
    TAF_hyc = modelTAF.Elements[1].Mat.TAF_hyc
    nodesinsideHyCs = nodes_HyCs(modelTAF)
    d0_VEGF[nodesinsideHyCs] .= TAF_hyc

    #=
    for node in 1:n
        coords = get_node_coordinates(node, MeshEC)
        d0_EC[node] = EC_distrib_PhaseField(coords[1], coords[2])
        d0_VEGF[node] = VEGF_distrib_PhaseField(coords[1]/MeshEC.EndPoint[1], coords[2]/MeshEC.EndPoint[2])
        d0_Fiber[node] = Fiber_distrib_PhaseField(coords[1]/MeshEC.EndPoint[1], coords[2]/MeshEC.EndPoint[2])
    end
    =#
    return d0_EC, d0_VEGF, HyC_vec, modelTAF
end

"""
    rho_distrib_TC(r_rel)

    Returns EC density distribution inside a TipCell, given the relative distance
    to its center (r_rel = distance/radius)
"""
function rho_distrib_TC(r_rel::Float64)
    genericTC = TipCell(0, _2D)
    x = r_rel

    x0 = 0.4
    x´ = x-x0 #relative distance to internal radius
    y0 = genericTC.rho

    x1´= 1. - x0
    y1 = 0.8

    c = y0
    b= 0. #dy/dx = 0. cuando x-x0 es 0.
    a = (y1-c-b*x1´)/(x1´^2) # y = y1 cuando x = 1

    if x <= x0
        y = y0
    else
        y = a*(x´)^2 + b*x´ + c
    end

    return y
end

"""
function FreeEnergyContribution(modelEC, f_EC, d_EC):
    • modelEC: FemMaodel made with Endothelium material
    • f_EC: Forcing term of the differential equation
    • d_EC: Endothelial cell distribution

Returns the contribution of the surface free energy functional (capillary-extracellular matrix interface) to the time derivative of EC density.
 This follows a diffusive behaviour:

    ∂ρ/∂t = M ⋅ λ² ⋅ (Δρ)

Where M and λ are parameters of the Endothelium material:
     • M: Motility coefficient
     • λ: Capillary wall width coefficient
"""
function FreeEnergyContribution(modelEC::FemModel, f_EC::Vector{Float64}, d_EC::Vector{Float64})
    return  sparse(modelEC.M) \(f_EC - sparse(modelEC.K)*d_EC)
end


"""
 ChemicalContribution

Returns the contribution of chemical free energy functional to the time derivative of EC density.
This determines wether a capillary grows or shrinks depending on theconditions.
Modelled by the following equation in each node:

    ∂ρ/∂t = - M ⋅ μ(ρ, f);
    μ = 0.5⋅(ρ²-1)⋅(ρ-3⋅α⋅γ)

Where:
    • ρ: EC density value in the point of study
    • f: TAF conentration value in the point of study
    • M: Motility coefficient
    • α: Phenotype switch parameter 1 of Endothelium material
    • γ: Tilting function, towards growth or regression depending on f value


--------------------------------------------------------------------------
ChemicalContribution(Mat, ρ, f):
        • Mat: Endothelium material containing the parameters involved
        • ρ: EC density value in the point of study
        • f: TAF conentration value in the point of study
"""
function ChemicalContribution(Mat::Endothelium, ρ::Float64, f::Float64)
    α = Mat.α_Cell
    M = Mat.M_Cell
    γ = tilting_function_TAF(Mat, f)

    μ = 0.5*(ρ^2-1)*(ρ-3*α*γ)

    return (-M) .* μ
end


function ChemicalContribution(Mat::Endothelium, d_EC::Vector{Float64}, d_VEGF::Vector{Float64})

    μ = zeros(size(d_EC)[1])
    for i in 1:size(d_EC)[1]
    μ[i] = ChemicalContribution(Mat, d_EC[i], d_VEGF[i])
    end

    return μ
end

""""
ChemicalContribution(ModelEC, d_EC, d_VEGF):
        • ModelEC: FemModel made with Endothelium material
        • d_EC: Endothelial cell distribution
        • d_VEGF: TAF distribution
"""
function ChemicalContribution(ModelEC::FemModel, d_EC::Vector{Float64}, d_VEGF::Vector{Float64})
    Mat = ModelEC.Elements[1].Mat

    return ChemicalContribution(Mat, d_EC, d_VEGF)
end



"""
tilting_function_TAF

Function that tilts capillary behaviour towards growth or regression depending on TAF concentration.

    γ(f) = exp(-exp(β*(f-f_act))) - exp(-1)

---------------------------------------------------------------------------
tilting_function_TAF(Mat, f):
    • Mat: Endothelium material containing the parameters involved
    • f: TAF conentration value in the point of study

"""
function tilting_function_TAF(Mat::Endothelium, f::Float64)
    β = Mat.β_Cell
    f_act = Mat.f_act
    γ = exp(-exp(β*(f-f_act))) - exp(-1)

    return γ
end

function tilting_function_TAF(Mat::Endothelium, f::Vector{Float64})
    γ = zeros(size(f)[1])
    for i in 1:size(f)[1]
        γ[i] = tilting_function_TAF(Mat, f[i])
    end

    return γ
end


"""
function time_derivative_Endothelium(ModelE, f_EC, d_EC, d_VEGF)

    Returns the time derivative of EC distribution
Input:

        • ModelEC: FemModel made with Endothelium material
        • f_EC: Forcing term of the differential equation
        • d_EC: Endothelial cell distribution
        • d_VEGF: TAF distribution
"""
function time_derivative_Endothelium(ModelEC::FemModel, f_EC::Vector{Float64},
    d_EC::Vector{Float64}, d_VEGF::Vector{Float64})

    v1 = FreeEnergyContribution(ModelEC, f_EC, d_EC)
    v2 = ChemicalContribution(ModelEC, d_EC,d_VEGF)

    return v1 + v2
end


"""
function UpdateTAFmodel!(modelTAF, d_EC, HyC_vector)

    Updates modelTAF with the info about Endothelial and Hypoxic Cell distributions
"""
function UpdateTAFmodel!(modelTAF::FemModel, d_EC::Vector{Float64}, HyC_vector::Vector{HypoxCell})
    Dismantle!(modelTAF)

    active_HyCs = filter(i -> i.active, HyC_vector)
    active_nodes = empty([0])
    for i in active_HyCs
        active_nodes = vcat(active_nodes, i.Nodes)
    end

    for el in modelTAF.Elements
        for i in 1:length(el.Nodes)
            node = el.Nodes[i]
            #HyC info update for the node
            if Base.isempty(setdiff(node, active_nodes)) #the node is in an active hypoxic cell
                el.Mat.HyCs[i] = 1
            else
                el.Mat.HyCs[i] = 0
            end

            #rho  info update for the node
            el.Mat.ECs[i] = d_EC[node]
        end
    end
    Assembly!(modelTAF)
    return modelTAF
end


"""
    nodes_HyCs(modelTAF)

    Returns the list of nodes in whic there are Hypoxic Cells, according to the
    info inside the model
"""
function nodes_HyCs(model::FemModel)
    nodes_HyCs = empty([0])
    for el in model.Elements
        local_nodes = SpatialNodesToIndices(el)
        local_HyCs = el.Mat.HyCs
        for i in 1:length(el.Nodes)
            if local_HyCs[i] == 1
                push!(nodes_HyCs, local_nodes[i])
            end
        end
    end
    nodes_no_HyCs = setdiff(1:model.DOF, nodes_HyCs)
    nodes_HyCs = setdiff(1:model.DOF, nodes_no_HyCs)

    return nodes_HyCs
end

"""
    nodes_ECs(modelTAF)

    Returns the list of nodes in whic there are Endothelial Cells, according to the
    info inside the model (EC density >= 0.)
"""
function nodes_ECs(model::FemModel)
    nodes_ECs = empty([0])
    for el in model.Elements
        local_nodes = SpatialNodesToIndices(el)
        local_ECs = el.Mat.ECs
        for i in 1:length(el.Nodes)
            if local_ECs[i] >= 0.
                push!(nodes_ECs, local_nodes[i])
            end
        end
    end
    nodes_no_ECs = setdiff(1:model.DOF, nodes_ECs)
    nodes_ECs = setdiff(1:model.DOF, nodes_no_ECs)
    return nodes_ECs
end

"""
    extract_EC_distribution(modelTAF)

    Returns the Endothelial Cell density distribution from the info inside the model (loaded with UpdateTAFmodel!)
"""
function extract_EC_distribution(model::FemModel)
    d_ECs =  spzeros(Float64, model.DOF)

    for el in model.Elements
        local_nodes = SpatialNodesToIndices(el)
        local_ECs = el.Mat.ECs
        for i in 1:length(el.Nodes)
            d_ECs[local_nodes[i]] = local_ECs[i]
        end
    end
    return d_ECs
end

"""
    function Update_TAF_evolution!

    Calculates TAF evolution after tmax has passed. Updates modelTAF and d0_TAF after that

"""
function Update_TAF_evolution!(modelTAF::FemModel, d0_TAF::Vector{Float64},d0_EC::Vector{Float64}, HyC_vec::Vector{HypoxCell}, tmax::Float64, dt::Float64)

    UpdateTAFmodel!(modelTAF, d0_EC, HyC_vec)
    f_TAF = Forces(modelTAF)
    #ratio = modelTAF.Elements[1].Mat.ratio
    ratio = 1. #using a time ratio may be a mistake

    for i in 1 : Int64(tmax*ratio÷dt)
        t = i*dt

        #Dismantle!(modelTAF)
        #Assembly!(modelTAF)

        #Backward Euler
        if i == 1
            global v_TAF
            v_TAF = sparse(modelTAF.M) \ (f_TAF - sparse(modelTAF.K)*d0_TAF) .* ratio
        end
        d0_TAF, v_TAF = dFormIntegration(modelTAF, dt, f_TAF, d0_TAF, v_TAF, alpha = 1)


        #Forward Euler
        #v_TAF = sparse(modelTAF.M) \(f_TAF - sparse(modelTAF.K)*d0_TAF)
        #d0_TAF += v_TAF*dt

        #Filter negative and values above the limit to avoid divergence of errors
        TAF_hyc = modelTAF.Elements[1].Mat.TAF_hyc
        for j in 1:length(d0_TAF)
            if d0_TAF[j] > TAF_hyc
                d0_TAF[j] = TAF_hyc
                #println("TAF concentration exceeds limit value")
            elseif d0_TAF[j] < 0.
                d0_TAF[j] = 0.
            end
        end

        #=
        i_save = 15
        if (i>=i_save)&(i%i_save==0)
            println("Operando el instante t = $t")

            SaveToParaview(Mesh,[d0_TAF], ["TAF"], "TAF$(i÷i_save)", dir; replace=true,append=true)
        end=#
    end
    return modelTAF, d0_TAF
end
