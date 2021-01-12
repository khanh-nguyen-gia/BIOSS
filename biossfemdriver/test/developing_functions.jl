"""
File for function development before including them in BIOSSFemDriver
"""

#TODO: Hacer el equivalente para el espacio 8-Adyacente
function init_VesselInfo(d0_EC::Vector{Float64})

    d_vessels = -1 .* ones(Int64, length(d0_EC))
    notch_viability = trues(length(d0_EC))
    active_nodes = zeros(Int64, length(d_vessels)) #active_nodes[i] = 0 if inactive, 1 if active
    nodes_BDs = zeros(Float64, length(d_vessels)) #Contains the birthdate of the vessel segment to which each node belongs
    fill!(nodes_BDs, Inf64)

    for i in 1:length(d0_EC)
        if d0_EC[i]>= 0.8
            d_vessels[i]= 0
            active_nodes[i] = 1
            nodes_BDs[i] = 0.
        end
    end
    return VesselInfo(d_vessels, active_nodes, nodes_BDs, notch_viability)
end

function Update_VesselInfo!(d0_EC::Vector{Float64}, d_vessels::Vector{Int64},
    active_nodes::Vector{Int64}, nodes_BDs::Vector{Float64},
    Adjacent_nodes::Dict{Int64,Array{Int64,1}})

    #=
    #first, find and deactivate nodes that are no longer in a vessel
    for i in 1:length(d0_EC)
        if (d0_EC[i]<0.8) & (d_vessels[i]!=-1) #It was part of a vessel, but now there is not enough ρ
            d_vessels[i] = -1
            active_nodes[i] = 0
            nodes_BDs[i] = Inf64
        end
    end
    =#

    #then, activate those that are just now part of a vessel
    for i in 1:length(d0_EC)
        if (d0_EC[i]>=0.8) & (d_vessels[i]==-1) #It hs enough ρ, but is not oficially part of any vessel
            adj = Adjacent_nodes[i]
            invessels_adj = filter(j-> d_vessels[j]!=-1, adj) #adjacent nodes belonging to a vessel
            #active_adj = filter(j-> active_nodes[j]==1, invessels_adj)
            if !Base.isempty(invessels_adj) #if there is a vessel nearby
                node = rand[invessels_adj] #randomly choose one of the nodes
                d_vessels[i] = d_vessels[node]
                active_nodes[i] = active_nodes[node]
                nodes_BDs[i] = nodes_BDs[node]
            end
        end
    end

    return d_vessels, active_nodes, nodes_BDs
end

function Update_VesselInfo!(d0_EC::Vector{Float64}, vessels::VesselInfo,
    Adjacent_nodes::Dict{Int64,Array{Int64,1}})

    d_vessels = vessels.distribution
    active_nodes = vessels.active_nodes
    nodes_BDs = vessels.nodes_BDs

    #=
    #first, find and deactivate nodes that are no longer in a vessel
    for i in 1:length(d0_EC)

        if (d0_EC[i]<0.8) & (d_vessels[i]!=-1) #It was part of a vessel, but now there is not enough ρ
            d_vessels[i] = -1
            active_nodes[i] = 0
            #nodes_BDs[i] = Inf64
        end
    end
    =#

    #then, activate those that are just now part of a vessel
    for i in 1:length(d0_EC)
        if (d0_EC[i]>=0.8) & (d_vessels[i]==-1) #It hs enough ρ, but is not oficially part of any vessel
            adj = Adjacent_nodes[i]
            invessels_adj = filter(j-> d_vessels[j]!=-1, adj) #adjacent nodes belonging to a vessel
            #active_adj = filter(j-> active_nodes[j]==1, invessels_adj)
            if !Base.isempty(invessels_adj) #if there is a vessel nearby
                node = rand(invessels_adj) #randomly choose one of the nodes
                d_vessels[i] = d_vessels[node]
                active_nodes[i] = active_nodes[node]
                if nodes_BDs[i] == Inf64
                    nodes_BDs[i] = nodes_BDs[node]
                end
            end
        end
    end

    return vessels
end

#=
function reaction_TAF(d_EC::Vector{Float64}, nodes_tumour::Vector{Int64})
    nnodes = length(d_EC)
    P = 350.0
    Uu = 21.875
    Ud = 0.35
    f_hyc = 1.

    K_reaction = zeros(nnodes, nnodes)
    forces_reaction = zeros(nnodes)
    for i in nodes_tumour
        K_reaction[i,i] += -P
        forces[i] = -P * f_hyc
    end

    for i in 1:nnodes
        if d_EC[i] >= 0.
            U = Uu*d_EC[i]
        else
            U = -Ud*d_EC[i]
        end

        K_reaction[i] -= U
    end

    return K_reaction, forces_reaction
end


function TAF_forces(el::AbstractElement{D,B,M},
    Point::Union{Vector{Float64},Float64}) where {D <: AbstractDim,
    B <: AbstractBasis, M <: TAF}

    Jac = ElementJacobian(el,Point)

    J = det(Jac)
    # Inverse of the jacobian, to calculate B matrix
    J_1 = inv(Jac)
    # Calc B matrix
    Bmat = getBmat(el, J_1, Point)
    Nmat = getNmat(el, Point)
    hycs = Nmat*el.Mat.HyCs #Concentracion de celulas hipoxicas en el punto

    if hycs >= 0.5
        P = el.Mat.Prod_hyc
    else
        P = 0.
    end
    f_i = -Nmat'*(P)*J
    return f_i
end

function calculo_f_TAF(el::Element{D,B,M}) where {D <: AbstractDim, B<: AbstractBasis, M<:TAF}

    # get IntPoints

    (IntPoints, Weights) = get_int_points(el)

    #number of integration points

    Npoints = size(IntPoints,1)

    N = get_basis_gdl(el)
    f = zeros(N*el.ndof)

    for i in 1:Npoints
        Point = IntPoints[i, :]
        Jac = ElementJacobian(el,Point)
        J = det(Jac)
        # Inverse of the jacobian, to calculate B matrix
        J_1 = inv(Jac)
        # Calc B matrix
        Bmat = getBmat(el, J_1, Point)
        Nmat = getNmat(el, Point)

        hycs = Nmat*el.Mat.HyCs #Concentracion de celulas hipoxicas en el punto
        if hycs >= 0.5
            P = el.Mat.Prod_hyc
        else
            P = 0.
        end

        f_i = -Nmat'*(P)*J
        # Evaluate K at point of integration
        f += f_i*Weights[i]
    end

    f = - K * u

    return f
end

=#
#import BIOSSFemDriver:
