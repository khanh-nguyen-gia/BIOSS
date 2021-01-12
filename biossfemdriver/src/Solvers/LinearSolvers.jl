"""
    function SolveProgram(m :: FemModel)

    Solve a finite element model.

    Input:

        m :: FemModel            -> Finite element model

    Output:

        u :: Vector{Float64}  -> Displacement vector
"""
function SolveProgram(p::FemModel)

    # calculate global stifness matrix
    Dismantle!(p)
    Assembly!(p)

    # calculate forces
    f = Forces(p)

    # calculate initial displacements (given by boundary conditions)
    u = uval(p)

    # get auxiliar stifness matrix
    Kaux = getKll(p.K, p.FreeNodes)


    # solve the linear system
    faux = f[p.FreeNodes]
    uaux = Kaux \ faux
    u[p.FreeNodes] = uaux


    # return the displacements
    return u

end

function SolveProgramNonStationary(Mesh, model ::FemModel,
    dt:: Float64, TotalTime::Float64, dir::String; replace=true, append=true)

    d0 = uval(model)

    SolveProgramNonStationary(Mesh, model, dt, TotalTime, d0, dir; replace=true, append=true)
end

# Solver given initial conditions
function SolveProgramNonStationary(Mesh, model ::FemModel,
    dt:: Float64, TotalTime::Float64, d0::Array{Float64,1}, dir::String; replace=true, append=true)

    Assembly!(model)
    AssemblyMass!(model)

    #comprobation that the initial conditions are properly given
    if size(d0) != model.DOF
        println("Condiciones iniciales introducidas incorrectamente")
        d0=uval(model)
    end
    #
    # for Element_i in model.Elements
    #     Element_i.Mat.Vel = grad(d0_VEGF) # gradiente de vegf
    # end

    f = Forces(model)

    v0 = sparse(model.M) \ (f - sparse(model.K)*d0)

    Nsteps = round(Int, TotalTime/dt)

    if replace && ispath(homedir()*dirSep*dir)
        rm(homedir()*dirSep*dir, recursive=true)
    end

    for i in 1:Nsteps

        d1, v1 = dFormIntegration(model, dt, f, d0, v0, alpha = 1)
        v0 = v1; d0 = d1

        if mod(i,5) == 0
            println("Guardando el paso $i")
            SaveToParaview(Mesh,[d1], ["Temp"], "temp$i",
            dir; replace=replace,append=append)
        end
    end

end


# function SolveProgram2(p::FemModel)
#
#     @time fastRegularAssembly!(p)
#
#     @time f = Forces(p)
#
#     u = uval(p)
#
#     Kaux = getKll(p.K, p.FreeNodes)
#
#     faux = f[p.FreeNodes]
#     @time uaux = Kaux \ faux
#     u[p.FreeNodes] = uaux
#
#     return u
#
# end

function SolveProgram3(p::FemModel)

    @time fastRegularAssembly!(p)

    @time f = Forces(p)

    u = uval(p)

    Kaux = getKll(p.K, p.FreeNodes)

    faux = f[p.FreeNodes]
    @time uaux = Kaux \ faux
    u[p.FreeNodes] = uaux

    return u

end
#
export SolveProgram3

# function SolveProgram2(p::FemModel,N,i,j,v)
#
#     Assembly!(p)
#
#     f = Forces(p)
#     u = uval(p)
#
#     f[(j-1)*(N+1) + i] += v
#     f[(j-1)*(N+1) + i+1] += v
#     f[(j)*(N+1) + i] += v
#     f[(j)*(N+1) + i+1] += v
#
#     Kaux = getKll(p.K, p.FreeNodes)
#
#     faux = f[p.FreeNodes]
#     uaux = Kaux \ faux
#     u[p.FreeNodes] = uaux
#
#     return u
#
# end


using LinearAlgebra
function CalcGradient(p::FemModel,u)
    C = [-1. -1.;
          1. -1.;
          1.  1.;
         -1   1.]

    du = zeros(Float64,3,p.DOF)

    for el in p.Elements
        for i in 1:4
            Jac = element_jacobian(el,C[i,:])
                # Calculate element jacobian
            J = det(Jac)
                    # Inverse of the jacobian, to calculate B matrix
            J_1 = inv(Jac)
                    # Calc B matrix
            Bmat = elements.getBmat(el,J_1,C[i,:])
            du[1:2,el.Nodes[i]] = -u[el.Nodes]'*Bmat'
        end
    end

    return du

end


function CalcGradient2(p::FemModel,u)
    C, code = get_basis_props(p.Elements[1])

    El_gdl = get_basis_gdl(p.Elements[1])
    dim = length(Size(u))
    du = zeros(length(Size(u)),p.DOF)

    for el in p.Elements
        for i in 1:El_gdl
            Jac = ElementJacobian(el,C[i,:])
                # Calculate element jacobian
            J = det(Jac)
                    # Inverse of the jacobian, to calculate B matrix
            J_1 = inv(Jac)
                    # Calc B matrix
            Bmat = getBmat(el,J_1,C[i,:])
            du[1:dim,el.Nodes[i]] = -u[el.Nodes]'*Bmat'
        end
    end

    return du

end

export CalcGradient2
