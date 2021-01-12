using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

clearconsole()


using Plots

function SolveNonLinear2(N, Niters)
    # elementos en x

    initJacobian(_1D, Lin2)
    init_getBmat(_1D, Lin2, TestNonElastic1D)

    u = zeros(N+1)
    u0 = copy(u)


    UU = zeros(Niters, N+1)

    F = collect(range(0.,25e3,length=Niters))

    for n in 1:Niters

        Mesh = RegularMesh(_1D, Lin2, TestNonElastic1D, [N], [1.]) # dimesiones de la malla

        for el in Mesh.Elements
            el.Mat.DMat .= 1e9*1e-4
        end

        U1 = [0.]
        Bc1 = CreateBoundary(Mesh, TestNonElastic1D, U1, [1], BCDirichlet) # Dirichlet

        F1 = [F[n]]

        Bc2 = CreateBoundary(Mesh, TestNonElastic1D, F1, [N+1], BCForces)
        Bc = [Bc1,Bc2]


        p = FemModel("nonlinear", "1D", Mesh, Bc)

        u = uval(p)

        for i in 1:length(p.Elements)
            el = p.Elements[i]
            nodes = SpatialNodesToIndices(el)
            BIOSSFemDriver.add!(p.K,nodes,nodes,LocalStiffness(el, u0[nodes]))
        end

        println("Iteración $n")

        # u[p.FreeNodes], it = SolveNewtonR(u -> Residual2(
        # p, TestNonElastic1D, u), u0[p.FreeNodes], tol = 1e-5)

        u[p.FreeNodes], it = SolveNewtonR(u -> Residual2(
        p, TestNonElastic1D, u), u0[p.FreeNodes], tol = 1e-3)

        # println(getKll(p.K,p.FreeNodes)*u[p.FreeNodes])

        # println(Forces(p))

        # println(sparse(p.K)*u)

        # println("Residuo", Residual2(p, TestNonElastic1D, u[p.FreeNodes]))
        # println(u)
        # println(it)

        u0 = copy(u)

        UU[n,:] = u

    end

    return u, UU

end

init_getBmat(_1D, Lin2, TestNonElastic1D)
initJacobian(_1D,Lin2)

function Residual2(p,::Type{TestNonElastic1D}, u)

    f = Forces(p)

    up = zeros(p.DOF)
    up[p.FreeNodes] = u

    Res = zeros(length(up))

    for i in 1:length(p.Elements)
        el = p.Elements[i]
        nodes = SpatialNodesToIndices(el)
        BIOSSFemDriver.add!(p.K,nodes,nodes,LocalStiffness(el, up[nodes]))
    end

    Res = f[p.FreeNodes] - getKll(p.K,p.FreeNodes)*u

    return Res
end


function Residual3(p,::Type{TestNonElastic1D}, u)

    f = Forces(p)

    up = zeros(p.DOF)
    up[p.FreeNodes] = u

    Res = f[p.FreeNodes] - getKll(p.K,p.FreeNodes)*u

    return Res
end

N = 3
Niters = 50
@time u, UU = SolveNonLinear2(N, Niters)
println(u)
F = collect(range(0.,1 , length=size(UU,1)))
pp = plot(legend=nothing)
for i in 1:N+1
    global pp
    plot!(pp,UU[:,i], F, label="u: Nodo $i")
end

display(pp)
