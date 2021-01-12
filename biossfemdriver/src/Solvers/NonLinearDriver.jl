function SolveNonLinear(p::FemModel)

    N = 3 # elementos en x

    Mesh = RegularMesh(_1D, Lin2, TestNonElastic1D, [N], [1.]) # dimesiones de la malla

    for el in Mesh.Elements
        el.Mat.DMat .= 1e9*1e-4
    end

    initJacobian(_1D, Lin2)
    init_getBmat(_1D, Lin2, TestNonElastic1D)

    u = uval(p)
    u0 = copy(u)

    Niters = 50

    U1 = [0.]
    Bc1 = CreateBoundary(Mesh, TestNonElastic1D, U1, [1], BCDirichlet) # Dirichlet

    F = collect(range(0.1,25e3 , length=Niters))

    for n in 1:Niters

        F1 = [F[n]]
        Bc2 = CreateBoundary(Mesh, TestNonElastic1D, F1, [N+1], BCForces)
        Bc = [Bc1,Bc2]
        p = FemModel("nonlinear", "1D", Mesh, Bc)

        println("Iteración $n")

        for i in 1:length(p.Elements)
            el = p.Elements[i]
            nodes = SpatialNodesToIndices(el)
            add!(p.K,nodes,nodes,LocalStiffness(el, u[nodes]))
        end

        u[p.FreeNodes], it = SolveNewtonR(u -> Residual(p, TestNonElastic1D, u), u0[p.FreeNodes], tol = 1e-5)

        println(Forces(p))
        println("Residuo", Residual(p, TestNonElastic1D, u[p.FreeNodes]))
        println(u)
        println(it)


        u0 = copy(u)

    end

    return u, p

end


function Residual(p,::Type{TestNonElastic1D}, u)

    f = Forces(p)

    Res = copy(u)

    Res= f[p.FreeNodes] - getKll(p.K,p.FreeNodes)*u

    return Res
end

# function NlForces(model::FemModel, u0)
#     f = spzeros(Float64, model.DOF)
#
#     for B in model.Bc
#         if B.Kind == BCDirichlet()
#             for i = 1:length(B.Elements)
#                 local_nodes = SpatialNodesToIndices(model.Elements[B.Elements[i]])
#                 f[local_nodes] += LocalForces(model.Elements[B.Elements[i]], B.Cond[i], u0[local_nodes])
#             end
#         elseif B.Kind == BCForces()
#             for i = 1:length(B.Elements)
#                 local_nodes = SpatialNodesToIndices(model.Elements[B.Elements[i]])
#                 f[local_nodes] += B.Cond[i]
#             end
#         end
#     end
#     return Array(f)
# end
