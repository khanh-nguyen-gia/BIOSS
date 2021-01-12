using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

init_getBmat(_2D,Quad4, Elastic2D) # inicio de la funcion cinematica
initJacobian(_2D,Quad4)

N = 50 # elementos en x
M = 50 # elementos en y

println("***")
@time Mesh = RegularMesh(_2D, Quad4, Elastic2D, [N, M], [1., .2]) # dimesiones de la malla
println("***")

U1 = [0., 0.]
Bc1 = CreateBoundary(Mesh, Elastic2D, U1, BCLeft, BCDirichlet) # Dirichlet


Nodes = vec(Mesh.Nodes[end,:])
F1 = [0.1, 0] / length(Nodes)
Bc2 = CreateBoundary(Mesh, Elastic2D, Nodes, F1, BCForces) # Fuerza


Bc = [Bc1, Bc2]

p = FemModel("Elastic", "2D", Mesh, Bc)
@time u = SolveProgram(p)

function exportu(u,N)
    vars = []
    nombres = []

    for i in 1:N
        push!(vars, u[i:N:end])
        push!(nombres, "cosa $i")
    end
    return vars, nombres
end

vars, nombres = exportu(u,2)

SaveToParaview(Mesh,vars,nombres,"elastic","/Desktop/")

OutputDisp(Mesh,vars,nombres,"elastic2",homedir()*"/Desktop/")
