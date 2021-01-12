using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

init_getBmat(_3D, Brick8, Elastic3D) # inicio de la funcion cinematica
initJacobian(_3D, Brick8)

N = 10 # elementos en x
M = 10 # elementos en y
O = 10

println("***")
@time Mesh = RegularMesh(_3D, Brick8, Elastic3D, [N, M, O], [2., 0.5, 4.]) # dimesiones de la malla
println("***")

U1 = [0., 0., 0.]
@time Bc1 = CreateBoundary(Mesh, Elastic3D, vec(Mesh.Nodes[:,:,1]), U1, BCDirichlet) # Dirichlet

F1 = [0., 0.02, 0.]
nodes = vec(Mesh.Nodes[end,:,Int(O/2):end])
@time Bc2 = CreateBoundary(Mesh, Elastic3D, nodes, F1./length(nodes),BCForces) # Fuerzas

nodes = vec(Mesh.Nodes[1,:,Int(O/2):end])
@time Bc3 = CreateBoundary(Mesh, Elastic3D, nodes, -F1./length(nodes), BCForces) # Fuerzas

F1 = [0., -0, 0.08]
nodes = vec(Mesh.Nodes[:,:,end])
@time Bc4 = CreateBoundary(Mesh, Elastic3D, nodes, F1./length(nodes), BCForces) # Fuerzas

Bc = [Bc1, Bc2, Bc3, Bc4]

p = FemModel("Elastic", "3D", Mesh, Bc)
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

vars, nombres = exportu(u,3)

SaveToParaview(Mesh,vars,nombres,"elastic","/Desktop")
#OutputDisp3D(Mesh,vars,nombres,"elastic3D",homedir()*"/Desktop")
