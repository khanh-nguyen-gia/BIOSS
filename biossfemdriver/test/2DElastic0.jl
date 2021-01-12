using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

init_getBmat(_2D,Quad4, Elastic2D) # inicio de la funcion cinematica
initJacobian(_2D,Quad4)

N = 20 # elementos en x
M = 20 # elementos en y

Nend = 100

println("***")
@time Mesh0 = RegularMesh(_2D, Quad4, Elastic2D, [N, M], [1., .2]) # dimesiones de la malla
println("***")

U1 = [0., 0.]
Bc1 = CreateBoundary(Mesh0, Elastic2D, U1, BCLeft, BCDirichlet) # Dirichlet

Bc3 = CreateBoundary(Mesh0, Elastic2D, U1, BCRight, BCDirichlet) # Dirichlet

NN = Int(N/2)
Nodes = vec(Mesh0.Nodes[end,:])

F0 = 2.6e-4
F1 = [0, F0*Nend]/length(Nodes)
Bc2 = CreateBoundary(Mesh0, Elastic2D, F1, Nodes, BCForces) # Fuerza

Bc = [Bc1,Bc2]

p = FemModel("Elastic", "2D", Mesh0, Bc)
@time u = SolveProgram(p)


u0 = deepcopy(u)
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

OutputDisp(Mesh0,vars,nombres,"elastic2",homedir()*"/Desktop/")


for i in 2:Nend

    global Bc1, Bc3, u, Mesh0, vars, Mesh, Nodes, F0, p, pp, N, M

    F1 = [0,F0*i]/length(Nodes)
    Bc2 = CreateBoundary(Mesh0, Elastic2D, F1, Nodes, BCForces) # Fuerza

    Bc = [Bc1, Bc2]

    Mesh = deepcopy(Mesh0)

    for i in 1:length(Mesh.Nodes)
        I = CartesianIndices(Mesh.Nodes)[i]
        Mesh.Coordinates[:,I] += [vars[1][i], vars[2][i]]
    end

    Mesh.Elements = BIOSSFemDriver.CreateElements(_2D, Quad4, Elastic2D, [N, M],Mesh.Nodes, Mesh.Coordinates)
    pp = FemModel("Elastic", "2D", Mesh, Bc)

    u = SolveProgram(pp)

    vars, nombres = exportu(u,2)

    OutputDisp(Mesh0,vars,nombres,"res/elast$i",homedir()*"/Desktop/")
end

using LinearAlgebra
println(norm(u0-u))

Output(Mesh,vars,nombres,"elastic",homedir()*"/Desktop/")
