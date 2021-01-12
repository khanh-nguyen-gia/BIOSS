using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

using SparseArrays


init_getNmat(_2D,Quad4, ConvecDif) # inicio de la funcion de forma
init_getBmat(_2D,Quad4, ConvecDif) # inicio de la funcion cinematica
initJacobian(_2D,Quad4)

N = 20 # elementos en x
M = 20 # elementos en y

Mesh = RegularMesh(_2D, Quad4, ConvecDif, [N, M], [1., 1.]) # dimesiones de la malla

for i in 1:Int(length(Mesh.Elements))
    Mesh.Elements[i].Mat.DMat .*= .5
end

for i in 1:Int(length(Mesh.Elements)/2)
    Mesh.Elements[i].Mat.Vel = [3,.3]*0
    Mesh.Elements[i+Int(length(Mesh.Elements)/2)].Mat.Vel = [1,.0]
    Mesh.Elements[i].Mat.DMat .*= 1
end

Phi1 = [0.]
Bc1 = CreateBoundary(Mesh, ConvecDif, Phi1, BCRight, BCDirichlet) # Dirichlet

Phi2 = [5.]
nodes = vec(Mesh.Nodes[1,round(Int,end/2-1):end])
Bc2 = CreateBoundary(Mesh, ConvecDif, Phi2, nodes, BCDirichlet)

Phi3 = [10.]
nodes = vec(Mesh.Nodes[1,1:round(Int,end/2-1)])
Bc3 = CreateBoundary(Mesh, ConvecDif, Phi3, nodes, BCDirichlet)

Phi4 = [10.]
nodes = vec(Mesh.Nodes[end,round(Int,end/2-1):end])
Bc4 = CreateBoundary(Mesh, ConvecDif, Phi4, nodes, BCDirichlet)

Bc5 = CreateBoundary(Mesh, Heat, [1.], [70], BCForces)

Bc = [Bc1, Bc2, Bc3, Bc4, Bc5]


model = FemModel("Convection", "2D", Mesh, Bc)

TotalTime = 2.
dt = 0.005

dir = "Desktop/Results/"

SolveProgramNonStationary(Mesh, model, dt, TotalTime, dir)
