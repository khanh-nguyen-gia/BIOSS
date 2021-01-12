using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

init_getNmat(_2D,Quad4, ConvecDif) # inicio de la funcion de forma
init_getBmat(_2D,Quad4, ConvecDif) # inicio de la funcion cinematica
initJacobian(_2D,Quad4)

N = 20 # elementos en x
M = 20 # elementos en y

Mesh = RegularMesh(_2D, Quad4, ConvecDif, [N, M], [1., 1.]) # dimesiones de la malla

for el in Mesh.Elements
    el.Mat.Vel = [cos(-pi/3),sin(-pi/3)]
    el.Mat.DMat .*= 1e-3
end

Phi1 = [0.]
nodes =  vcat(vec(Mesh.Nodes[:,1]),
vec(Mesh.Nodes[1,2:round(Int,end/2)]),
vec(Mesh.Nodes[end,:]))
# nodes = vec(vec(Mesh.Nodes[1,2:140]))
Bc1 = CreateBoundary(Mesh, ConvecDif, Phi1, nodes, BCDirichlet) # Dirichlet

Phi2 = [1.]
nodes = vcat(vec(Mesh.Nodes[:,end]),
vec(Mesh.Nodes[1,1:end-1]))
Bc2 = CreateBoundary(Mesh, ConvecDif, Phi2, nodes, BCDirichlet)

Bc = [Bc1,Bc2]

p = FemModel("Convection", "2D", Mesh, Bc)

@time u = SolveProgram(p)

SaveToParaview(Mesh,[u],["cosas"],"datos","/Desktop")
