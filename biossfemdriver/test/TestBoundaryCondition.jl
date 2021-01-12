using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

N = 400 # elementos

@time FemMesh = RegularMesh(_2D, Quad4, Difusion, [N, N], [1., 1.]) # dimesiones de la malla

Phi = [1.]
Nodes = vec(FemMesh.Nodes[:,1])
@time Bc1 = CreateBoundary(FemMesh, Difusion, Nodes, Phi, BCDirichlet)
