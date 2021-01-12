using Revise

using BIOSSFemCore
using BIOSSFemMaterial

using BIOSSFemElements
using BIOSSFemDriver


init_getNmat(_2D, Quad4, Difusion) # inicio de la funcion de forma
init_getBmat(_2D, Quad4, Difusion) # inicio de la funcion cinematica
initJacobian(_2D, Quad4)

N = 1000 # elementos en x
M = 200 # elementos en y

@time FemMesh = RegularMesh(_2D, Quad4, Difusion, [N, M], [1., 1.]) # dimesiones de la malla

Phi1 = [1.]
@time Bc1 = CreateBoundary(FemMesh, Difusion, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = [0.2]
@time Bc2 = CreateBoundary(FemMesh, Difusion, Phi2, BCRight, BCDirichlet)

F = [0.4]
@time Bc3 = CreateBoundary(FemMesh, Difusion, [(round(Int,N/4))*round(Int,M)+ round(Int,N/2)], F, BCForces)

Bc = [Bc1, Bc2]

model = FemModel("Difusion", "2D", FemMesh, Bc)
@time u = SolveProgram(model)

SaveToParaview(FemMesh,[u],["Calor"],"Test","/Desktop/")
