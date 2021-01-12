using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

using SparseArrays

ClearStoredCache()
init_getNmat(_2D,Quad4, CellVegf2) # inicio de la funcion de forma
init_getNmat(_2D,Quad4, Difusion)

init_getBmat(_2D,Quad4, CellVegf2) # inicio de la funcion cinematica
init_getBmat(_2D,Quad4, Difusion)

initJacobian(_2D,Quad4)

N = 50 # elementos en x
M = 50 # elementos en y

MeshVEGF = RegularMesh(_2D, Quad4, Difusion, [N, M], [20., 20.]) # Malla comportamiento VEGF

MeshEC = RegularMesh(_2D, Quad4, CellVegf2, [N, M], [20., 20.]) # Malla comportamiento celulas endoteliales

Phi1 = [0.5]
BcVEGF1 = CreateBoundary(MeshVEGF, Difusion, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = [2.]
BcVEGF2 = CreateBoundary(MeshVEGF, Difusion, Phi2, BCRight, BCDirichlet)

Phi3 = [10.]
BcEC1 = CreateBoundary(MeshEC, CellVegf2, Phi3, BCLeft, BCDirichlet) # Dirichlet

Phi4 = [0.5]
BcEC2 = CreateBoundary(MeshEC, CellVegf2, collect(2:51:1532), Phi4, BCForces) # Manantial

BcVEGF = [BcVEGF1,BcVEGF2]
BcEC = [BcEC1, BcEC2]

modelVEGF = FemModel("VEGF", "2D", MeshVEGF, BcVEGF)
modelEC = FemModel("Endothelial cells", "2D", MeshEC, BcEC)

#d0_VEGF=SolveProgram(modelVEGF)

dir= "Desktop/Results"
#dir=homedir()* BIOSSFemCore.dirSep*dir
SolveProgramAngiogenesis(MeshVEGF,modelVEGF,MeshEC, modelEC, 1., 50.,dir; replace=true, append=true)

#SaveToParaview(MeshVEGF,[d0_VEGF], ["VEGF"], "VEGF0", dir; replace=true,append=true)
