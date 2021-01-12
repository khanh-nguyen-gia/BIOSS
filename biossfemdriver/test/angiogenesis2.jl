using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

using SparseArrays

dir= "Desktop/Results2"

init_getNmat(_2D,Quad4, CellVegf2) # inicio de la funcion de forma
init_getNmat(_2D,Quad4, VEGF)
init_getNmat(_2D, Quad4, Difusion)

init_getBmat(_2D,Quad4, CellVegf2) # inicio de la funcion cinematica
init_getBmat(_2D,Quad4, VEGF)
init_getBmat(_2D,Quad4, Difusion)

initJacobian(_2D,Quad4)

N = 50 # elementos en x
M = 50 # elementos en y

x̂ = 1. # x_max en mms adimensionalizada con x0 = 2mm
ŷ = 1. # y_max en mms adimensionalizada con y0 = 2mm

tmax= 50. # every time unit corresponds to 1.5 days

MeshVEGF_dif = RegularMesh(_2D, Quad4, Difusion, [N, M], [x̂, ŷ]) # Malla comportamiento VEGF

Phi1 = [0.5]
BcVEGF_dif_1 = CreateBoundary(MeshVEGF_dif, Difusion, Phi1, BCLeft, BCDirichlet) # Dirichlet
Phi2 = [15.]
BcVEGF_dif_2 = CreateBoundary(MeshVEGF_dif, Difusion, Phi2, BCRight, BCDirichlet)
BcVEGF_dif =[BcVEGF_dif_2]

modelVEGF_dif = FemModel("VEGF_dif", "2D", MeshVEGF_dif, BcVEGF_dif)
d0_VEGF_dif = SolveProgram(modelVEGF_dif)

SaveToParaview(MeshVEGF_dif,[d0_VEGF_dif], ["VEGF_diffusion"], "VEGF_dif",
dir; replace=true,append=true)




MeshVEGF = RegularMesh(_2D, Quad4, VEGF, [N, M], [x̂, ŷ]) # Malla comportamiento VEGF

Phi1 = [0.5]
BcVEGF1 = CreateBoundary(MeshVEGF, VEGF, Phi1, BCLeft, BCDirichlet) # Dirichlet
Phi2 = [15.]
BcVEGF2 = CreateBoundary(MeshVEGF, VEGF, Phi2, BCRight, BCDirichlet)
BcVEGF =[BcVEGF1, BcVEGF2]

modelVEGF = FemModel("VEGF", "2D", MeshVEGF, BcVEGF)
d0_VEGF = SolveProgram(modelVEGF)
SaveToParaview(MeshVEGF,[d0_VEGF], ["VEGF_start"], "VEGF_start",
dir; replace=true,append=true)




MeshEC = RegularMesh(_2D, Quad4, CellVegf2, [N, M], [x̂, ŷ]) # Malla comportamiento celulas endoteliales

Phi3 = [15.]
BcEC1 = CreateBoundary(MeshEC, CellVegf2, Phi3, BCLeft, BCDirichlet) # Dirichlet
Phi4 = [1.]
BcEC2 = CreateBoundary(MeshEC, CellVegf2, collect(2:51:1532), Phi4, BCForces) # Manantial
BcEC = [BcEC2]

modelEC = FemModel("Endothelial cells", "2D", MeshEC, BcEC)

#dir=homedir()* BIOSSFemCore.dirSep*dir
SolveProgramAngiogenesis2(MeshVEGF,modelVEGF,MeshEC, modelEC, 1., tmax, dir; replace=true, append=true)
