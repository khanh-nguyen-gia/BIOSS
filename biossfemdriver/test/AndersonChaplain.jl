using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

using SparseArrays

dir= "Desktop/AndersonChaplain"

init_getNmat(_2D,Quad4, CellVegf2) # inicio de la funcion de forma
init_getNmat(_2D,Quad4, VEGF)
init_getNmat(_2D, Quad4, Fibronectin)

init_getBmat(_2D,Quad4, CellVegf2) # inicio de la funcion cinematica
init_getBmat(_2D,Quad4, VEGF)
init_getBmat(_2D,Quad4, Fibronectin)

initJacobian(_2D,Quad4)

N = 50 # elementos en x
M = 50 # elementos en y

x̂ = 1. # x_max en mms adimensionalizada con x0 = 2mm
ŷ = 1. # y_max en mms adimensionalizada con y0 = 2mm

tmax= 15. # every time unit corresponds to 1.5 days

# Since the equations have been adimentionalized, there is only two values we study for BCs
Phi1 = [1.]
Phi2 = [0.]

#Creation of Fibronectin model
MeshFiber = RegularMesh(_2D, Quad4, Fibronectin, [N, M], [x̂, ŷ]) # Malla

Bc_Fiber_1 = CreateBoundary(MeshFiber, Fibronectin, [0.8], BCLeft, BCDirichlet) # Dirichlet
Bc_Fiber_2 = CreateBoundary(MeshFiber, Fibronectin, Phi2, BCRight, BCDirichlet)
Bc_Fiber =[Bc_Fiber_1, Bc_Fiber_2]

modelFiber = FemModel("Fiber", "2D", MeshFiber, Bc_Fiber)


# Creation of VEGF model
MeshVEGF = RegularMesh(_2D, Quad4, VEGF, [N, M], [x̂, ŷ]) # Malla comportamiento VEGF

BcVEGF1 = CreateBoundary(MeshVEGF, VEGF, Phi1, BCRight, BCDirichlet) # Dirichlet
BcVEGF2 = CreateBoundary(MeshVEGF, VEGF, Phi2, BCLeft, BCDirichlet)
BcVEGF =[BcVEGF1, BcVEGF2]

modelVEGF = FemModel("VEGF", "2D", MeshVEGF, BcVEGF)



# Creation of EC model

MeshEC = RegularMesh(_2D, Quad4, CellVegf2, [N, M], [x̂, ŷ]) # Malla comportamiento celulas endoteliales

BcEC1 = CreateBoundary(MeshEC, CellVegf2, Phi2, BCLeft, BCForces) # Flujo nulo
BcEC = [BcEC1]

modelEC = FemModel("Endothelial cells", "2D", MeshEC, BcEC)


SolverAngiogenesisAC(modelVEGF, modelFiber, MeshEC, modelEC, 0.01, tmax, dir; replace=true, append=true)
