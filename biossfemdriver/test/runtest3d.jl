using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

init_getNmat(_3D,Brick8, Difusion) # inicio de la funcion de forma
init_getBmat(_3D,Brick8, Difusion) # inicio de la funcion cinematica
initJacobian(_3D,Brick8)

N = 40 # elementos en x
M = 40 # elementos en y
O = 40 # elementos en z

println("Creacion de malla")
@time Mesh = RegularMesh(_3D, Brick8, Difusion, [N, M, O], [1., 1., 0.1]) # dimesiones de la malla

Phi1 = [4.]
Bc1 = CreateBoundary(Mesh, Difusion, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = [10.]
Bc2 = CreateBoundary(Mesh, Difusion, Phi2, BCRight, BCDirichlet)

Bc = [Bc1,Bc2]

p = FemModel("Difusion", "3D", Mesh, Bc)

println("Tiempo total:")
@time u = SolveProgram(p)
println(minimum(u))

Output(Mesh,[u],["cosas"],"datos3D",homedir()*"/Desktop/results/")
