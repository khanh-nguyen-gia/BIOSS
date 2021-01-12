using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver


## Difusion/Heat test

initJacobian(_1D, Lin2)
init_getNmat(_1D, Lin2, Difusion)
init_getBmat(_1D, Lin2, Difusion)

N = 10 # elementos en x

Mesh = RegularMesh(_1D, Lin2, Difusion, [N], [1.]) # dimesiones de la malla

# getBmat(Mesh.Elements[1], 1/sqrt(3), 1.)

Phi1 = [4.]
Bc1 = CreateBoundary(Mesh, Difusion, Phi1, [1], BCDirichlet) # Dirichlet

Phi2 = [10.]
Bc2 = CreateBoundary(Mesh, Difusion, Phi2, [N+1], BCDirichlet)
Bc = [Bc1,Bc2]

p = FemModel("Difusion", "1D", Mesh, Bc)

@time u = SolveProgram(p)
println(u)


## Elastic Test

N = 3 # elementos en x

Mesh = RegularMesh(_1D, Lin2, Elastic1D, [N], [1.]) # dimesiones de la malla

for el in Mesh.Elements
    el.Mat.DMat .= 1e9*1e-4
end

initJacobian(_1D, Lin2)
init_getBmat(_1D, Lin2, Elastic1D)

U1 = [0.]
Bc1 = CreateBoundary(Mesh, Elastic1D, U1, [1], BCDirichlet) # Dirichlet

F1 = [25e3]
Bc2 = CreateBoundary(Mesh, Elastic1D, F1, [N+1], BCForces)
Bc = [Bc1,Bc2]

p = FemModel("Difusion", "1D", Mesh, Bc)

@time u = SolveProgram(p)

# display(Matrix(sparse(p.K)))
println(u)
