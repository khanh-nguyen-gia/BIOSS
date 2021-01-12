using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

using SparseArrays
using LinearAlgebra

init_getNmat(_2D,Quad4, Heat) # inicio de la funcion de forma
init_getBmat(_2D,Quad4, Heat) # inicio de la funcion cinematica
initJacobian(_2D,Quad4)

N = 40 # elementos en x
M = 40 # elementos en y

Mesh = RegularMesh(_2D, Quad4, Heat, [N, M], [2., 1.]) # dimesiones de la malla
Phi1 = [5.]
Bc1 = CreateBoundary(Mesh, Heat, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = [10.]
Bc2 = CreateBoundary(Mesh ,Heat, Phi2, BCRight, BCDirichlet)

Bc3 = CreateBoundary(Mesh, Heat,  [9], [1.], BCForces)

Bc = [Bc1,Bc2,Bc3]

model = FemModel("Difusion", "2D", Mesh, Bc)

Assembly!(model)
AssemblyMass!(model)

empty!(model.K)

d0 = uval(model)
f = Forces(model)

d0[model.FreeNodes] += rand(length(model.FreeNodes))*3

v0 = sparse(model.M) \ (f - sparse(model.K)*d0)

# println(v0)

dt = 0.01

if ispath(homedir()*"/Desktop/results/")
    rm(homedir()*"/Desktop/results/", recursive=true)
end
mkdir(homedir()*"/Desktop/results/")

for i in 1:500
    global d0, v0, v1, d1

    d1, v1 = dFormIntegration(model, dt, f, d0, v0, alpha = 1.)
    v0 = v1; d0 = d1

    if mod(i,5) == 0
        SaveToParaview(Mesh,[d1], ["Temp"], "temp$i", "/Desktop/results")
    end
end
