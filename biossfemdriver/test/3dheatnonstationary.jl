using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

using SparseArrays
using LinearAlgebra

init_getNmat(_3D,Brick8, Heat) # inicio de la funcion de forma
init_getBmat(_3D,Brick8, Heat) # inicio de la funcion cinematica
initJacobian(_3D,Brick8)

N = 20 # elementos en x
M = 20 # elementos en y
O = 5

Mesh = RegularMesh(_3D, Brick8, Heat, [N, M, O], [2., 1., 2.]) # dimesiones de la malla
Phi1 = [0.]
Bc1 = CreateBoundary(Mesh, Heat, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = [10.]
Bc2 = CreateBoundary(Mesh ,Heat, Phi2, BCRight, BCDirichlet)

Bc = [Bc1,Bc2]

model = FemModel("Difusion", "2D", Mesh, Bc)

Assembly!(model)
AssemblyMass!(model)

d0 = uval(model)

d0[model.FreeNodes] += rand(length(model.FreeNodes))*10
f = Forces(model)

v0 = sparse(model.M) \ (f - sparse(model.K)*d0)

# println(v0)

dt = 0.01


if ispath(homedir()*"/Desktop/results/")
    rm(homedir()*"/Desktop/results/", recursive=true)
end
mkdir(homedir()*"/Desktop/results/")
for i in 1:500
    global d0, v0, v1, d1

    d1, v1 = vFormIntegration(model, dt, f, d0, v0, alpha = 0.)
    v0 = v1; d0 = d1

    if mod(i,5) == 0
        println("Guardando el paso $i")
        Output(Mesh,[d1], ["Temp"], "temp$i", homedir()*"/Desktop/results/")
    end
end
