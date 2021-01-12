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
O = 10 # elementos en z

Endpoint = [5., 5., 1.]
@time Mesh = RegularMesh(_3D, Brick8, Difusion, [N, M, O], Endpoint) # dimesiones de la malla

Mesh.Elements[1].Mat.DMat

for el in Mesh.Elements
    el.Mat.DMat .*= 1e-5
end

Phi1 = [7.]
Bc1 = CreateBoundary(Mesh, Difusion, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = [7.]
Bc2 = CreateBoundary(Mesh, Difusion, Phi2, BCRight, BCDirichlet)

Fval = zeros([N, M, O].+1...)

# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002533

BrainDens = 1.046 # g/cm3
BrainWeigth = 1400 # g
Glialcells = 8.4e10

BrainVolume = BrainWeigth/BrainDens
volume = prod(Endpoint) #cm3

EstCells = Glialcells/BrainVolume*volume

Ncells = 1e5
for i in 1:Ncells
    pos = [rand(1:N), rand(1:M), rand(1:O) ]
    Fval[CartesianIndex(pos...)] .-=  1e-12*EstCells/Ncells
end


# Fval = zeros([N, M, O].+1...) .-1

println(sum(Fval)/length(Fval))

Bc3 = CreateBoundary(Mesh, Fval)

Bc = [Bc1,Bc2,Bc3]

p = FemModel("Difusion", "3D", Mesh, Bc)

f = Forces(p)

println("Tiempo total:")
@time u = SolveProgram(p)

O2 = reshape(u,(N+1,M+1,O+1))



# y = zeros(N+1)

y = O2[:,round(Int, M+1/2), round(Int, O+1/2)]
# for j in 1:M+1, k in 1:O+1
#     global O2, y
#     y += O2[:,j,k]
# end

using Plots

plot(y)

contourf(O2[:,:,round(Int, O+1/2)])
