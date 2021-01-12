using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver


NSpecies = 2

ClearStoredCache()
init_getNmat(_3D,Brick8, MultiDifusion{NSpecies}) # inicio de la funcion de forma

# init_getNmat(_2D,Quad4, MultiDifusion{1}) # inicio de la funcion de forma
init_getBmat(_3D,Brick8, MultiDifusion{NSpecies}) # inicio de la funcion cinematica
initJacobian(_3D,Brick8)

N = 20 # elementos en x
M = 20 # elementos en y
O = 10
println("***")
@time Mesh = RegularMesh(_3D,Brick8, MultiDifusion{NSpecies}, [N, M, O], [1., 1., .2]) # dimesiones de la malla
println("***")

Phi1 = ones(NSpecies)
Bc1 = CreateBoundary(Mesh, MultiDifusion{NSpecies}, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = ones(NSpecies)*0.2
Bc2 = CreateBoundary(Mesh, MultiDifusion{NSpecies}, Phi2, BCRight, BCDirichlet)

Bc = [Bc1,Bc2]

p = FemModel("Difusion", "2D", Mesh, Bc)
println("***")
@time u = SolveProgram(p)
println("***")

function exportu(u,N)
    vars = []
    nombres = []

    for i in 1:N
        push!(vars, u[i:N:end])
        push!(nombres, "Especie $i")
    end
    return vars, nombres
end

vars, nombres = exportu(u,NSpecies)

# using LinearAlgebra
# for i in 2:length(vars)
#     println(vars[i-1] == vars[i])
#     println(norm(vars[i-1] -vars[i]))
# end

SaveToParaview(Mesh,vars,nombres,"datos2","Desktop/results/")
