using Revise

using BIOSSFemCore
using BIOSSFemMaterial

using BIOSSFemElements
using BIOSSFemDriver

NSpecies = 4
init_getNmat(_2D, Quad4, MultiDifusion{NSpecies}) # inicio de la funcion de forma
init_getBmat(_2D, Quad4, MultiDifusion{NSpecies}) # inicio de la funcion cinematica
initJacobian(_2D, Quad4)


N = 10 # elementos en x
M = 10 # elementos en y

println("***")
@time Mesh = RegularMesh(_2D, Quad4, MultiDifusion{NSpecies}, [N, M], [1., 1.]) # dimesiones de la malla
println("***")

Phi1 = ones(NSpecies)
Bc1 = CreateBoundary(Mesh, MultiDifusion{NSpecies}, Phi1, BCLeft, BCDirichlet) # Dirichlet

Phi2 = ones(NSpecies)*0.2
Bc2 = CreateBoundary(Mesh, MultiDifusion{NSpecies}, Phi2, BCRight, BCDirichlet)

Bc = [Bc1,Bc2]

p = FemModel("Difusion", "2D", Mesh, Bc)
@time u = SolveProgram(p)

function exportu(u,N)
    vars = []
    nombres = []

    for i in 1:N
        push!(vars, u[i:N:end])
        push!(nombres, "cosa $i")
    end
    return vars, nombres
end

vars, nombres = exportu(u, NSpecies)

using LinearAlgebra
for i in 2:length(vars)
    println(norm(vars[i-1] - vars[i])< 1e-14)
end

SaveToParaview(Mesh,vars,nombres,"Test2","/Desktop/")
