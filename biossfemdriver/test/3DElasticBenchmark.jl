using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

# ClearStoredCache()
#
# init_getBmat(_3D, Brick8, Elastic3D) # inicio de la funcion cinematica
# initJacobian(_3D, Brick8)

println("Iniciando test")
println("***")

using DelimitedFiles

NofTest = 6

file = open(homedir()*"/Desktop/Benchresults.txt","w")

@time for kk in collect(1:NofTest)

    global NofTest, file

    N = 2
    M = 2
    O = 2

    Mesh = RegularMesh(_3D, Brick8, Elastic3D, [N, M, O]*kk, [2., 0.5, 4.]) # dimesiones de la malla

    U1 = [0., 0., 0.]
    Bc1 = CreateBoundary(Mesh, Elastic3D, U1, vec(Mesh.Nodes[:,:,1]), BCDirichlet) # Dirichlet

    F1 = [0., 0.02, 0.]
    nodes = vec(Mesh.Nodes[end,:,Int(O/2):end])
    Bc2 = CreateBoundary(Mesh, Elastic3D, F1./length(nodes), nodes, BCForces) # Fuerzas

    nodes = vec(Mesh.Nodes[1,:,Int(O/2):end])
    Bc3 = CreateBoundary(Mesh, Elastic3D, -F1./length(nodes), nodes, BCForces) # Fuerzas


    Bc = [Bc1, Bc2, Bc3]

    p = FemModel("Elastic", "3D", Mesh, Bc)

    time = @elapsed u = SolveProgram(p)

    println("N = $(N*kk), time = $time")
    writedlm(file, [N*kk time])
end

close(file)

function exportu(u,N)
    vars = []
    nombres = []

    for i in 1:N
        push!(vars, u[i:N:end])
        push!(nombres, "cosa $i")
    end
    return vars, nombres
end

vars, nombres = exportu(u,3)

Output(Mesh,vars,nombres,"elastic",homedir()*"/Desktop")
OutputDisp3D(Mesh,vars,nombres,"elastic3D",homedir()*"/Desktop")
