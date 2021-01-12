using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver


using Test

@testset "Test de tracción lineal" begin

    println("***")
    println("Iniciando ensayo a tracción")
    println("***")

    init_getBmat(_3D, Brick8, Elastic3D) # inicio de la funcion cinematica
    initJacobian(_3D, Brick8) # inicio del jacobiano

    N = 1 # elementos en x
    M = 1 # elementos en y
    O = 1 # elementos en z

    println("***")
    println("Creando malla")
    print("Tiempo total:")

    # Mesh creation
    @time Mesh = RegularMesh(_3D, Brick8, Elastic3D, [N, M, O], [1., 1., 1.]) # dimesiones de la malla
    println()

    A = 1 # area
    E = 1e9 # elasticity module

    for el in Mesh.Elements
        el.Mat.DMat .*= E
    end

    println("Creando condiciones de contorno")
    println()

    # Simulation of Traction Test

    Bc1 = CreateBoundary(Mesh, Elastic3D, [1], [0., 0., 0.], BCDirichlet) # Dirichlet
    Bc2 = CreateBoundary(Mesh, Elastic3D, [2], [missing, 0., 0.], BCDirichlet) # Dirichlet
    Bc3 = CreateBoundary(Mesh, Elastic3D, [4], [missing, missing, 0.], BCDirichlet) # Dirichlet
    Bc4 = CreateBoundary(Mesh, Elastic3D, [3], [0., missing, 0.], BCDirichlet) # Dirichlet
    Bc5 = CreateBoundary(Mesh, Elastic3D, [5], [0., 0., missing], BCDirichlet) # Dirichlet
    Bc6 = CreateBoundary(Mesh, Elastic3D, [6], [missing, 0., missing], BCDirichlet) # Dirichlet
    Bc8 = CreateBoundary(Mesh, Elastic3D, [7], [0.,missing, missing], BCDirichlet) # Dirichlet

    # vector of the different forces

    forces = collect(range(0,10e5,length = 21))

    function exportu(u,N)
        vars = []
        nombres = []

        for i in 1:N
            push!(vars, u[i:N:end])
            push!(nombres, "cosa $i")
        end
        return vars, nombres
    end

    for f in forces

        F1 = [0., 0., f]
        BForc = CreateBoundary(Mesh, Elastic3D, reshape(Mesh.Nodes[:,:,end] ,length(Mesh.Nodes[:,:,end])), F1, BCForces) # Fuerzas

        Bc = [Bc1, Bc2, Bc3, Bc4, Bc5, Bc6, Bc8, BForc]
        print("Fuerza total aplicada: ")

        FAplicada = 0.
        for ff in BForc.Cond[1]
            if ff !== missing
                FAplicada += ff
            end
        end

        println(FAplicada, "N")

        println("Creando modelo")
        p = FemModel("Elastic", "3D", Mesh, Bc)
        println("Resolviendo modelo")
        print("Tiempo total:")
        @time u = SolveProgram(p)

        vars, nombres = exportu(u,3)

        dis = sum(vars[3][5:end])/length(vars[3][5:end])

        print("Desplazamientos medios: ")
        println("ΔL = $dis")

        print("Desplazamientos estimados: ")
        println("ΔL = $(FAplicada/E/A)")

        println(dis)

        @test (dis-FAplicada/E/A) < 1e-12
        println()

    end

end
