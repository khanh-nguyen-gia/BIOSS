using BIOSSFem

using Test

@testset "Example" begin

    init_getBmat(_2D,Quad4, Difusion) # inicio de la funcion de forma

    N = 10 # elementos en x
    M = 10 # elementos en y

    Mesh = RegularMesh(_2D, Quad4, Difusion, [N, M], [1., 1.]) # dimesiones de la malla

    Phi1 = 0.7
    Bc1 = CreateBoundary(Mesh, Phi1, "Left") # Dirichlet

    Phi2 = 0.2
    Bc2 = CreateBoundary(Mesh, Phi2, "Right")

    Bc = [Bc1,Bc2]

    p = FemModel("Difusion", "2D", Mesh, Bc)

    u = SolveProgram(p)

    println(u)

    @test u[1] == 0.7

end
