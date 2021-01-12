using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver


using LinearAlgebra
using Test


"""
    How to use Test:

        1) Create function that you want to test
        2) Create comparation data
        3) Create the macro @testset: it will include all the tests
        4) Run each test wich the macro @test. You have to pass
        a boolean result to test (true or false).

"""

# linear distribution heat equation
function LinearHeatTest(N, M, xfin, Phi1, Phi2)

    yfin = 1. # posición y de la esquina

    FemMesh = RegularMesh(_2D, Quad4, Difusion, [N, M], [xfin, yfin]) # dimesiones de la malla

    Bc1 = CreateBoundary(FemMesh, Difusion, Phi1, BCLeft, BCDirichlet) # Dirichlet

    Bc2 = CreateBoundary(FemMesh, Difusion, Phi2, BCRight, BCDirichlet) # Dirichlet

    Bc = [Bc1,Bc2]

    p = FemModel("Difusion", "2D", FemMesh, Bc)

    @time u = SolveProgram(p)

    SaveToParaview(FemMesh,[u],["Calor"],"Test","/Desktop/")

    return u
end
# pseudo analytic solution for linear distribution in heat equation
function AnalyticSol(N, M, xfin, Phi1, Phi2)

    x = collect(range(0, xfin, length=N+1))

    ResComp = Phi1[1] .+ (Phi2[1]-Phi1[1]).*x/xfin
    TestRes = []

    for j in 1:M+1
        TestRes = vcat(TestRes, ResComp)
    end

    return TestRes
end

@testset "Test de la ecuación del calor" begin

    init_getNmat(_2D,Quad4, Difusion) # inicio de la funcion de forma
    init_getBmat(_2D,Quad4, Difusion) # inicio de la funcion cinematica
    initJacobian(_2D,Quad4)

    # SaveToParaview(FemMesh,[u],["cosas"],"datos2","Desktop/Test"; replace=false, append=false)
    N = 10; M = 10; xfin = 2.
    Phi1 = [1.]; Phi2 = [2.]

    u = LinearHeatTest(N, M, xfin, Phi1, Phi2)
    Testu = AnalyticSol(N, M, xfin, Phi1, Phi2)

    @test maximum(u-Testu) < 1e-12

    N = 100; M = 10; xfin = 10.
    Phi1 = [1.]; Phi2 = [-2.]

    u = LinearHeatTest(N, M, xfin, Phi1, Phi2)
    Testu = AnalyticSol(N, M, xfin, Phi1, Phi2)

    @test maximum(u-Testu) < 1e-12

    N = 100; M = 100; xfin = 10.
    Phi1 = [100.]; Phi2 = [-2.]

    u = LinearHeatTest(N, M, xfin, Phi1, Phi2)
    Testu = AnalyticSol(N, M, xfin, Phi1, Phi2)

    @test maximum(u-Testu) < 1e-10
end
