using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

init_getNmat(_2D,Quad4, ConvecDif) # inicio de la funcion de forma
init_getBmat(_2D,Quad4, ConvecDif) # inicio de la funcion cinematica
initJacobian(_2D,Quad4)

N = 400 # elementos en x
M = 400 # elementos en y

Mesh = RegularMesh(_2D, Quad4, ConvecDif, [N, M], [1., 1.]) # dimesiones de la malla

for el in Mesh.Elements
    el.Mat.Vel = [1,1]
    el.Mat.DMat .*= 1
end

Phi1 = [1.]
Bc1 = CreateBoundary(Mesh, ConvecDif, Phi1, BCRight, BCDirichlet) # Dirichlet

Phi2 = [0.]
Bc2 = CreateBoundary(Mesh, ConvecDif, Phi2, BCLeft, BCDirichlet)

Bc = [Bc1,Bc2]

p = FemModel("Convection", "2D", Mesh, Bc)

@time u = SolveProgram(p)

Output(Mesh,[u],["cosas"],"datos",homedir()*"/Desktop")

# using Plots
#
# pyplot()
# surface(reshape(u,(N+1),(M+1)),camera=(0,30))
