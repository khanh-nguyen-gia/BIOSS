using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

init_getBmat(_2D,Quad4,Difusion)
initJacobian(_2D,Quad4)
clearconsole()

@time Content = ImportNastranMesh(2, homedir()*"/Desktop/meshtest/mesh2.bdf")

@time FemMesh = NastranMesh(_2D, Quad4, Difusion, Content)

Nodes = FemMesh.Elements[1].Nodes

Nodes = vcat(FemMesh.Elements[1].Nodes, FemMesh.Elements[end].Nodes)
Bc1 = CreateBoundary(FemMesh, Difusion, Nodes, [0.], BCDirichlet)

Nodes = vcat(FemMesh.Elements[round(Int,end/2)+1].Nodes, FemMesh.Elements[round(Int,end/2)].Nodes)
Bc2 = CreateBoundary(FemMesh, Difusion, Nodes, [10.], BCForces)

for i in 1:length(Bc2.Cond)
    Bc2.Cond[i] .+= 1.
end

Bc = [Bc1,Bc2]

p = FemModel("Difusion", "2D", FemMesh, Bc)
@time u = SolveProgram(p)

using Plots

x,y = Content.GridCoords[1,:],Content.GridCoords[2,:]
p1 = Plots.surface(x,y,u,cam=(30,35),lw=0,c=cgrad([:orange,:yellow,:white]))

p2 = Plots.surface(x,y,u,cam=(0,90),lw=0,c=cgrad([:orange,:yellow,:white]))


display(p1)
display(p2)

println(maximum(u))
#
#p2 = Plots.contourf(x,y,u,c=cgrad([:orange,:yellow,:white]))
