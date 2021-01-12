using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver

using Plots
using SparseArrays
init_getNmat(_2D,Quad4, Difusion)
init_getBmat(_2D,Quad4, Difusion)
initJacobian(_2D,Quad4)
clearconsole()

@time Content = ImportNastranMesh(2, homedir()*"/Desktop/meshtest/mesh2.bdf")

for i in 1:length(Content.ElementNodes)
    Content.ElementNodes[i] .-= minimum(Content.GridNodes)-1
end
Content.GridNodes .-= minimum(Content.GridNodes)-1

Content.ID .-= minimum(Content.ID) -1

# display(Content.GridNodes)
# scatter(Content.GridCoords[1,:].*10, Content.GridCoords[2,:].*10)

@time FemMesh = CreateMesh(_2D, Quad4, Difusion, Content)

Nodes = FemMesh.Elements[1].Nodes


Nodes = rand(1:length(Content.GridCoords),20)
Bc1 = CreateBoundary(FemMesh, Difusion, [0.], Nodes, BCDirichlet)

Nodes = rand(1:length(Content.GridCoords),20)
Bc2 = CreateBoundary(FemMesh, Difusion, [10.], Nodes, BCForces)

for i in 1:length(Bc2.Cond)
    Bc2.Cond[i] .+= 1.
end

Bc = [Bc1,Bc2]

model = FemModel("Difusion", "2D", FemMesh, Bc)

Assembly!(model)
AssemblyMass!(model)

d0 = uval(model)
f = Forces(model)

d0[model.FreeNodes] += rand(length(model.FreeNodes))*3

v0 = sparse(model.M) \ (f - sparse(model.K)*d0)

# println(v0)

dt = 0.0008

if ispath(homedir()*"/Desktop/results/")
    rm(homedir()*"/Desktop/results/", recursive=true)
end
mkdir(homedir()*"/Desktop/results/")

for i in 1:500
    global d0, v0, v1, d1

    d1, v1 = dFormIntegration(model, dt, f, d0, v0, alpha = 1.)
    v0 = v1; d0 = d1

    if mod(i,5) == 0
        x,y = Content.GridCoords[1,:],Content.GridCoords[2,:]
        p1 = Plots.surface(x,y,d0,cam=(30,90),lw=0,c=cgrad([:orange,:yellow,:white]))

        # p2 = Plots.surface(x,y,d0,cam=(0,90),lw=0,c=cgrad([:orange,:yellow,:white]))
        display(p1)
        # display(p2)
        println(maximum(d0))
    end
end





#
#p2 = Plots.contourf(x,y,u,c=cgrad([:orange,:yellow,:white]))
