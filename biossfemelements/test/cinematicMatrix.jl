using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements

# init_getBmat(_3D,Brick8,MultiDifusion{1})
# getBmat(_3D,Brick8,MultiDifusion{1},ones(3,3),zeros(3))

init_getNmat(_2D,Quad4,Difusion)
N = getNmat(_2D,Quad4,Difusion,[-1,-1.])

init_getNmat(_2D,Quad4,MultiDifusion{3})
N = getNmat(_2D,Quad4,MultiDifusion{3},[-1.,-1.])
display(N)


ClearStoredCache()
init_getBmat(_2D,Quad4,MultiDifusion{2})
B = getBmat(_2D,Quad4,MultiDifusion{2},[1 0; 0. 2],zeros(2))
display(Matrix(B))


ClearStoredCache()
init_getBmat(_2D,Quad4,CellVegf)
B = getBmat(_2D,Quad4,CellVegf,[1 0; 0. 2],zeros(2))
display(Matrix(B))

init_getBmat(_2D,Quad4,CellVegf2)
B = getBmat(_2D, Quad4,CellVegf2, [1 0; 0. 2], zeros(2))
display(Matrix(B))

init_getBmat(_2D,Quad4,Difusion)
B2 = getBmat(_2D,Quad4,Difusion,[1 0; 0. 2],zeros(2))
display(Matrix(B2))


init_getBmat(_3D,Brick8,Difusion)
B2 = getBmat(_3D,Brick8,Difusion,[1 0 0; 0. 1 0; 0 0 1],zeros(3))
display(Matrix(B2))


ClearStoredCache()
init_getBmat(_1D,Lin2,CellVegf)
B = getBmat(_1D,Lin2,CellVegf,[1],1)
println(B*[0., 0., 2., 0.])
display(Matrix(B))

init_getNmat(_2D,Quad4,CellVegf)
N = getNmat(_2D,Quad4,CellVegf, [1, 1])
display(N)

init_getNmat(_2D,Quad4,CellVegf2)
N = getNmat(_2D,Quad4,CellVegf2, [1, 1])
display(N)

init_getBmat(_1D,Lin2,Difusion)
B = getBmat(_1D,Lin2,Difusion, [1] , 0.)
println(B*[0., 1.])
