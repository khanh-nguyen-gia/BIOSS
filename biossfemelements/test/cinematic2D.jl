using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements

ClearStoredCache()
init_getBmat(_2D,Quad4,Elastic2D)
B = getBmat(_2D,Quad4,Elastic2D,[1. 0. ; 0. 1.],[.4, .5])
display(B)


ClearStoredCache()
init_getBmat(_3D,Brick8,Elastic3D)
B = getBmat(_3D,Brick8,Elastic3D,[1. 0. 0 ; 0. 1. 0; 0. 0. 1.],[0., .0, .0])
display(Array(B))
