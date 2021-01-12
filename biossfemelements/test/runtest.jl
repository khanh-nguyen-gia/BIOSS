using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements

init_getNmat(_2D,Quad4,Heat)

nodes = collect(1:8)

el3d = Element(_3D,Brick8,Heat,1,nodes)

f = getShape(el3d)
Df = getShapeGrad(el3d)

nodes = collect(1:4)

el2d = Element(_2D,Quad4,Heat,1,nodes)

el2d.Coords = [0. 0.; 1. 0.; 1. 1.; 0. 1.]

el3d.Coords = [0. 0. 0.; 1. 0. 0.; 1. 1. 0.; 0. 1. 0.;
               0. 0. 1.; 1. 0. 1.; 1. 1. 1.; 0. 1. 1.]

coords, w = get_basis_props(el3d)
initJacobian(el3d)
initJacobian(el2d)
println(ElementJacobian(el3d,[-1.,-1.,-1]))

coords, w = get_basis_props(el2d)

# @time for i in 1:10000
#     ElementJacobian(el2d,[-1.,-1.])
# end
#
# @time for i in 1:10000
#     ElementJacobian(el3d,[-1.,-1.,-1.])
# end

Mmat = LocalMass(el2d)

display(Mmat)
