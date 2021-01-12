using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements

El1 = Element(_1D, Lin2, Difusion, 1, [1,2])
El1.Coords = [0. 1.]

clearconsole()
println(typeof(El1.Coords))
c, b = get_basis_props(El1)
display(c)

init_getNmat(_1D, Lin2, Difusion)

initJacobian(El1)

init_getBmat(_1D,Lin2,Difusion)
getBmat(El1, 1, 1)
getNmat(El1, -1.)
