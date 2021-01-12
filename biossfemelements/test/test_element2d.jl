using Test

using material
using elements

@testset "Element 2d" begin

    el1 = Element(_2D,Quad4,Difusion,[1,2,3,4])

    el1.Coords = [ -1. -1.; # node 1
                1. -1.; # node 2
                1.  1.; # node 3
               -1.  1.] # node 4

    k = LocalStiffness(el1)

    println(k)

end
