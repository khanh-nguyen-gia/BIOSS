struct Brick8 <: AbstractBasis
end


function get_basis_gdl(::Type{Brick8})
	return 8
end

function get_basis_props(::Type{Brick8})
	# coordenadas nodales
    #        | x | y |  z |
    #        |------------|
    coords = [ -1. -1. -1.; # node 1
                1. -1. -1.; # node 2
                1.  1. -1.; # node 3
               -1.  1. -1.; # node 4
               -1. -1.  1.; # node 5
                1. -1.  1.; # node 6
                1.  1.  1.; # node 7
               -1.  1.  1.] # node 8

	# Brick8 shape function:
	code = Meta.parse("1 + ξ + η + ν + ξ*η + η*ν + ν*ξ + ξ*η*ν")

    return coords, code
end

# return the coordinates and weigths of
# the integration points
function get_int_points(::Type{Brick8}, NPoints)

	# coordenadas de los puntos
	# de integracion
    #             | x  | y |
    #             |--------|
	IntPoints = [ -1. -1. -1.; #  1
                   1. -1. -1.; #  2
                   1.  1. -1.; #  3
                  -1.  1. -1.; #  4
                  -1. -1.  1.; #  5
                   1. -1.  1.; #  6
                   1.  1.  1.; #  7
                  -1.  1.  1.] #  8

	IntPoints *= 1. /√3

	Weights = [1,1,1,1,1,1,1,1]

	return IntPoints,Weights
end
