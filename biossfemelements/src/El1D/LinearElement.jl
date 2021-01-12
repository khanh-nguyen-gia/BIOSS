"""

    This file  contains 1D linear elements

"""


struct Lin2 <: AbstractBasis
end


struct Lin3 <: AbstractBasis
end


function get_basis_gdl(::Type{Lin2})
	return 2
end

function get_basis_gdl(::Type{Lin3})
	return 3
end


function get_basis_props(::Type{Lin2})
	# coordenadas nodales
    #        | x  | y |
    #        |--------|
    coords = [ -1. # node 1
                1.] # node 2


	# Quad4 shape function:
	code = Meta.parse("1 + ξ")

    return coords, code
end

function get_basis_props(::Type{Lin3})
	# coordenadas nodales
    #        | x  | y |
    #        |--------|
    coords = [ -1. # node 1
				0. # node 2
                1.] # node 3


	# Quad4 shape function:
	code = Meta.parse("1 + ξ + ξ*ξ")

    return coords, code
end


function get_int_points(::Type{Lin2},NPoints)

	# coordenadas de los puntos
	# de integracion
    #             | x |
    #             |---|
	IntPoints = [  -1.; # 1
                    1.] # 2
	IntPoints *= 1. /√3

	# pesos de cada punto: w

	Weights = [1,1]

	return IntPoints,Weights
end

function get_int_points(::Type{Lin3},NPoints)

	# coordenadas de los puntos
	# de integracion
    #             | x |
    #             |---|
	IntPoints = [  -sqrt(3/5); # 1
					0.; 	   # 2
                    sqrt(3/5)] # 3

	# pesos de cada punto: w

	w=[5/9., 8/9., 5/9.];

	return IntPoints,Weights
end
