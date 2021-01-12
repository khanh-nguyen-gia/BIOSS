struct Quad4 <: AbstractBasis
end

struct Quad8 <: AbstractBasis
end

struct Quad9 <: AbstractBasis
end


function get_basis_gdl(::Type{Quad4})
	return 4
end

function get_basis_gdl(::Type{Quad8})
	return 8
end

function get_basis_gdl(::Type{Quad9})
	return 9
end



function get_basis_props(::Type{Quad4})
	# coordenadas nodales
    #        | x  | y |
    #        |--------|
    coords = [ -1. -1.; # node 1
                1. -1.; # node 2
                1.  1.; # node 3
               -1.  1.] # node 4

	# Quad4 shape function:
	code = Meta.parse("1 + ξ + η + ξ*η")

    return coords, code
end

function get_basis_props(::Type{Quad8})
	# coordenadas nodales
    #        | x  | y |
    #        |--------|
    coords = [ -1. -1.; # node 1
                1. -1.; # node 2
                1.  1.; # node 3
               -1.  1.; # node 4
                0. -1.; # node 5
                1.  0.; # node 6
                0.  1.; # node 7
               -1.  0.] # node 8

	# Quad8 shape function:

	code = Meta.parse("1 + ξ + η + ξ*η + ξ*ξ + η*η + ξ*ξ*η + η*η*ξ ")

    return coords, code
end

function get_basis_props(::Type{Quad9})
	# coordenadas nodales
    #        | x  | y |
    #        |--------|
	coords = [ -1. -1.; # node 1
                1. -1.; # node 2
                1.  1.; # node 3
               -1.  1.; # node 4
                0. -1.; # node 5
                1.  0.; # node 6
                0.  1.; # node 7
               -1.  0.; # node 8
               -0.  0.] # node 9

	# Quad9 shape function:
	code = Meta.parse("1 + ξ + η + ξ*η + ξ*ξ + η*η + ξ*ξ*η + η*η*ξ + ξ*ξ*η*η ")

	n_points = 9

    return coords, code
end


function get_int_points(::Union{Type{Quad8},Type{Quad9}},Npoints)


	w=[5/9, 8/9, 5/9];
	pG=[-sqrt(3/5), 0, sqrt(3/5)]

	IntPoints = zeros(9,2)
	Weights = zeros(9)

	k = 0

	for i in 1:length(w)
		for j in 1:length(w)
			k += 1
			IntPoints[k,:] = [pG[i] pG[j]]
			Weights[k] = w[i]*w[j]
		end
	end

	return IntPoints,Weights
end


# TODO: Mejorar esta estructura, para hacerla más eficiente
# (Hacer benchmark con las opciones)

function get_int_points(::Type{Quad4},NPoints)

	# coordenadas de los puntos
	# de integracion
    #             | x  | y |
    #             |--------|
	IntPoints = [  -1  -1.; #  1
                    1. -1.; #  2
                    1.  1.; #  3
                   -1.  1.] #  4

	IntPoints *= 1. /√3

	# pesos de cada punto: w_x[i]*w_y[i]

	Weights = [1,1,1,1]


	return IntPoints,Weights
end


# int points en el corntorno
function get_BD_int_points(::Type{Quad4},NPoints)

	# coordenadas de los puntos
	# de integracion
    #             | Γ  |
    #             |----|
	IntPoints = [  -1.,
				    1.]

	IntPoints *= 1. /√3

	Weights = [1,1]

	return IntPoints,Weights
end
