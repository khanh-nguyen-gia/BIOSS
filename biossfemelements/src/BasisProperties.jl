"""
	function get_basis_gdl(::Type{B})

	For a given basis, returns the number of
	nodes

	Input:

		Type{B}  -> Basis Type

	Output:

		Npoints :: Int -> numer of nodes of each element
"""
function get_basis_gdl(::Type{B}) where {B<: AbstractBasis}
	println("todavia no implementado")
end


"""

	function get_basis_props(::Type{B})

	For a given basis, returns the coords of
	the nodes and the basic form of the
	shape functions

	Input:

		Type{B}  -> Basis Type

	Output:

		Coords :: Array{Float64,2}   -> Natural coordinates of the nodes
		Code :: String             -> Shape functions definition
"""
function get_basis_props(::Type{B}) where {B<: AbstractBasis}
	println("todavia no implementado")
end

"""
	function get_int_points(Type{B})

	For a given element, returns the coords of
	the integration points and the corresponding
	weights

	Input:

		Type{B}  -> Basis Type

	Output:

		IntPoints :: Array{Float64,2}   -> Natural coordinates
										 of the integration points

		Weights :: Vector{Float64}      -> Weights of each int. point
"""
function get_int_points(::Type{B},NPoints) where {B<: AbstractBasis}
	println("todavia no implementado")
	return nothing
end
