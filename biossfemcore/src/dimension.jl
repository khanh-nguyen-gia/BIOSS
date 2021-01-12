"""
	This file contains the dimension-related
	DataTypes and the asociated functions.
"""

# Struct that indicates the spatial dimension
abstract type AbstractDim end

struct _1D <: AbstractDim
end

struct _2D <: AbstractDim
end

struct _3D <: AbstractDim
end


"""
	function getSymDim(::Type{D})

	Returns the dimension number (int) asociated
	to each type
"""
function getSymDim(::Type{D}) where {D <: AbstractDim}
	println("Todavía no implementado")
	return nothing
end

# 1D
function getSymDim(::Type{_1D})
	return 1
end

# 2D
function getSymDim(::Type{_2D})
	return 2
end

# 3D
function getSymDim(::Type{_3D})
	return 3
end


DimMap = Dict([(1, _1D), (2, _2D), (3,_3D)])
"""
	function getDimType(N::Int)

	Return the dimension type asociated to
	the integer N
"""
function getDimType(N::Int)
	global DimMap
	return DimMap[N]
end

"""
	function getSymVars(::Type{D})

	Return an expresion (::Expr) containing
	the variables for each dimension
"""
function getSymVars(::Type{D}) where {D <: AbstractDim}
	println("Todavía no implementado")
	return nothing
end


# variables for 1D
function getSymVars(::Type{_1D})
	return :(ξ)
end

# variables for 2D
function getSymVars(::Type{_2D})
	return :(ξ,η)
end

# variables for 3D
function getSymVars(::Type{_3D})
	return :(ξ,η,ν)
end


"""
	function getSymSymbol(::Type{D})

	Returns an array containing
	the variables symbol for each dimension.

	It is used to indicate the variables that
	will be differentiated
"""
function getSymSymbol(::Type{D}) where {D <: AbstractDim}
	println("Todavía no implementado")
	return nothing
end

# array for 1D
function getSymSymbol(::Type{_1D})
	return [:ξ]
end

# array for 2D
function getSymSymbol(::Type{_2D})
	return [:ξ,:η]
end

# array for 3D
function getSymSymbol(::Type{_3D})
	return [:ξ,:η,:ν]
end
