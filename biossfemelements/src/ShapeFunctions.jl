"""
	function VandermondeMatrix(f::Exprnode::Int, ::Type{B}, ::Type{D}, coords::Array{Float64,2})

	Get the Vandermonde matrix, used to calculate the shape functions

	Input:

		f :: Expr                 -> Symbolic expresion of the shape functions
		Type{D}					  -> Dimension (_1D, _2D, _3D)
		Type{B}					  -> Basis of the element
		coords::Array{Float64,2}  -> Coordinates of the element nodes

	Output:

		VMat :: Array{Float64,2}    -> Vandermonde matrix
"""
function VandermondeMatrix(f::Expr,
				 ::Type{D}, ::Type{B},
				 coords::Array{Float64,2}) where {D <: AbstractDim, B <: AbstractBasis}

	N = get_basis_gdl(B)

	# init the Vandermode Matrix
    VMat = zeros(N,N)

	# get the symbolic variables, depending of the dimension
	vars = getSymVars(D)

	# construct the matrix:
    for i in 1:N
        for j in 1:N
            eval(:( $(vars) = $(coords[i,:]) ))
            VMat[i,j] = eval(f.args[j+1])
        end
    end

    return VMat
end

function VandermondeMatrix(f::Expr, ::Type{D},
				 ::Type{B}, coords::Vector{Float64}) where {B <: AbstractBasis, D <: AbstractDim}

	N = get_basis_gdl(B)

 	# init the Vandermode Matrix
    VMat = zeros(N,N)

 	# get the symbolic variables, depending of the dimension
 	vars = getSymVars(D)

 	# construct the matrix:
    for i in 1:N
        for j in 1:N
            eval(:( $(vars) = $(coords[i]) ))
            VMat[i,j] = eval(f.args[j+1])
        end
    end

    return VMat
end

"""
	function NodeShapeFunction(f::Expr, Node_i, v::Array{Float64,2})

	Get the shape function at a given node (Node_i). The shape function must
	equal 1 in the Node_i and 0 in the rest of the nodes.

	Input:

		f :: Expr              -> Symbolic expresion of the shape functions
		Node_i :: Int          -> Node where the shape function is calculated
		VMat::Array{Float64,2} -> Vandermode Matrix

	Output:
		f_i :: Expr 			-> Symbolic expresion of the shape function
								   at Node_i
"""
function NodeShapeFunction(f::Expr, Node_i::Int, VMat::Array{Float64,2})

	# Init the system for the Node_i
    N = length(VMat[1,:])
    b = zeros(N);
    b[Node_i] = 1.

    fi = Expr(:call,:+)

	# Solve the system to get the shape function coeficients
    a = VMat\b

	# Assing the coeficients to each given component
    for j in 1:N
        push!(fi.args, Expr(:call,:*, f.args[j+1],a[j]))
    end

    return simplify(fi)
end

"""
	function getShape(::Type{D},::Type{B})

	Get the symbolic shape functions of an element, which basis is B
	and wich dimension is D.

	Input:

		Type{D}					  -> Dimension (_1D, _2D, _3D)
		Type{B}					  -> Basis of the element

	Output:
		ShapeFunctions :: Array   -> Array that contains the symbolic expresion
									 of the shape functions
"""

function getShape(::Type{D},::Type{B}) where {D <: AbstractDim, B <: AbstractBasis}

	# Init the array that stores all the shape functions
	ShapeFunctions = []

	# Get the basis degrees of freedom
	BasisNDof = get_basis_gdl(B)

	# Get the basis properties: Coordinates and symbolic shape functions
    (coords, f) = get_basis_props(B)

	# Calculate the Vandermode Matrix
	VMat = VandermondeMatrix(f, D, B, coords)

    for Node_i in 1:BasisNDof
        code = NodeShapeFunction(f, Node_i, VMat)
        push!(ShapeFunctions,code)
    end

    return ShapeFunctions
end

"""
	function getShapeGraded(::Type{D},::Type{B})

	Get the symbolic shape functions derivates (gradient) of an element,
	which basis is B and wich dimension is D.

	Input:

		Type{D}					  -> Dimension (_1D, _2D, _3D)
		Type{B}					  -> Basis of the element

	Output:

		ShapeGradient :: Array   -> Array that contains the symbolic expresion
									   of the shape functions gradient
"""
function getShapeGrad(::Type{D},::Type{B}) where {D <: AbstractDim, B <: AbstractBasis}
    ShapeFunctions = getShape(D,B)
    return getShapeGrad(ShapeFunctions,D)
end


"""
	function getShapeGraded(ShapeFunctions::Array, ::Type{D},::Type{B})

	Get the given the symboluc shape functions, calculate the symbolic
	derivates (gradient) of an element, which basis is B and wich dimension is D.

	Input:

		ShapeFunctions :: Array   -> Array that contains the symbolic expresion
								     of the shape functions
		Type{D}					  -> Dimension (_1D, _2D, _3D)
		Type{B}					  -> Basis of the element

	Output:

		ShapeGradient :: Array   -> Array that contains the symbolic expresion
								   of the shape functions gradient
"""
function getShapeGrad(ShapeFunctions::Array, ::Type{D}) where {D <: AbstractDim}

	# Init the array that stores all the gradients
	ShapeGradient = []

	vars = getSymSymbol(D)

    for ShapeFun_i in ShapeFunctions
		# calculate the gradient of the ShapeFunction i
        Gradient_i = differentiate(ShapeFun_i, vars)
		# Simplify the expresion
        for i in 1:getSymDim(D)
            Gradient_i[i] = simplify(Gradient_i[i])
        end
		# Add the gradient to the array
        push!(ShapeGradient, Gradient_i)
    end

    return ShapeGradient
end

# TODO: Explicar bien esto

function get_Dmat(::Type{B},::Type{M}) where {B<:AbstractBasis,M<:AbstractMaterial}
    println("Todavía no se ha implementado")
end

function get_Dmat(el::AbstractElement{D,B,M}) where {D,B,M}
    return el.Mat.DMat
end
