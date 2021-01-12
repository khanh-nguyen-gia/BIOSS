function init_getBmat(::Type{D}, ::Type{B},
	::Type{M}) where {D <: AbstractDim, B<:AbstractBasis,M<:AbstractMaterial}

	file = getFileID(B,M)

    if fileInCache(file)

		@debug "Encontrada la funcion en el cache"
    	include(homeDir*cacheDir*file)
		eval(:(export getBmat))
		return

	else

		@debug "No encontrada la funcion"

		b = SymbolicBmat(D,B,M)

		# creation of the file that constains the function

		vars = getSymVars(D)
    	code = quote
        function getBmat(::Type{$D},::Type{$B},::Type{$M},J_1,$vars)
            $b
            return Bmat
        end

        function getBmat(el::Element{D,B,M},J_1,x) where {D,B,M}
            $(vars) = x
            return getBmat(D,B,M,J_1,$vars)
        end

		export getBmat
    	end

    	@debug formatMetaProg("This is the generated shape funcion",code)
    	@debug "Saving file to cache"*file

    	SaveToCache(file, code)

    	eval(code)

    end
end

function getBmat()
end

export getBmat

#obtain the cinematic matrix for a given material
function SymbolicBmat(::Type{D}, ::Type{B},
	::Type{M}) where {D <: AbstractDim, B<:AbstractBasis,M<:ScalarMaterial}

	Dim = getSymDim(D)

	dfuns = getShapeGrad(D,B)
	b = Expr(:block)
	push!(b.args,:(Bmat = zeros($(getSymDim(D)),$(get_basis_gdl(B)))))

	# evaluation of symbolic function

	for i in 1:get_basis_gdl(B)
		for j in 1:Dim
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$j]))
			end
			push!(b.args,:(Bmat[$j,$i] = $RowMulti))
		end
	end

	return b
end



function SymbolicBmat(::Type{_1D}, ::Type{B},
	::Type{M}) where {B<:AbstractBasis,M<:SolidMaterial}

	Dim = getSymDim(_1D)

	dfuns = getShapeGrad(_1D,B)
	b = Expr(:block)
	push!(b.args,:(Bmat = zeros($(1),$(get_basis_gdl(B)))))

	# evaluation of symbolic function

	for i in 1:get_basis_gdl(B)
		# inicio del sumatorio
		for j in 1:1
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$j]))
			end
		push!(b.args,:(Bmat[$j,$i] = $RowMulti))
		end
	end


	return b
end


function SymbolicBmat(::Type{_2D}, ::Type{B},
	::Type{M}) where {B<:AbstractBasis,M<:SolidMaterial}

	Dim = getSymDim(_2D)

	dfuns = getShapeGrad(_2D,B)
	b = Expr(:block)
	push!(b.args,:(Bmat = zeros($(3),$(get_basis_gdl(B)*2))))

	# evaluation of symbolic function

	for i in 1:get_basis_gdl(B)
		# inicio del sumatorio
		for j in 1:2
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$j]))
			end
		push!(b.args,:(Bmat[$j,$(2*(i-1)+j)] = $RowMulti))
		end
	end
	for i in 1:get_basis_gdl(B)
		# inicio del sumatorio
		for j in 2:-1:1
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$(3-j)]))
			end
		push!(b.args,:(Bmat[3,$(2*(i-1)+j)] = $RowMulti))
		end
	end

	return b
end

function SymbolicBmat(::Type{_3D}, ::Type{B},
	::Type{M}) where  {B<:AbstractBasis,M<:SolidMaterial}

	Dim = getSymDim(_3D)

	dfuns = getShapeGrad(_3D,B)
	b = Expr(:block)
	push!(b.args,:(Bmat = @MMatrix zeros($(6),$(get_basis_gdl(B)*3))))
	#push!(b.args,:(Bmat = @MMatrix zeros($(6),$(get_basis_gdl(B)*3))))
	# evaluation of symbolic function

	for i in 1:get_basis_gdl(B)
		# inicio del sumatorio
		for j in 1:3
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$j]))
			end
		push!(b.args,:(Bmat[$j,$(3*(i-1)+j)] = $RowMulti))
		end
	end
	for i in 1:get_basis_gdl(B)
		# inicio del sumatorio
		for j in [[2,3], [3,2]]
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$(j[2])]))
			end
		push!(b.args,:(Bmat[4,$(3*(i-1)+j[1])] = $RowMulti))
		end

		for j in  [[1,3], [3,1]]
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$(j[2])]))
			end
		push!(b.args,:(Bmat[5,$(3*(i-1)+j[1])] = $RowMulti))
		end

		for j in  [[1,2], [2,1]]
			RowMulti = Expr(:call,:+)
			for k in 1:Dim
				push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$(j[2])]))
			end
		push!(b.args,:(Bmat[6,$(3*(i-1)+j[1])] = $RowMulti))
		end
	end

	return b
end

#obtain the cinematic matrix for MultiDifusion material
function SymbolicBmat(::Type{D}, ::Type{B},
	::Type{MultiDifusion{N}}) where {D <: AbstractDim, B<:AbstractBasis,N}

	# Number of spatial dimension
	NDim = getSymDim(D)

	# Number of substancess
	NSpecies = getSpecies(MultiDifusion{N})

	# gradient of shape functions
	dfuns = getShapeGrad(D,B)
	b = Expr(:block)

	#push!(b.args,:(Bmat = zeros($(NDim*NSpecies),$(get_basis_gdl(B)*NSpecies))))
	push!(b.args,:(Bmat = spzeros($(NDim*NSpecies),$(get_basis_gdl(B)*NSpecies))))
	# evaluation of symbolic function
	for m in 1:NSpecies
		for i in 1:get_basis_gdl(B)
				for j in 1:NDim
					RowMulti = Expr(:call,:+)
					for k in 1:NDim
						push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$j]))
					end
					I = NSpecies*(i-1)+m
					J = (NDim)*(m-1)+j
					#println("J = $J, I = $I")
					push!(b.args,:(Bmat[$J,$I] = $RowMulti))
				end
			end
	end

	return b
end

# Obtain B matrix for CellVegf material
function SymbolicBmat(::Type{D}, ::Type{B},
	::Type{M}) where {D <: AbstractDim, B<:AbstractBasis, M<:AngiogenesisMat}

	# Number of spatial dimension
	NDim = getSymDim(D)

	# Number of substancess
	NSpecies = getMatDof[M]

	# gradient of shape functions
	dfuns = getShapeGrad(D,B)
	b = Expr(:block)

	#push!(b.args,:(Bmat = zeros($(NDim*NSpecies),$(get_basis_gdl(B)*NSpecies))))
	push!(b.args,:(Bmat = spzeros($(NDim*NSpecies),$(get_basis_gdl(B)*NSpecies))))
	# evaluation of symbolic function
	for m in 1:NSpecies
		for i in 1:get_basis_gdl(B)
				for j in 1:NDim
					RowMulti = Expr(:call,:+)
					for k in 1:NDim
						push!(RowMulti.args,:($(dfuns[i][k])*J_1[$k,$j]))
					end
					I = NSpecies*(i-1)+m
					J = (NDim)*(m-1)+j
					#println("J = $J, I = $I")
					push!(b.args,:(Bmat[$J,$I] = $RowMulti))
				end
			end
	end

	return b
end
