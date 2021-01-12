function init_getNmat(::Type{D}, ::Type{B},
	::Type{M}) where {D <: AbstractDim, B<:AbstractBasis,M<:AbstractMaterial}

	file = getFileID2(B,M)

    if fileInCache(file)

		@debug "Encontrada la funcion en el cache"
    	include(homeDir*cacheDir*file)
		eval(:(export getNmat))
		return

	else

		@debug "No encontrada la funcion"

		b = SymbolicNmat(D,B,M)

		# creation of the file that contains the function

		vars = getSymVars(D)
    	code = quote
        function getNmat(::Type{$D},::Type{$B},::Type{$M},$vars)
            $b
            return Nmat
        end

        function getNmat(el::Element{D,B,M},x) where {D,B,M}
            $(vars) = x
            return getNmat(D,B,M,$vars)
        end

		export getNmat
    	end

    	@debug formatMetaProg("This is the generated shape funcion",code)
    	@debug "Saving file to cache"*file

    	SaveToCache(file, cacheFormat(code))

    	eval(code)

    end
end

#obtain the nodal matrix for a given material
function SymbolicNmat(::Type{D}, ::Type{B},
	::Type{M}) where {D <: AbstractDim, B<:AbstractBasis,M<:ScalarMaterial}

	Dim = getSymDim(D)

	funs = getShape(D,B)
	b = Expr(:block)
	push!(b.args,:(Nmat = zeros(1,$(get_basis_gdl(B)))))

	# evaluation of symbolic function

	for i in 1:get_basis_gdl(B)
		RowMulti = Expr(:call,:+)
		push!(b.args,:(Nmat[$i] = $(funs[i])))
	end

	return b
end

function SymbolicNmat(::Type{D}, ::Type{B},
	::Type{MultiDifusion{N}}) where {D <: AbstractDim, B<:AbstractBasis, N}

	# Number of spatial dimension
	NDim = getSymDim(D)
	# Number of substancess
	NSpecies = getSpecies(MultiDifusion{N})

	funs = getShape(D,B)
	b = Expr(:block)
	push!(b.args,:(Nmat = zeros($NSpecies,$(get_basis_gdl(B)*NSpecies))))

	# evaluation of symbolic function

	for i in 1:get_basis_gdl(B)
		for m in 1:NSpecies
			RowMulti = Expr(:call,:+)
			push!(b.args,:(Nmat[$(m),$((i-1)*NSpecies+m)] = $(funs[i])))
		end
	end

	return b
end

function SymbolicNmat(::Type{D}, ::Type{B},
	::Type{M}) where {D <: AbstractDim, B<:AbstractBasis, M<:AngiogenesisMat}


	# Number of spatial dimension
	NDim = getSymDim(D)
	# Number of substancess
	NSpecies = getMatDof[M]

	funs = getShape(D,B)
	b = Expr(:block)
	push!(b.args,:(Nmat = zeros($NSpecies,$(get_basis_gdl(B)*NSpecies))))

	# evaluation of symbolic function

	for i in 1:get_basis_gdl(B)
		for m in 1:NSpecies
			RowMulti = Expr(:call,:+)
			push!(b.args,:(Nmat[$(m),$((i-1)*NSpecies+m)] = $(funs[i])))
		end
	end

	return b
end


function getNmat()
end

export getNmat
