# functions that save the code to a cache file,
# no recalculation is needed

const cacheDir = ".cache"*dirSep*"elements"*dirSep

# Cache will be saved in ~/.julia/.cache
const homeDir = homedir()*dirSep*".julia"*dirSep

# Create general cache direrctory if not found
if !isdir(homeDir*".cache")
	@debug "No se ha encontrado el directorio de cache, creando uno"
	mkdir(homeDir*".cache")
end

#create cache directory for the elements
if !isdir(homeDir*cacheDir)
	@debug "creando cache para los elementos"
	mkdir(homeDir*cacheDir)
end

"""
	function getFileID(B,M)

	Generate the ID of the cache file for the element

	Input:
		Type{B}    -> Shape function basis
		Type{M}    -> Material
"""
function getFileID(::Type{B}, ::Type{M}) where {B <: AbstractBasis, M <: AbstractMaterial}
	return "_cache_B_$B$M.jl"
end

# renormbrar
function getFileID2(::Type{B}, ::Type{M}) where {B <: AbstractBasis, M <: AbstractMaterial}
	return "_cache_N_$B$M.jl"
end

"""
	function fileInCache(file::String)

	Check if a file is in the cache

	Input:

		file :: String

	Output:

		result :: Bool
"""
function fileInCache(file::String)
	@debug "Buscando archivo en "*homeDir*cacheDir
	return isfile(homeDir*cacheDir*file)
end


"""
	SaveToCache(file::String, code::Expr)

	Save a function (code) in the cache

	Input:

		file :: String  -> File to save the cache
		code :: Expr    -> Function to save

"""
function SaveToCache(file::String,code)
	io = open(homeDir*cacheDir*file,"w")
	println(io,code)
	close(io)
end


"""
	function ClearStoredCache()

	Deletes all the saved cache
"""
function ClearStoredCache()
	rm(homeDir*cacheDir, recursive=true)
	mkdir(homeDir*cacheDir)
end
