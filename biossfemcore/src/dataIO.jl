# Check te OS to add dir sep

if Sys.iswindows()
	const dirSep= "\\"
elseif Sys.isapple() | Sys.islinux()
	const dirSep = "/"
end



"""
	function cacheFormat(code)

	Transform machine-generated code to
	a more readable format, used to save it
	in the cache directory
"""
function cacheFormat(code)

    tempStr = string(code)
    tempStr = split(tempStr,"\n")

    formatedCode = ""

    for i in 2:length(tempStr)-1
        if !occursin("#",tempStr[i])
		temp= replace(tempStr[i], "\n" => "")
		formatedCode *= temp[5:end]*"\n "
        end
    end

    return formatedCode
end



"""
	function DirCreation(dir; replace)

	This function creates a directory, used for
	saving files later. If the directory exists,
	there are two options:

		# replace = false -> creates new directory
		with different name (ej: dir2)

		# replace = true -> deletes the previoues

	The directory created is relative to user's homedir

"""
function DirCreation(dir::String; replace = false, append=false)

	NDirMax = 10000 # maximun number of directories

	CurrentDir = homedir()*dirSep*dir #absolute dir path

	# check replace variable
	if replace
		# remove directory if exists

		if ispath(CurrentDir)
			if append
				return
			else
	        	rm(CurrentDir, recursive=true)
				mkdir(CurrentDir)
			end
	    else
			# create directory if it doesnt exist
		    mkdir(CurrentDir)
			return
		end

	else
		# check if directory exists
		i = 2
		if ispath(CurrentDir)
			# iterate through all the directories.
			while i < NDirMax
				if ispath(CurrentDir*"$i")
					# if directory i exists, go to i+1
					i+= 1
				else
					# chreate directory i if not exists
					mkdir(CurrentDir*"$i")
					return
				end
			end

			println("Maximum number of directories reached")
		else
			# if it does not exists, create dir
		    mkdir(CurrentDir)
			return

		end
	end
end
