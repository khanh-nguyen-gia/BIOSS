"""
    function SaveToParaview(Mesh,u,name,fname,dir)

    Save results to a paraview-compatible format file

    # if replace, it will delete previews directory

    # if append, it will add

"""
function SaveToParaview(Mesh::RegularMesh{D},
    u::Array, name::Vector,
    fname::String, dir::String; replace = true, append = true) where {D <: AbstractDim}


    LocalDir = homedir()*dirSep*dir
    # create directory if not exists
    DirCreation(dir; replace=replace, append=append)

    vtkfile = vtk_grid(fname, Mesh.Coordinates)

    for i in 1:length(u)
        vtkfile[name[i]] = u[i]
    end

    prev_dir = pwd()

    cd(LocalDir)

    outfiles = vtk_save(vtkfile)

    cd(prev_dir)

end

function SaveToParaview(Mesh::RegularMesh{D},fname,dir) where {D <: AbstractDim}
    SaveToParaview(Mesh, [vec(Mesh.Nodes)], ["Nodos"], fname, dir)
end


function OutputDisp(Mesh::RegularMesh{D},u,name,fname,dir) where {D <: AbstractDim}
    println("Exportando tus datos")

    CDef = deepcopy(Mesh.Coordinates)

    for I in CartesianIndices(Mesh.Coordinates[1,:,:])
        i = LinearIndices(Mesh.Coordinates[1,:,:])[I]
        CDef[:,I] += [u[1][i], u[2][i]]
    end


    vtkfile = vtk_grid(fname, CDef)  # 2-D
    for i in 1:length(u)
        vtkfile[name[i]] = u[i]
    end

    cd(dir)

    outfiles = vtk_save(vtkfile)
end


function OutputDisp3D(Mesh::RegularMesh{D},u,name,fname,dir) where {D <: AbstractDim}

    if isfile(fname)
        rm(fname)
    end


    CDef = deepcopy(Mesh.Coordinates)

    for I in CartesianIndices(Mesh.Coordinates[1,:,:,:])
        i = LinearIndices(Mesh.Coordinates[1,:,:,:])[I]
        CDef[:,I] += [u[1][i], u[2][i], u[3][i]]
    end


    vtkfile = vtk_grid(fname, CDef)  # 2-D
    for i in 1:length(u)
        vtkfile[name[i]] = u[i]
    end

    prev_dir = pwd()

    cd(dir)

    outfiles = vtk_save(vtkfile)

    println("Exportando tus datos")

    cd(prev_dir)

end

export OutputDisp, OutputDisp3D
