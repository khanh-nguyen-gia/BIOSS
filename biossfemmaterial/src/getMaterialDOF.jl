"""
    This file contains the number of degrees of freedom
    for each material

    It is a map that returns the material DOF
"""
getMatDof = Dict{DataType,Int}()

getMatDof[Difusion] = 1
for Nindex in 1:10
    getMatDof[MultiDifusion{Nindex}] = Nindex
end


getMatDof[Heat] = 1
getMatDof[ConvecDif] = 1
getMatDof[DifusionReac] = 1

getMatDof[Elastic1D] = 1
getMatDof[Elastic2D] = 2
getMatDof[Elastic3D] = 3

getMatDof[TestNonElastic1D] = 1

# concentracion de celulas y de vegf
getMatDof[CellVegf] = 2
getMatDof[CellVegf2] = 1
getMatDof[VEGF] = 1
getMatDof[Fibronectin] = 1
getMatDof[Endothelium] = 1
getMatDof[TAF] = 1
