module BIOSSFem

using BIOSSFemCore

export AbstractDim
export _1D, _2D, _3D


using BIOSSFemMaterial

export Difusion, Heat

using BIOSSFemElements

export Quad4, Quad8, Quad9
export Qmix4, Qmix8, Qmix9
export Element

export init_getBmat, initJacobian, init_getNmat


using BIOSSFemDriver

export RegularMesh
export CreateBoundary
export FemModel
export SolveProgram

end # module
