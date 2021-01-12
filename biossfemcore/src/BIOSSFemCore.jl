"""
    Core module of the FEMX program. It contains
    the structures and functions that will be used
    on the other modules

"""
module BIOSSFemCore

    using Calculus
    using LinearAlgebra

    using SparseArrays

    # declaration of the dimension struct
    # and asociated functions
    include("dimension.jl")

    #export the structs and functions
    export AbstractDim
    export _1D, _2D, _3D

    export getSymVars, getSymSymbol, getSymDim
    export getDimType


    # debug related files
    include("debug.jl")
    export formatMetaProg
    export initDebug


    # data in and out functions
    include("dataIO.jl")
    export dirSep
    export cacheFormat

    export DirCreation


    # sparse methods for sparse matrix
    include("SparseMethods.jl")
    export SparseMatCOO, sparse, sparsevec
    export add!, empty!, append!, isempty, Matrix


    # math algorithms to solve systems of equations
    include("math.jl")
    export SolveNewtonR
    export CalcJacobian

    # other useful algorithms and functions

    include("Other.jl")
    export findall2

end
