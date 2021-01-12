"""
    BIOSSFemElements module. It contains the defininition
    of the element structure.

    Each element has a given dimension, material
    and a basis, which detemines the shape functions
    that it uses.

"""
module BIOSSFemElements

# Load Julia modules:

using Calculus              # Symbolic Calculus
using LinearAlgebra         # Inverse and determinant

using StaticArrays          # static arrays, faster performance
using SparseArrays          # sparse matrix store


# Load BIOSS dependencies

using BIOSSFemCore
using BIOSSFemMaterial

# Supertype of the basis struct
abstract type AbstractBasis end

export AbstractBasis


# Load general basis properties
include("BasisProperties.jl")

export get_basis_gdl

# Load the basis for 1D, 2D, 3D elements

# 1D Element basis
include("El1D/LinearElement.jl")

export Lin2, Lin3

# 2D Element basis

include("El2D/QuadElement.jl")

export Quad4, Quad8, Quad9

include("El2D/QmixElement.jl")

export Qmix4, Qmix8, Qmix9

# 3D Element basis

include("El3D/BrickElement.jl")
export Brick8

# Definition of the element supertype
# The element has 3 types asociated: Dimension, Basis, Material

abstract type AbstractElement{D <: AbstractDim, B <: AbstractBasis, M <: AbstractMaterial} end

export AbstractElement


# Creation of the element struct
"""
    struct  Element{D,B,M} <: AbstractElement{D,B,M}

    Id :: Int                            -> ID of the element
    Nodes :: Array{Int,1}                -> Correspondence with global nodes
    Coords :: Array{Float64,2}           -> Global coordinates
    IntP :: Int                          -> Number of integration points
    ndof :: Int                          -> Element degrees of freedom

    Dim :: D                             -> Dimensions of the element
    Basis :: B                           -> Basis of the element
    Mat :: M                             -> Material of the element

"""
mutable struct Element{D,B,M} <: AbstractElement{D,B,M}  # -----------------------
    Id :: Int                                          # ID del elemento
    Nodes :: Array{Int,1}                              # Nodos globales
    Coords :: Array{Float64,2}                         # Coordenadas globales
    IntP :: Int                                        # Number of integration points
    ndof :: Int                                        # element degrees of freedom

    Dim :: D                                           #
    Basis :: B                                         # Tipo de elemento
    Mat :: M                                           # Tipo de material
end                                                    # -----------------------

# Load element constructors #

# 1D Element
include("El1D/Element1D.jl")

# 2D Element
include("El2D/Element2D.jl")

# 3D Element

include("El3D/Element3D.jl")

export Element


# Load functions to save the shape functions
include("Cache.jl")

export ClearStoredCache

include("ElementTasks.jl")

export getElementDim, get_basis_props, get_int_points

export SpatialNodesToIndices


# Load the constructor for the shape functions
include("ShapeFunctions.jl")
export getShape, getShapeGrad

include("Bmat.jl")

export init_getBmat

include("Nmat.jl")

export init_getNmat

include("interpolation.jl")
export CalcB


# init_getBmat(_2D,Quad4, Difusion)
# init_getBmat(_2D,Quad9, Difusion)
# export getBmat

export get_Dmat

include("jacobians.jl")
export initJacobian
export element_jacobian


export homeDir, cacheDir

end
