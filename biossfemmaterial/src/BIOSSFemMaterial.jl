"""

    This module contains the properties of the
    different materials

    It is divided in:

        * scalar fields: concentration, temperature

        * vector fields: displacements

    Each material contains information of
    the number of degrees of freedom of each node
    (1 scalar, 2 displacement 2D...)

"""
module BIOSSFemMaterial

using LinearAlgebra
using BIOSSFemCore
#using BIOSSFemElements
using SparseArrays # sparse arrays

# declaration of abstract type
abstract type AbstractMaterial end

export AbstractMaterial

# difusitive element

abstract type ScalarMaterial <: AbstractMaterial end

abstract type DifusionEquation <: ScalarMaterial end

export ScalarMaterial
export DifusionEquation

include("scalarfield/difusion.jl")
export Difusion, MultiDifusion
export getSpecies

include("scalarfield/convection-difusion.jl")
export ConvecDif

# disuion reaction
include("scalarfield/difusion-reaction.jl")
export DifusionReac

#thermal element
include("scalarfield/heat.jl")
export Heat

include("solid/solid.jl")
include("solid/elastic/elasto.jl")
export SolidMaterial
export Elastic1D, Elastic2D, Elastic3D, Elastic


include("solid/nonlinear/test1d.jl")
export TestNonElastic1D


# Loading angiogenesis material
include("angiogenesis/angiogenesis.jl")
export AngiogenesisMat
export CellVegf, CellVegf2, VEGF, Fibronectin
export TipCell, TipCellList, expandTipCellList!, TipCellDeactivation!, mothervessel
export Endothelium, TAF
export HypoxCell, HypoxCellDeactivation!
export getChemotacticCoefficient

include("getMaterialDOF.jl")
export getMatDof


end
