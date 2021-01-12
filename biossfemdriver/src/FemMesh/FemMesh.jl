# Definition of the AbstractMesh abstract type

abstract type AbstractMesh{D <: AbstractDim} end


# Definition of the RegularMesh type
"""
    struct RegularMesh

    Struct that contains the information of a regular mesh

    DimNodes :: Array{Int,1}                     -> Nodes in each dimension

    EndPoint :: Array{Float64,1}                 -> Coordinates of the corner

    Elements :: Array{AbstractElement,1}         -> Array that contains the elements

    NofNodes :: Int                              -> number of nodes
    DOF :: Int                                   -> number od degrees of freedom
    Nodes :: Array{Int}                          -> Array with node disposition

    Coordinates :: Array{Float64}
    Dim :: D
"""
mutable struct RegularMesh{D} <: AbstractMesh{D}
    DimNodes :: Array{Int,1}                     # Nodes in each dimension

    EndPoint :: Array{Float64,1}                 # Coordinates of the corner

    Elements :: Array{AbstractElement,1}         # Array that contains the elements

    NofNodes :: Int                              # number of nodes
    DOF :: Int                                   # number od degrees of freedom
    Nodes :: Array{Int}                          # Array with node disposition

    Coordinates :: Array{Float64}
    Dim :: D
end

# regular mesh constructor
include("RegularMesh.jl")

"""
    struct NastranMesh

    Struct that contains the information of a Fem mesh
    imported from Nastran

    Elements :: Array{AbstractElement,1}         -> Array that contains the elements

    NofNodes :: Int                              -> number of nodes
    DOF :: Int                                   -> number od degrees of freedom
    Nodes :: Array{Int}                          -> Array with node disposition

    Coordinates :: Array{Float64}
    Dim :: D
"""
mutable struct NastranMesh{D} <: AbstractMesh{D}

    Elements :: Array{AbstractElement,1}         # Array that contains the elements

    NofNodes :: Int                              # number of nodes
    DOF :: Int                                   # number od degrees of freedom
    Nodes :: Array{Int}                          # Array with node disposition

    Coordinates :: Array{Float64}
    Dim :: D
end

# NastranMesh Constructor
include("NastranMesh.jl")

include("BoundaryConditions.jl")
