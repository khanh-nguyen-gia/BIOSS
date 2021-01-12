"""

    BIOSSFemDriver is the driver of the program BIOSSFem.

    It contains the functions to calculate the local stifness,
    do the Assembly and solve a given model.


"""
module BIOSSFemDriver

using BIOSSFemCore
using BIOSSFemMaterial

using BIOSSFemElements

using LinearAlgebra
using SparseArrays
using StaticArrays
using WriteVTK
using Random, Distributions

#= StiffnessModel is the file that contains the
Garlekin implementation of the Physical Ecuations
for a given physical problem

For example, the Laplace equation:

    (1)       -Δu = 0

has associated the following weak formulation of the
local stiffness

    (2)     K_el = ∫ B^T D B dx

Therefore, StiffnessModel returns the RHS of (2)
=#

include("StiffnessModel.jl")
export StiffnessModel

# File to calculate element stiffness matrix
include("ElementStiffness.jl")
export LocalStiffness
export LocalForces, calculo_fsup

# file that contains Garlekin implementarion of element mass
include("MassModel.jl")
export MassModel
# File to calculate element mass matrix
include("ElementMass.jl")
export LocalMass

include("nonlinearstifness.jl")

include("DataImport/ImportNastran.jl")
export NastranContent, ImportNastranMesh

include("FemMesh/FemMesh.jl")
# export mesh types and constructors
export AbstractMesh, RegularMesh, NastranMesh

# export Boundary Condition kinds
export AbstractBC, BCDirichlet, BCNeumann, BCForces

export BCLeft, BCRight

# export Boundary Condition struct
export BoundaryConditions
#export Boundary Condition constructor
export CreateBoundary

export Forces!
export uval!

export get_node_coordinates, get_local_coordinates, getElement, get_value, get_gradient
export distances_Matrix, adjacent4

include("DataExport/ParaviewExport.jl")
export SaveToParaview

export getKll

include("FemModel.jl")
export FemModel

include("ModelAssembly.jl")
export Assembly!, Dismantle!
export Forces, uval, getKll

export AssemblyMass!, DismantleMass!

include("NonStationaryIntegration.jl")
export vFormIntegration
export dFormIntegration

include("Solvers/LinearSolvers.jl")
export SolveProgram, SolveProgram2
export SolveProgramNonStationary
export CalcGradient

include("Solvers/AngiogenesisSolver.jl")
export SolveProgramAngiogenesis, SolveProgramAngiogenesis2
export SolverAngiogenesisAC, initial_cond_AndersonChaplain

include("PhaseField/PhaseField.jl")
#Export new TipCell constructor method
export TipCell, HypoxCell, VesselInfo
#Export Continuous behaviour functions
export initial_cond_PhaseField, time_derivative_Endothelium
export UpdateTAFmodel!, Update_TAF_evolution!
#Export miscellaneous functions
export checkifinsight, distance_nodetoTC
export Find_irrigated_HypoxCells!
#Export TipCell behaviour functions
export Update_TipCell_chemotaxis!, Update_TipCell_anastomosis!
export update_rho_TipCells!, TipCell_displacement!
export filter_nodes_PF, branching_nodes_PF, Branching!
export Merge!


# include("Solvers/NonLinearDriver.jl")
# export SolveNonLinear

end
