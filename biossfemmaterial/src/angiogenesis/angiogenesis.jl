"""
    AngiogenesisMat

    This material is used for angiogenesis simulation
"""
abstract type AngiogenesisMat <: AbstractMaterial end


include("continuous/cellsvegf.jl")

include("continuous/cellvegf2.jl")

include("continuous/VEGF.jl")

include("continuous/fibronectin.jl")

L0 = 1.25 # nondimensionalization length in  μm
T0 = 5460. #nondimensionalization time in seconds

include("phase-field/tipcell.jl")

include("phase-field/endothelium.jl")

include("phase-field/TAF.jl")

include("phase-field/hypox_cell.jl")
