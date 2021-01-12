using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver



using TimerOutputs

import BIOSSFemDriver:CreateElements, getElementNodes, getCorner

N = 20

to = TimerOutput()

function CreateElements(::Type{D},::Type{B},::Type{M}, DimNodes, Nodes,
    Coordinates) where {D <: AbstractDim, B <: AbstractBasis, M<:AbstractMaterial}

    Els::Array{AbstractElement,1} = []
    NofElements = prod(DimNodes)
    ElementTable = collect(reshape(1:NofElements, DimNodes...))

    for I in CartesianIndices(ElementTable)
        @timeit to "get nodes" elNodes = getElementNodes(I.I,DimNodes)
        element_i = Element(D,B,M,ElementTable[I],elNodes)
        for k in 1:length(elNodes)
            index = CartesianIndices(Nodes)[elNodes[k]]
            element_i.Coords[k,:] = Coordinates[:,CartesianIndices(Nodes)[elNodes[k]]]
        end
        push!(Els, element_i)
    end
    return Els
end

function getElementNodes(I,DimNodes)
    Ndim = length(DimNodes)
    dx = zeros(Int,Ndim)
    for i in 1:Ndim
        dx[i] = prod(DimNodes[1:i-1].+1)
    end

    Corner = getCorner(Ndim)
    elNodes = zeros(Int,2^Ndim)

    for i in 1:2^Ndim
        elNodes[i] = (Corner[i] + collect(I).-1 )'*dx + 1
    end

    return elNodes
end

function getCorner(N::Int)
    if N > 3
        PrevCorner = getCorner(N-1)
        return vcat(push!.(deepcopy(PrevCorner),0),push!.(deepcopy(PrevCorner),1))
    elseif N ==3
        return [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    elseif N == 2
        return [[0, 0], [1, 0], [1, 1], [0, 1]]
    elseif N == 1
        return [[0], [1]]
    end
end

@timeit to "Mesh creation" mesh = RegularMesh(_3D,Brick8,Difusion, [N, N, N], [1., 1., 1.])

print_timer(to)
