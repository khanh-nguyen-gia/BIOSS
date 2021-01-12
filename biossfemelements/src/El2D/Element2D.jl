"""

    function Element(D,T,M,nodes)

    Creates a 2d element of the type T, and material M,
    which global nodes correspond to the variable nodes.

    Default settings:

        4 integration points for Quad4
        9 integration points for Quad8, Quad9

"""
function Element(::Type{_2D}, ::Type{B}, ::Type{M}, Id::Int,
                   nodes::Array{Int,1}) where {B <:AbstractBasis, M<:AbstractMaterial}

    Dim = _2D()
    Basis = B()
    Mat = M(_2D)

    N = get_basis_gdl(B)

    return Element(Id,nodes,zeros(N,2),N,getMatDof[M],Dim,Basis,Mat)
end

function Element(::Type{_2D}, ::Type{B}, ::Type{M},
                   nodes::Array{Int,1}) where {B<:AbstractBasis, M<:AbstractMaterial}


    Dim = _2D()
    Basis = B()
    Mat = M(_2D)

    N = get_basis_gdl(B)

    return Element(-1,nodes,zeros(N,2),N,getMatDof[M],Dim,Basis,Mat)
end

function Element(::Type{_2D},::Type{B}, Mat::M,
                   nodes::Array{Int,1}) where {B<:AbstractBasis, M<:AbstractMaterial}

    Basis = B()
    Dim = _2D()

    N = get_basis_gdl(B)

    return Element(-1,nodes,zeros(N,2),N,getMatDof[M],Dim,Basis,Mat)
end
