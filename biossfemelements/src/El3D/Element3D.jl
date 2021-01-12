# Declaration of 3d Element


function Element(::Type{_3D}, ::Type{B}, ::Type{M}, Id::Int,
                   nodes::Array{Int,1}) where {B <:AbstractBasis, M<:AbstractMaterial}

    Dim = _3D()
    Basis = B()
    Mat = M(_3D)

    N = get_basis_gdl(B)

    return Element(Id,nodes,zeros(N,3),N,getMatDof[M],Dim,Basis,Mat)
end

function Element(::Type{_3D}, ::Type{B}, ::Type{M},
                   nodes::Array{Int,1}) where {B<:AbstractBasis, M<:AbstractMaterial}


    Dim = _3D()
    Basis = B()
    Mat = M(_3D)

    N = get_basis_gdl(B)

    return Element(-1,nodes,zeros(N,3),N,getMatDof[M],Dim,Basis,Mat)
end

function Element(::Type{_3D},::Type{B}, Mat::M,
                   nodes::Array{Int,1}) where {B<:AbstractBasis, M<:AbstractMaterial}

    Basis = B()
    Dim = _3D()

    N = get_basis_gdl(B)

    return Element(-1,nodes,zeros(N,3),N,getMatDof[M],Dim,Basis,Mat)
end
