function Element(::Type{_1D}, ::Type{B}, ::Type{M}, Id::Int,
                   nodes::Array{Int,1}) where {B <:AbstractBasis, M<:AbstractMaterial}

    Dim = _1D()
    Basis = B()
    Mat = M(_1D)

    N = get_basis_gdl(B)

    return Element(Id,nodes,zeros(N,1),N,getMatDof[M],Dim,Basis,Mat)
end


function Element(::Type{_1D}, ::Type{B}, ::Type{M},
                   nodes::Array{Int,1}) where {B<:AbstractBasis, M<:AbstractMaterial}


    Dim = _1D()
    Basis = B()
    Mat = M(_1D)

    N = get_basis_gdl(B)

    return Element(-1,nodes,zeros(N,1),N,getMatDof[M],Dim,Basis,Mat)
end

function Element(::Type{_1D},::Type{B}, Mat::M,
                   nodes::Array{Int,1}) where {B<:AbstractBasis, M<:AbstractMaterial}

    Dim = _1D()
    Basis = B()

    N = get_basis_gdl(B)

    return Element(-1,nodes,zeros(N,1),N,getMatDof[M],Dim,Basis,Mat)
end
