"""

    function findall2

    Find all indices that match a certain condicion (f)
"""
function findall2(f, a::Array{T, N}) where {T, N}
    j = 1
    b = Vector{Int}(undef, length(a))
    @inbounds for i in eachindex(a)
        b[j] = i
        j = ifelse(f(a[i]), j+1, j)
    end
    resize!(b, j-1)
    sizehint!(b, length(b))
    return b
end
