# Mixed quad elements, 2d
struct Qmix4 <: AbstractBasis
    N :: Int
end

# Init the Qmix4 element
function Qmix4()
    return Qmix4(4)
end


struct Qmix8 <: AbstractBasis
    N :: Int
end

# Init the Qmix8 element
function Qmix8()
    return Qmix8(8)
end


struct Qmix9 <: AbstractBasis
    N :: Int
end

# Init the Qmix9 element
function Qmix9()
    return Qmix9(9)
end
