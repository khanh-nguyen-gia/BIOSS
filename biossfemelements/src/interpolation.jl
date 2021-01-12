"""
    Calculates gradient at a given point
"""
function CalcB(el::Element, u, point)
    J = ElementJacobian(el, point)
    Bmat = getBmat(el, inv(J), point)
    
    B0 = Bmat*u

    return B0

end
