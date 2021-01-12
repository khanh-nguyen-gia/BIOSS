# # TODO Buscar un nombre, este es provisional,
# # y continuar con las formulas
#
# if !@isdefined(Mat2d)
#     const Mat2d = [1 -1 -1  1;
#                    1  1 -1 -1;
#                    1  1  1  1;
#                    1 -1  1 -1]
# end
#
#
# """
#     element_jacobian(element,nodal_coords)
#
#     Calculate the  Jacobian  matrix of the element
#     at the given nodal_coords.
#
#
#     The  Jacobian matrix of a trapezoidal  element
#     from a bilineal transformation can be computed
#     solving the following linear system:
#
#     ------------------------------------------------
#                 η               y        4o
#         4'o    |    o 3'        |
#                |                |               o 3
#          ------|------ ξ  <--   | 1 o
#                |                |------- x   o 2
#         1'o    |    o 2'
#
#         1 -> 1'; 2 -> 2'; 3 -> 3'; 4 -> 4'
#
#     ------------------------------------------------
# """
# function element_jacobian(el::Element{_2D,B,M}) where {B<:AbstractBasis, M<:AbstractMaterial}
#     return element_jacobian(el::Element{_2D,B,M},[0., 0.])
# end
#
# function element_jacobian(el::Element{_2D,B,M},point::Array{Float64,1}) where {B<:AbstractBasis, M<:AbstractMaterial}
#
#     a = Mat2d\el.Coords[1:4,1]
#     a2 = Mat2d\el.Coords[1:4,2]
#
#     J = zeros(2,2)
#     J[1,1]= a[2] + a[4]*point[2]
#     J[1,2]= a[3] + a[4]*point[1]
#     J[2,1]= a2[2] + a2[4]*point[2]
#     J[2,2]= a2[3] + a2[4]*point[1]
#
#     return J
# end
#
# # TODO: pensar como calcular esto de un modo mejor, para que vaya más rapidp

function element_jacobian(el::Element{_3D,B,M},point::Array{Float64,1}) where {B<:AbstractBasis, M<:AbstractMaterial}

    return "todavia nada"
end

function initJacobian(::Type{D},::Type{B}) where {D <: AbstractDim, B <: AbstractBasis}
    # obtención de la dimensión del problema
    Ndim = getSymDim(D)

    # obtención de la derivada de las funciones de forma
    Df = getShapeGrad(D,B)
    #bloque de código solución
    b = Expr(:block)

    vars = getSymVars(D)

    for i in 1:Ndim
        for j in 1:Ndim
            # creamos la expresion que almacena la suma en ij
            jac_ij = Expr(:call,:+)
            for k in 1:length(Df)
                push!(jac_ij.args,:($(Df[k][j])*el.Coords[$k,$i]))
            end
            #igualamos el jacobiano en ij a su valor
            push!(b.args,:(J[$i,$j] = $jac_ij))
        end
    end
    code = quote

    function ElementJacobian(el::Element{$D,B,M},$vars::Vector{Float64}) where {B,M}
        J = @MMatrix zeros($Ndim,$Ndim)
        $b
        return J
    end

    export ElementJacobian
    end

    eval(code)
end

function initJacobian(el::Element{D,B,M}) where {D <: AbstractDim, B<:AbstractBasis, M<:AbstractMaterial}
    return initJacobian(D,B)
end


"""
    function ElementJacobian(el,Point)

    Calculate the Jacobian matrix of an element at a given point


    Input:

        el :: Element              -> Element
        Point :: Vector{Float64}   -> Point to calculate jacobian

    Output:

        J :: Array{Float64,2}      -> Jacobian Matrix

"""
function ElementJacobian()
end

export ElementJacobian
