#=
    Math.jl


    This folder contains various functions to solve
    systems of equations among other complementary functions



=#

"""
    function SolveNewtonR(F::Fun1, DF::Fun2,
                            x0::Vector{Float64})

    Compute the solution of the system of equations F, where the
    gradient is given in DF. Returns the solution and the number
    of iterations

    Input:

        F :: Function           -> System of equations to solve
        F2 :: Function          -> Gradient of the system
        x0 :: Vector{Float64}   -> Initial estimate

    Optional Input:

        MaxIter :: Int          -> Maximum number of iterations
        tol :: Float64          -> Tolerance of the solution

    Output:
        sol :: Vector{Float64}  -> Solution of the system
        i :: Int                -> Number of iterations
"""
function SolveNewtonR(F::Fun1, DF::Fun2,
    x0::Vector{Float64};
    MaxIter = 1e4,
    tol = 1e-6) where {Fun1 <: Function, Fun2 <: Function}
    sol = x0.*0

    for i in 0:Int(MaxIter)

        sol = x0 - inv(DF(x0))*F(x0)

        if norm(x0 - sol) < tol
            return sol, i
        end

        x0 = sol
    end
end


function CalcJacobian(F::Fun1,x::Vector{Float64}) where {Fun1 <: Function}
    ndims = length(x)
    funs::Vector{Function} = []
    for i in 1:ndims
        push!(funs,x -> F(x)[1])
    end
    return funs
end


"""
    function SolveNewtonR(F::Fun1, x0::Vector{Float64})

    Compute the solution of the system of equations.
    The gradient is calculated numerically.
    Returns the solution and the number of iterations

    Input:

        F :: Function           -> System of equations to solve
        x0 :: Vector{Float64}   -> Initial estimate

    Optional Input:

        MaxIter :: Int          -> Maximum number of iterations
        tol :: Float64          -> Tolerance of the solution

    Output:
        sol :: Vector{Float64}  -> Solution of the system
        i :: Int                -> Number of iterations
"""
function SolveNewtonR(F::Fun1,
    x0::Vector{Float64};
    MaxIter = 1e2,
    tol = 1e-6) where {Fun1 <: Function}
    sol = x0.*0

    ndims = length(x0)
    DF = zeros(ndims,ndims)

    for i in 1:Int(MaxIter)

        for k in 1:ndims
            DF[k,:] = Calculus.gradient(x -> F(x)[k],x0)
        end

        sol = x0 - DF \ F(x0)

        #sol = x0 - inv(DF)*F(x0)

        if norm(x0 - sol) < tol
            return sol, i
        end

        # println("Iteración: $i, error: $(norm(x0 - sol))")

        x0 = copy(sol)
    end

    return sol, MaxIter


end


function SolveNewtonR(F::Fun1,
    x0::Float64;
    MaxIter = 1e6,
    tol = 1e-6) where {Fun1 <: Function}
    sol = x0.*0

    ndims = length(x0)


    for i in 1:MaxIter

        DF = Calculus.gradient(F,x0)

        sol = x0 - F(x0)/DF

        if norm(x0 - sol) < tol
            break
        end

        x0 = sol
    end

    return sol
end
