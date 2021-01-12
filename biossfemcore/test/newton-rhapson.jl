using Revise
using BIOSSFemCore


using Calculus



function f1(x::Vector{Float64})
    return cos(x[1]) + cos(x[2])
end



function f2(x::Vector{Float64})
    return x[1] + x[2]
end

function F(x::Vector{Float64})
    F = zeros(2)

    F[1] = cos(x[1]) + cos(x[2])
    F[2] = x[1] + x[2]

    return F
end


function Dfuns(x::Vector{Float64})

    DF = zeros(2,2)

    DF[1,:] = [-sin(x[1]),-sin(x[2])]
    DF[2,:] = [1.,1.]

    return DF
end

F([3.,3.])
Dfuns([3.,3.])

println(supertype(typeof(Dfuns)))

@time x = SolveNewtonR(F, Dfuns,[3., 2.])
@time xp = SolveNewtonR(F,[3., 2.])


function F2(x::Vector{Float64})
    F = zeros(3)

    F[1] = 3x[1] - cos(x[2]*x[3]) -1/2
    F[2] = x[1]^2 - 81*(x[2]+.1)^2  + sin(x[3])+ 1.06
    F[3] = exp(-x[1]*x[2]) + 20*x[3] + (10π -3)/3.

    return F
end

@time xp = SolveNewtonR(F2,[0.1, 0.1,-0.1])
println(xp,F2(xp))

function Res(x::Float64,P::Float64)

    F = P - x/cos(x)

    return F
end

Pes = 0:.5:30
x = zeros(length(Pes))

@time for i in 2:length(Pes)
    global x, Pes
    x[i] = SolveNewtonR(x -> Res(x,float(Pes[i])),x[i-1])
end

p= plot(x,Pes)
display(p)
