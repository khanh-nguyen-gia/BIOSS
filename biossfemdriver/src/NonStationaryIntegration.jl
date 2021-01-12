"""
    function vFormIntegration(model, dt, f, d0, v0; alpha)

    This function performs the integration of a time step
    in the velocity form.

    Implementation from Hughes,  page 460

    <
            The finite element method:
    linear static and dynamic finite element analysis

    >

    Input:

        model :: FemModel           -> Fem model to solve
        dt :: Float64            -> Time step
        f :: Vector{Float64}     -> Forces
        d0 ::Vector{Float64}     -> Initial displacement
        v0 ::Vector{Float64}     -> Initial velocity

    Optional Input:

        alpha :: Float64         -> Parameter to select time-integration method,
                                    alpha in [0,1]. Default: alpha = 1

        alpha   |          method
        --------------------------
        0            Forward Euler
        1/2         Crank-Nicolson
        1           Backward Euler

    Output:

    d1 ::Vector{Float64}     -> Final displacement
    v1 ::Vector{Float64}     -> Final velocity
"""
function vFormIntegration(model::FemModel, dt::Float64,
    f::Vector{Float64}, d0::Vector{Float64},
    v0::Vector{Float64}; alpha=0.)

    v1 = copy(v0)
    d1 = copy(d0)

    predictor = d0[model.FreeNodes] + (1. - alpha)*dt.*v0[model.FreeNodes]

    if alpha != 0.
        Kaux  = sparse(model.M) + alpha*dt.* sparse(model.K)
    else
        Kaux  = sparse(model.M)
    end

    faux = f[model.FreeNodes] - sparse(model.K)[model.FreeNodes,model.FreeNodes]*predictor

    v1[model.FreeNodes] = Kaux[model.FreeNodes,model.FreeNodes] \ faux

    d1[model.FreeNodes] = predictor + alpha*dt.*v1[model.FreeNodes]

    return d1, v1
end

"""
    function dFormIntegration(model, dt, f, d0, v0; alpha)

    This function performs the integration of a time step
    in the displacement form.

    Implementation from Hughes,  page 460

    <
            The finite element method:
    linear static and dynamic finite element analysis

    >

    Input:

        model :: FemModel           -> Fem model to solve
        dt :: Float64            -> Time step
        f :: Vector{Float64}     -> Forces
        d0 ::Vector{Float64}     -> Initial displacement
        v0 ::Vector{Float64}     -> Initial velocity

    Optional Input:

        alpha :: Float64         -> Parameter to select time-integration method,
                                    alpha in [0,1]. Default: alpha = 1

        alpha   |          method
        --------------------------
        0            Forward Euler
        1/2         Crank-Nicolson
        1           Backward Euler

    Output:

    d1 ::Vector{Float64}     -> Final displacement
    v1 ::Vector{Float64}     -> Final velocity
"""
function dFormIntegration(model::FemModel, dt::Float64,
    f::Vector{Float64}, d0::Vector{Float64},
    v0::Vector{Float64}; alpha=1.)

    v1 = copy(v0)
    d1 = copy(d0)

    predictor = d0[model.FreeNodes] + (1. - alpha)*dt.*v0[model.FreeNodes]

    Kaux  = (sparse(model.M) + alpha*dt.* sparse(model.K))./(alpha*dt)

    faux = f[model.FreeNodes] + sparse(model.M)[model.FreeNodes,model.FreeNodes]./(alpha*dt)*predictor

    d1[model.FreeNodes] = Kaux[model.FreeNodes,model.FreeNodes] \ faux

    v1[model.FreeNodes] = (d1[model.FreeNodes] - predictor)./(alpha*dt)

    return d1, v1
end
