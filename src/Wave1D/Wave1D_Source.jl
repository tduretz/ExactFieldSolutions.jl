
function Wave1D_Source_u_fwd(X, params)
    x, t = X[1], X[2]
    @unpack ρ, k, σ, Gbg, G0 = params
    G = Gbg
    c = sqrt(G/ρ)
    return real(exp(im*k*(x - c*t)))
end

@doc raw"""
    sol = Wave1D_Source(x; params)  

Evaluates the manufactured solution of an initial 1D wave:

    x      : is the coordinate and time vector, x[1] is space and x[2] is time 
    params : optional parameter tuple, default: (ρ=1.0, k=8.0, σ=0.02, Gbg=1, G0=9.0) 
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Wave1D_Source([0 0])
(u = 1.0, ∇u = [0.0 -0.0], s = 0.0, G = 5.5, ρ = 1.0)
```

```julia-repl
julia> params = (ρ=10.0, k=3.0, σ=0.002, Gbg=10, G0=100.0) 
(ρ = 10.0, k = 3.0, σ = 0.002, Gbg = 10, G0 = 100.0)

julia> Wave1D_Source([0 2]; params)
(u = -0.5309925994349535, ∇u = [-2925.093952476743 -6.226919816783157], s = -3.384131664429493e8, G = 60.0, ρ = 10.0)
```
"""
function Wave1D_Source(X;
    params = (ρ=1.0, k=8.0, σ=0.2, Gbg=1, G0=9.0) )
    x, t   = X[1], X[2]
    u      = Wave1D_Source_u_fwd(X, params)
    f_cl   = X -> Wave1D_Source_u_fwd(X, params)
    gradu  = ForwardDiff.gradient(f_cl, X)
    hessu  = ForwardDiff.hessian(f_cl, X)

    ∂σ∂x   = params.Gbg*gradu[1]
    s      = 0.01*exp(-x.^2/params.σ^2)*exp(-t*100)

    G      =  params.Gbg
    # hessu  = ForwardDiff.hessian(f_cl, X)
    # s      = hessu[2,2] - params.Gbg/params.ρ*(hessu[1,1])
    # @show s
    return (u=u, ∇u=gradu, s=s, G=G, ρ=params.ρ)
end