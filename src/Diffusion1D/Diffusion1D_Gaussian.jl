
function Diffusion1D_Gaussian_u_fwd(x, params)
    @unpack T0, K, σ  = params
    return T0 / sqrt(1 + 2*x[2]*K/σ^2) * exp.(-x[1]^2 / (2*(σ^2 + 2*x[2]*K)) )
end

@doc raw"""
    sol = Diffusion1D_Gaussian(x; params)  

Evaluates the manufactured solution of an initial 1D Gaussian anomaly:

    x      : is the coordinate and time vector, x[1] is space and x[2] is time 
    params : optional parameter array, default: (T0 = 100., K = 1e-6, σ = 0.1 ) 
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Diffusion1D_Gaussian([0 0])
(u = 100.0, ∇u = [0.0 0.0 -0.06283185307179585], s = 0.0)
```

```julia-repl
julia> params = (T0 = 1., K = 1.0, σ = 0.1 ) 
(T0 = 1.0, K = 1.0, σ = 0.1)

julia> Diffusion1D_Gaussian([0 0])
(u = 1.0, ∇u = [0.0 -99.99999999999999], s = 0.0)
```
"""
function Diffusion1D_Gaussian(x;
    params = (T0 = 100., K = 1e-6, σ = 0.1 ) )
    u      = Diffusion1D_Gaussian_u_fwd(x, params)
    f_cl   = x -> Diffusion1D_Gaussian_u_fwd(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = gradu[2] - params.K*(hessu[1,1])
    return (u=u, ∇u=gradu, s=s)
end