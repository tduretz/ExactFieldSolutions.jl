
function Diffusion2D_Gaussian_u_fwd(x, params)
    @unpack T0, K, σ  = params
    return  T0/(1 + 2*π*K*x[3]/σ^2)*exp(-(x[1]^2 + x[2]^2) /( 2*σ^2/π + 4*K*x[3]))
end

@doc raw"""
    sol = Diffusion2D_Gaussian(x; params)  

Evaluates the manufactured solution of an initial 2D Gaussian anomaly:

    x      : is the coordinate and time vector, x[1:2] is space and x[3] is time
    params : optional parameter array
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Diffusion2D_Gaussian([0 0 0])
(u = 100.0, ∇u = [0.0 0.0 -0.06283185307179585], s = 0.0)
```
"""
function Diffusion2D_Gaussian(x;
    params = (T0 = 100., K = 1e-6, σ = 0.1 ) )
    u      = Diffusion2D_Gaussian_u_fwd(x, params)
    f_cl   = x -> Diffusion2D_Gaussian_u_fwd(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = gradu[3] - params.K*(hessu[1,1] + hessu[2,2])
    return (u=u, ∇u=gradu, s=s)
end