function Wave1D_dAlembert_u_fwd(x, params)
    @unpack ρ, k, G  = params
    return real(exp(im*k*(x[1] - sqrt(G/ρ)*x[2])))
end

@doc raw"""
    sol = Wave1D_dAlembert(x; params)  

Evaluates the manufactured solution of an initial 1D wave:

    x      : is the coordinate and time vector, x[1] is space and x[2] is time 
    params : optional parameter tuple, default: (ρ=2.0, k=8.0, G=1.0) 
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Wave1D_dAlembert([0 0])
(u = 1.0, ∇u = [0.0 -0.0], s = 0.0, G = 1.0, ρ = 2.0)
```

```julia-repl
julia> params = (ρ = 200.0, k = 80.0, G = 0.5)
(ρ = 200.0, k = 80.0, G = 0.5)

julia> Wave1D_dAlembert([0 2]; params)
(u = -0.14550003380861354, ∇u = [79.14865972987054 -3.957432986493527], s = -4.440892098500626e-16, G = 0.5, ρ = 200.0)
```
"""
function Wave1D_dAlembert(x;
    params = (ρ=2.0, k=8.0, G=1.0) )
    c      = sqrt(params.G/params.ρ)
    u      = Wave1D_dAlembert_u_fwd(x, params)
    f_cl   = x -> Wave1D_dAlembert_u_fwd(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = hessu[2,2] - c^2*(hessu[1,1])
    return (u=u, ∇u=gradu, s=s, G=params.G, ρ=params.ρ)
end