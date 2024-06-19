function Wave1D_dAlembert_u_fwd(x, params)
    @unpack c, k  = params
    return real(exp(im*k*(x[1] - c*x[2])))
end

@doc raw"""
    sol = Diffusion1D_Gaussian(x; params)  

Evaluates the manufactured solution of an initial 1D wave:

    x      : is the coordinate and time vector, x[1] is space and x[2] is time 
    params : optional parameter array, default: (c = 1.0, k = 8.0) 
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Wave1D_dAlembert([0 0])
(u = 1.0, ∇u = [0.0 -0.0], s = 0.0)
```

```julia-repl
julia> Wave1D_dAlembert([0 0]; params)
(u = 1.0, ∇u = [0.0 0.0], s = 0.0)

julia> Wave1D_dAlembert([0 2]; params)
(u = 0.424179007336997, ∇u = [10.866940344079486 10.866940344079486], s = 0.0)
```
"""
function Wave1D_dAlembert(x;
    params = (c=1.0, k=8.0) )
    u      = Wave1D_dAlembert_u_fwd(x, params)
    f_cl   = x -> Wave1D_dAlembert_u_fwd(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = hessu[2,2] - params.c^2*(hessu[1,1])
    return (u=u, ∇u=gradu, s=s, G=params.c)
end