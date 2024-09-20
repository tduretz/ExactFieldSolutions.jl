
function Poisson3D_Sevilla2018_u(x, params)
    @unpack α, β, a, b, c, d, e, f = params
    return exp(α*sin(a*x[1] + c*x[2] + e*x[3]) + β*cos(b*x[1] + d*x[2] + f*x[3]))
end

@doc raw"""
    sol = Poisson3D_Sevilla2018(x; params)  

Evaluates the manufactured solution of [Sevilla et al. (2018)](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.5833):

    x      : is the coordinate vector or a tuple
    params : optional parameter array, default = (α = 0.1, β = 0.3, a = 5.1, b = 4.3, c = -6.2, d = 3.4, e = 1.8, f = 1.7)
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Poisson3D_Sevilla2018( [0; 0; 0])
(u = 1.3498588075760032, ∇u = [0.6884279918637617, -0.8369124606971221, 0.2429745853636806], s = 12.425585309617865)
```
```julia-repl
julia> Poisson3D_Sevilla2018( (0, 0, 0) )
(u = 1.3498588075760032, ∇u = (x = 0.6884279918637617, y = -0.8369124606971221, z = 0.2429745853636806), s = 12.425585309617865)
```
"""
function Poisson3D_Sevilla2018(x;
    params = (α = 0.1, β = 0.3, a = 5.1, b = 4.3, c = -6.2, d = 3.4, e = 1.8, f = 1.7) )
    u      = Poisson3D_Sevilla2018_u(x, params)
    f_cl   = x -> Poisson3D_Sevilla2018_u(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = -(hessu[1,1] + hessu[2,2] + hessu[3,3])
    return (u=u, ∇u=gradu, s=s)
end

function Poisson3D_Sevilla2018(coords::Union{Tuple, NamedTuple};
    params = (α = 0.1, β = 0.3, a = 5.1, b = 4.3, c = -6.2, d = 3.4, e = 1.8, f = 1.7) )
    X = SVector(values(coords)...)
    sol = Poisson3D_Sevilla2018(X)
    return (u=sol.u, ∇u =(x=sol.∇u[1], y=sol.∇u[2], z=sol.∇u[3]), s=sol.s )
end