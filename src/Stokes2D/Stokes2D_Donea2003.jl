function Stokes2D_Donea2003_p(x, params)
    return x[1]*(1-x[1])
end

function Stokes2D_Donea2003_V(x, params)
    vx      =  x[1]^2*(1 - x[1])^2*(4*x[2]^3 - 6*x[2]^2 + 2*x[2])
    vy      = -x[2]^2*(1 - x[2])^2*(4*x[1]^3 - 6*x[1]^2 + 2*x[1])  
    return [vx; vy] 
end
 
@doc raw"""
    sol = Stokes2D_Donea2003(x; params)  

Evaluates a manufactured solution from [Donea & Huerta (2003)](https://onlinelibrary.wiley.com/doi/book/10.1002/0470013826):

    x      : is the coordinate vector 
    params : optional parameter array
and returns:

    sol    : tuple containing the solution fields p (pressure), V (velocity vector), L (velocity gratdient tensor), ε̇ (deviatoric strain rate tensor) and τ (deviatoric stress tensor)

# Examples
```julia-repl
julia> Stokes2D_Donea2003([1/4;1/4])
(p = 0.1875, V = [0.006591796875, -0.006591796875], L = [0.03515625 -0.0087890625; 0.0087890625 -0.03515625], ε̇ = [0.03515625 0.0; 0.0 -0.03515625], τ = [0.0703125 0.0; 0.0 -0.0703125])
```
"""
function Stokes2D_Donea2003(x;
    params = (η=1,) )
    p = Stokes2D_Donea2003_p(x, params)
    v = Stokes2D_Donea2003_V(x, params)
    f_cl = x -> Stokes2D_Donea2003_V(x, params)
    L = ForwardDiff.jacobian(f_cl, x)
    # Postprocess deviatoric stress and strain rate
    ε̇ = 1/2*(L + L')
    # Check position
    τ = 2*params.η*ε̇
    return (p=p, V=v, L=L, ε̇=ε̇, τ=τ)
end
