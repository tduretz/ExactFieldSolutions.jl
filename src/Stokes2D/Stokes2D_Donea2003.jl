function Stokes2D_Donea2003_p(x, params)
    return x[1]*(1-x[1])
end

function Stokes2D_Donea2003_V(x, params)
    vx      =  x[1]^2*(1 - x[1])^2*(4*x[2]^3 - 6*x[2]^2 + 2*x[2])
    vy      = -x[2]^2*(1 - x[2])^2*(4*x[1]^3 - 6*x[1]^2 + 2*x[1])  
    return [vx; vy] 
end

function TotalStress( x,;
        params = (η=1,) )
    p = Stokes2D_Donea2003_p(x, params)
    v = Stokes2D_Donea2003_V(x, params)
    f_cl = x -> Stokes2D_Donea2003_V(x, params)
    L = ForwardDiff.jacobian(f_cl, x)
    # Deviatoric strain rate
    ε̇ = 1/2*(L + L')
    # Deviatoric stress
    τ = 2*params.η*ε̇
    # Total stress
    return τ - [1 0; 0 1].*p
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
function Stokes2D_Donea2003(X;
    params = (η=1,) )
    p = Stokes2D_Donea2003_p(X, params)
    v = Stokes2D_Donea2003_V(X, params)
    f_cl = X -> Stokes2D_Donea2003_V(X, params)
    L = ForwardDiff.jacobian(f_cl, X)
    # Deviatoric strain rate
    ε̇ = 1/2*(L + L')
    # Deviatoric stress
    τ = 2*params.η*ε̇
    # Total stress
    σ = TotalStress( X )
    g_cl = X -> TotalStress(X; params)
    ∂σ∂x = ForwardDiff.jacobian(g_cl, X)
    # Stokes residual is the source for this solution
    s = -[∂σ∂x[1] + ∂σ∂x[3,2]; ∂σ∂x[2,1] + ∂σ∂x[4,2]]
    # # Just to check
    # x, y = X[1], X[2]
    # sx  = -p^2*(24*y - 12) - 4*x^2*(4*y^3 - 6*y^2 + 2*y) - 8*x*(2*x - 2)*(4*y^3 - 6*y^2 + 2*y) - 2*x + 1.0*y^2*(2*y - 2)*(12*x^2 - 12*x + 2) + 2.0*y*(1 - y)^2*(12*x^2 - 12*x + 2) - 4*(1 - x)^2*(4*y^3 - 6*y^2 + 2*y) + 1
    # sy  = -2*p*(1 - x)*(12*y^2 - 12*y + 2) - x^2*(2*x - 2)*(12*y^2 - 12*y + 2) + 1.0*y^2*(1 - y)^2*(24*x - 12) + 4*y^2*(4*x^3 - 6*x^2 + 2*x) + 8*y*(2*y - 2)*(4*x^3 - 6*x^2 + 2*x) + 4*(1 - y)^2*(4*x^3 - 6*x^2 + 2*x)
    # @show sy-s[2]
    return (p=p, V=v, L=L, ε̇=ε̇, τ=τ, s=s)
end
