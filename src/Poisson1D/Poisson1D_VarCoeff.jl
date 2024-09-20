function Poisson1D_VarCoeff_u_fwd(x)
    u = cos(x[1])
    return u
end

function Poisson1D_VarCoeff_coeff(x)
    β    = 1.0 
    # r    = 0.5
    # σ    = 0.1
    # a    = sqrt(x[1]^2) - r
    # logβ = 2*erf(-(a/σ))  
    # β    = 10^logβ
    return β
end

function Poisson1D_Flux(X)
    β      = Poisson1D_VarCoeff_coeff(X)
    f_cl   = X -> Poisson1D_VarCoeff_u_fwd(X)
    ∇u     = ForwardDiff.gradient(f_cl, X)
    return -β.*∇u
end

@doc raw"""
    sol = Poisson1D_VarCoeff(x)  

Evaluates the manufactured solution of [Helgadottir et al. (2018)](https://onlinelibrary.wiley.com/doi/10.1155/2018/9216703):

    x      : is the coordinate vector or a tuple
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Poisson1D_VarCoeff( [0, 0.8] )
(u = 0.7173560908995228, ∇u = [0.0, 0.6967067093471654], s = 1.4347121817990456)
```
```julia-repl
julia> Poisson1D_VarCoeff( (0, 0.8) )
(u = 1.0, q = (x = 0.0, y = 0.0), s = -1.0, β = 1.0)
```
"""
function Poisson1D_VarCoeff(x)
    u      = Poisson1D_VarCoeff_u_fwd(x)
    ∂u∂x   = ForwardDiff.gradient( x->Poisson1D_VarCoeff_u_fwd(x), x)
    q_cl   = x -> Poisson1D_Flux(x)
    q      = Poisson1D_Flux(x)
    ∇q     = ForwardDiff.jacobian(q_cl, x)
    s      = -∇q[1,1]
    return (u=u, ∇u=∂u∂x, q=q, s=s, β=Poisson1D_VarCoeff_coeff(x))
end

function Poisson1D_VarCoeff(coords::Tuple)
    X = SVector(values(coords)...)
    sol = Poisson1D_VarCoeff(X)
    return (u=sol.u, q =(x=sol.q[1], y=sol.q[2]), s=sol.s, β=sol.β )
end