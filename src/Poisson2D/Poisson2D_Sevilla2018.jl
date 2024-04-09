
function Poisson2D_Sevilla2018_u_fwd(x, params)
    @unpack α, β, a, b, c, d = params
    return exp(α*sin(a*x[1] + c*x[2]) + β*cos(b*x[1] + d*x[2]))
end

@doc raw"""
    sol = Poisson2D_Sevilla2018(x; params)  

Evaluates the manufactured solution of [Sevilla et al. (2018)](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.5833):

    x      : is the coordinate vector 
    params : optional parameter array
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Poisson2D_Sevilla2018( [0, 0] )
(u = 1.3498588075760032, ∇u = [0.6884279918637617, -0.8369124606971221], s = 11.298993148814933)
```
"""
function Poisson2D_Sevilla2018(x;
    params = (α = 0.1, β = 0.3, a = 5.1, b = 4.3, c = -6.2, d = 3.4) )
    u      = Poisson2D_Sevilla2018_u_fwd(x, params)
    f_cl   = x -> Poisson2D_Sevilla2018_u_fwd(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
     s1      = -(hessu[1,1] + hessu[2,2])
    T = u
    @unpack α, β, a, b, c, d = params
    s = T*(-a*α*cos(a*x[1] + c*x[2]) + b*β*sin(b*x[1] + d*x[2]))*(a*α*cos(a*x[1] + c*x[2]) - b*β*sin(b*x[1] + d*x[2])) + T*(a^2*α*sin(a*x[1] + c*x[2]) + b^2*β*cos(b*x[1] + d*x[2])) + T*(-α*c*cos(a*x[1] + c*x[2]) + β*d*sin(b*x[1] + d*x[2]))*(α*c*cos(a*x[1] + c*x[2]) - β*d*sin(b*x[1] + d*x[2])) + T*(α*c^2*sin(a*x[1] + c*x[2]) + β*d^2*cos(b*x[1] + d*x[2]))
    @show s, s1
    return (u=u, ∇u=gradu, s=s)
end

# I've tried to use Enzyme BUT:
# 1) outputs are difficult gradient components end nested in complicated data structures (arrays of arrays or tuples of arrays)
# 2) never managed to compute higher order derivatives
# 3) overall package compilation seemed much slower than with ForwardDiff

# Poisson2D_Sevilla2018_gradu(x, y, params)   = autodiff(Reverse, Poisson2D_Sevilla2018_u, Active(x), Active(y), Const(params))
# Poisson2D_Sevilla2018_dudx( x, y, params)   = autodiff(Reverse, Poisson2D_Sevilla2018_u, Active(x), Const(y), Const(params))
# Poisson2D_Sevilla2018_d2udx2( x, y, params) = autodiff_deferred(Reverse, Poisson2D_Sevilla2018_dudx, Active(x), Const(y), Const(params))
# Poisson2D_Sevilla2018_dudx(x, y, params)   = autodiff(Reverse, Poisson2D_Sevilla2018_u, Active(x), Const(y), Const(params))
# Poisson2D_Sevilla2018_d2udx2(x, y, params) = autodiff_deferred(Reverse, Poisson2D_Sevilla2018_dudx, Active(x), Const(y), Const(params))
# Poisson2D_Sevilla2018_d2udy2(x, y, params) = autodiff_deferred(Reverse, Poisson2D_Sevilla2018_dudx, Const(x), Active(y), Const(params))
# function Poisson2D_Sevilla2018_s(x, y, params)
# end

# function Poisson2D_Sevilla2018_u(x, y, params)
#     @unpack α, β, a, b, c, d = params
#     return exp(α*sin(a*x + c*y) + β*cos(b*x + d*y))
# end

# function Poisson2D_Sevilla2018_enz(x, y;
#     params = (α = 0.1, β = 0.3, a = 5.1, b = 4.3, c = -6.2, d = 3.4) )
#     u      = Poisson2D_Sevilla2018_u(x, y, params)
#     gradu  = Poisson2D_Sevilla2018_gradu(x, y, params)
#     # ∂2u∂x2 = Poisson2D_Sevilla2018_d2udx2(x, y, params)
#     # ∂2u∂y2 = Poisson2D_Sevilla2018_d2udy2(x, y, params)
#     # Poisson2D_Sevilla2018_d2udx2(x, y, params)
#     s      = 1
#     # @show ∂2u∂y2[1][1]
#     return (u=u, ∂u∂x=gradu[1][1], ∂u∂y=gradu[1][2], s=s)
# end