@doc raw"""
    sol = Elasticity2D_Hole(x; params)  

Evaluates the solution for an elastic medium with a circular hole under far field stress (derived from Kirsch solution, Jaeger & Cook etc...):

    x      : is the coordinate vector or tuple
    params : optional parameter array, default = (r0 = 0.2, P_inf = 1e-2, P0 = -1e-2, tau_inf = -1e-2)
and returns:

    sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# Examples
```julia-repl
julia> Elasticity2D_Hole([0,0])
(p = -0.01, τ = [0.0 -0.0; -0.0 -0.0])
```
```julia-repl
julia> Elasticity2D_Hole((0,0))
(p = -0.01, τ = (xx = 0.0, xy = -0.0, yx = -0.0, yy = -0.0))
```
"""
function Elasticity2D_Hole(x;
    params = (r0 = 0.2, P_inf = 1e-2, P0 = -1e-2, tau_inf = -1e-2))
    @unpack r0, P_inf, P0, tau_inf = params
    rad      = sqrt(x[1]^2 + x[2]^2)
    irad     = 1.0/max(rad,r0/2); 
    irad>1/r0 ? irad = 0 : nothing
    r0_r     = r0*irad
    cost     = x[1]*irad
    sint     = x[2]*irad
    C1       = -3*r0_r^4 + 2*r0_r^2
    C2       = -3*r0_r^4 + 2*r0_r^2
    cos2t    = cost^2 - sint^2
    sin2t    = 2*cost*sint
    Pt       = P_inf             - 2*tau_inf* r0_r^2.0.*cos2t
    tau_rr   = (P_inf-P0)*r0_r^2 + tau_inf*(C1 - 1)*cos2t
    tau_rt   =                     tau_inf*(C2 + 1)*sin2t
    tau_tt   = -tau_rr
    Txx      =  cos2t*tau_rr - sin2t*tau_rt
    Txy      = -sin2t*tau_tt + tau_rt*cos2t
    Tyy      = -Txx
    rad.<r0 ? Pt = P0 : nothing
    return (p= Pt, τ=[Txx Txy; Txy Tyy])
end

function Elasticity2D_Hole(coords::Tuple;
    params = (r0 = 0.2, P_inf = 1e-2, P0 = -1e-2, tau_inf = -1e-2) )
    X = SVector(values(coords)...)
    sol = Elasticity2D_Hole(X; params)
    return (p=sol.p, τ=(xx=sol.τ[1,1], xy=sol.τ[1,2], yx=sol.τ[2,1], yy=sol.τ[2,2]) )
end