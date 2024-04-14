function Stokes2D_Schmid2003_p(x, params)
    @unpack mm, mc, rc, gr, er  = params
    # Check position
    (x[1]^2 + x[2]^2) < rc^2 ? in=true : in=false
    # Outside
    pmat   = (4*mm*rc^2*(mm - mc)*(er*x[1]^2 + gr*x[1]*x[2] - er*x[2]^2))/((x[1]^2 + x[2]^2)^2*(mm + mc))
    # Inside
    pinc   = 0.0
    # Where it's needed
    p      = (in==0)*pmat + (in==1)*pinc 
    return p
end

function Stokes2D_Schmid2003_V(x, params)
    @unpack mm, mc, rc, gr, er  = params
    # Check position
    (x[1]^2 + x[2]^2) < rc^2 ? in=true : in=false
    # Outside
    A = mm*(mc-mm)/(mc+mm);
    v1mat = (2*A*er*rc^4*x[1]^3 + 3*A*gr*rc^4*x[1]^2*x[2] - 6*A*er*rc^4*x[1]*x[2]^2 - A*gr*rc^4*x[2]^3 - 4*A*er*rc^2*x[1]^5 - 
        4*A*gr*rc^2*x[1]^4*x[2] - 4*A*gr*rc^2*x[1]^2*x[2]^3 + 4*A*er*rc^2*x[1]*x[2]^4 + 2*er*mm*x[1]^7 + 2*gr*mm*x[1]^6*x[2] + 
        6*er*mm*x[1]^5*x[2]^2 + 6*gr*mm*x[1]^4*x[2]^3 + 6*er*mm*x[1]^3*x[2]^4 + 6*gr*mm*x[1]^2*x[2]^5 + 2*er*mm*x[1]*x[2]^6 + 
        2*gr*mm*x[2]^7)/(2*mm*(x[1]^2 + x[2]^2)^3)
    v2mat = -(A*gr*rc^4*x[1]^3 - 6*A*er*rc^4*x[1]^2*x[2] - 3*A*gr*rc^4*x[1]*x[2]^2 + 2*A*er*rc^4*x[2]^3 + 4*A*er*rc^2*x[1]^4*x[2] + 
        4*A*gr*rc^2*x[1]^3*x[2]^2 + 4*A*gr*rc^2*x[1]*x[2]^4 - 4*A*er*rc^2*x[2]^5 + 2*er*mm*x[1]^6*x[2] + 6*er*mm*x[1]^4*x[2]^3 + 
        6*er*mm*x[1]^2*x[2]^5 + 2*er*mm*x[2]^7)/(2*mm*(x[1]^2 + x[2]^2)^3)
    # Inside
    v1inc  = (4*er*mm*x[1] + 3*gr*mm*x[2] + gr*mc*x[2])/(2*(mm + mc));
    v2inc  = -(4*er*mm*x[2] - gr*mm*x[1] + gr*mc*x[1])/(2*(mm + mc));
    # Where it's needed
    vx      = (in==0)*v1mat + (in==1)*v1inc  
    vy      = (in==0)*v2mat + (in==1)*v2inc  
    return [vx; vy] 
end
 
function Stokes2D_Schmid2003_L(x, params) # just to check if automatic derivatives are correct
    @unpack mm, mc, rc, gr, er  = params
    # Check position
    (x[1]^2 + x[2]^2) < rc^2 ? in=true : in=false
    # Outside
    A = mm*(mc-mm)/(mc+mm);
    dvxdxmat = (- 3*A*er*rc^4*x[1]^4 - 6*A*gr*rc^4*x[1]^3*x[2] + 18*A*er*rc^4*x[1]^2*x[2]^2 + 6*A*gr*rc^4*x[1]*x[2]^3 - 3*A*er*rc^4*x[2]^4 + 2*A*er*rc^2*x[1]^6 + 4*A*gr*rc^2*x[1]^5*x[2] - 10*A*er*rc^2*x[1]^4*x[2]^2 - 10*A*er*rc^2*x[1]^2*x[2]^4 - 4*A*gr*rc^2*x[1]*x[2]^5 + 2*A*er*rc^2*x[2]^6 + er*mm*x[1]^8 + 4*er*mm*x[1]^6*x[2]^2 + 6*er*mm*x[1]^4*x[2]^4 + 4*er*mm*x[1]^2*x[2]^6 + er*mm*x[2]^8)/(mm*(x[1]^2 + x[2]^2)^4);
    dvxdymat = (3*A*gr*rc^4*x[1]^4 - 24*A*er*rc^4*x[1]^3*x[2] - 18*A*gr*rc^4*x[1]^2*x[2]^2 + 24*A*er*rc^4*x[1]*x[2]^3 + 3*A*gr*rc^4*x[2]^4 - 4*A*gr*rc^2*x[1]^6 + 24*A*er*rc^2*x[1]^5*x[2] + 8*A*gr*rc^2*x[1]^4*x[2]^2 + 16*A*er*rc^2*x[1]^3*x[2]^3 + 12*A*gr*rc^2*x[1]^2*x[2]^4 - 8*A*er*rc^2*x[1]*x[2]^5 + 2*gr*mm*x[1]^8 + 8*gr*mm*x[1]^6*x[2]^2 + 12*gr*mm*x[1]^4*x[2]^4 + 8*gr*mm*x[1]^2*x[2]^6 + 2*gr*mm*x[2]^8)/(2*mm*(x[1]^2 + x[2]^2)^4);
    dvydxmat = (A*rc^2*(3*gr*rc^2*x[1]^4 - 24*er*rc^2*x[1]^3*x[2] - 18*gr*rc^2*x[1]^2*x[2]^2 + 24*er*rc^2*x[1]*x[2]^3 + 3*gr*rc^2*x[2]^4 + 8*er*x[1]^5*x[2] + 12*gr*x[1]^4*x[2]^2 - 16*er*x[1]^3*x[2]^3 + 8*gr*x[1]^2*x[2]^4 - 24*er*x[1]*x[2]^5 - 4*gr*x[2]^6))/(2*mm*(x[1]^2 + x[2]^2)^4);
    dvydymat = -(- 3*A*er*rc^4*x[1]^4 - 6*A*gr*rc^4*x[1]^3*x[2] + 18*A*er*rc^4*x[1]^2*x[2]^2 + 6*A*gr*rc^4*x[1]*x[2]^3 - 3*A*er*rc^4*x[2]^4 + 2*A*er*rc^2*x[1]^6 + 4*A*gr*rc^2*x[1]^5*x[2] - 10*A*er*rc^2*x[1]^4*x[2]^2 - 10*A*er*rc^2*x[1]^2*x[2]^4 - 4*A*gr*rc^2*x[1]*x[2]^5 + 2*A*er*rc^2*x[2]^6 + er*mm*x[1]^8 + 4*er*mm*x[1]^6*x[2]^2 + 6*er*mm*x[1]^4*x[2]^4 + 4*er*mm*x[1]^2*x[2]^6 + er*mm*x[2]^8)/(mm*(x[1]^2 + x[2]^2)^4);
    # Inside
    dvxdxinc = (2*er*mm)/(mm + mc);
    dvxdyinc = (gr*(3*mm + mc))/(2*(mm + mc));
    dvydxinc = (gr*(mm - mc))/(2*(mm + mc));
    dvydyinc = -(2*er*mm)/(mm + mc);
    # Where it's needed
    dvxdx      = (in==0)*dvxdxmat + (in==1)*dvxdxinc  
    dvxdy      = (in==0)*dvxdymat + (in==1)*dvxdyinc 
    dvydx      = (in==0)*dvydxmat + (in==1)*dvydxinc  
    dvydy      = (in==0)*dvydymat + (in==1)*dvydyinc 
    return [dvxdx dvxdy; dvydx dvydy]
end

@doc raw"""
    sol = Stokes2D_Schmid2003(x; params)  

Evaluates the manufactured solution of [Schmid & Podladchikov (2003)](https://academic.oup.com/gji/article/155/1/269/713923):

    x      : is the coordinate vector 
    params : optional parameter array
and returns:

    sol    : tuple containing the solution fields p (pressure), V (velocity vector), L (velocity gratdient tensor), ε̇ (deviatoric strain rate tensor) and τ (deviatoric stress tensor)

# Examples
```julia-repl
julia> Stokes2D_Schmid2003( [0, 0] )
(p = 0.0, V = [0.0, 0.0], L = [-0.019801980198019802 0.0; 0.0 0.019801980198019802], ε̇ = [-0.019801980198019802 0.0; 0.0 0.019801980198019802], τ = [-3.9603960396039604 0.0; 0.0 3.9603960396039604])
```
"""
function Stokes2D_Schmid2003(x;
    params = (mm = 1.0, mc = 100, rc = 0.2, gr = 0.0, er =-1.0) )
    p = Stokes2D_Schmid2003_p(x, params)
    v = Stokes2D_Schmid2003_V(x, params)
    # L = Stokes2D_Schmid2003_L(x, params) # just to check if automatic derivatives are correct
    f_cl = x -> Stokes2D_Schmid2003_V(x, params)
    L = ForwardDiff.jacobian(f_cl, x)
    # Postprocess deviatoric stress and strain rate
    ε̇ = 1/2*(L + L')
    # Check position
    (x[1]^2 + x[2]^2) < params.rc^2 ? η=params.mc : η=params.mm
    τ = 2*η*ε̇
    return (p=p, V=v, L=L, ε̇=ε̇, τ=τ, η=η)
end
