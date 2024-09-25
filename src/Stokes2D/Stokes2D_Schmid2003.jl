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
    xi, yi = x[1], x[2]
  
    # Check position
    isinside = (x[1]^2 + x[2]^2) < rc^2
    # Outside
    A = mm*(mc-mm)/(mc+mm);

    v1mat = @fastpow (2*A*er*rc^4*xi^3 + 3*A*gr*rc^4*xi^2*yi - 6*A*er*rc^4*xi*yi^2 - A*gr*rc^4*yi^3 - 4*A*er*rc^2*xi^5 - 
        4*A*gr*rc^2*xi^4*yi - 4*A*gr*rc^2*xi^2*yi^3 + 4*A*er*rc^2*xi*yi^4 + 2*er*mm*xi^7 + 2*gr*mm*xi^6*yi + 
        6*er*mm*xi^5*yi^2 + 6*gr*mm*xi^4*yi^3 + 6*er*mm*xi^3*yi^4 + 6*gr*mm*xi^2*yi^5 + 2*er*mm*xi*yi^6 + 
        2*gr*mm*yi^7)/(2*mm*(xi^2 + yi^2)^3)
    v2mat = @fastpow -(A*gr*rc^4*xi^3 - 6*A*er*rc^4*xi^2*yi - 3*A*gr*rc^4*xi*yi^2 + 2*A*er*rc^4*yi^3 + 4*A*er*rc^2*xi^4*yi + 
        4*A*gr*rc^2*xi^3*yi^2 + 4*A*gr*rc^2*xi*yi^4 - 4*A*er*rc^2*yi^5 + 2*er*mm*xi^6*yi + 6*er*mm*xi^4*yi^3 + 
        6*er*mm*xi^2*yi^5 + 2*er*mm*yi^7)/(2*mm*(xi^2 + yi^2)^3)
    # Inside
    v1inc  = (4*er*mm*xi + 3*gr*mm*yi + gr*mc*yi)/(2*(mm + mc));
    v2inc  = -(4*er*mm*yi - gr*mm*xi + gr*mc*xi)/(2*(mm + mc));
    # Where it's needed
    vx      = (isinside==0)*v1mat + (isinside==1)*v1inc  
    vy      = (isinside==0)*v2mat + (isinside==1)*v2inc  
    return @SVector [vx; vy] 
end
 
function Stokes2D_Schmid2003_L(x, params) # just to check if automatic derivatives are correct
    @unpack mm, mc, rc, gr, er  = params
    xi, yi = x[1], yi

    # Check position
    isinside = (xi^2 + yi^2) < rc^2 
    # Outside
    A = mm*(mc-mm)/(mc+mm);
    dvxdxmat = @fastpow (- 3*A*er*rc^4*xi^4 - 6*A*gr*rc^4*xi^3*yi + 18*A*er*rc^4*xi^2*yi^2 + 6*A*gr*rc^4*xi*yi^3 - 3*A*er*rc^4*yi^4 + 2*A*er*rc^2*xi^6 + 4*A*gr*rc^2*xi^5*yi - 10*A*er*rc^2*xi^4*yi^2 - 10*A*er*rc^2*xi^2*yi^4 - 4*A*gr*rc^2*xi*yi^5 + 2*A*er*rc^2*yi^6 + er*mm*xi^8 + 4*er*mm*xi^6*yi^2 + 6*er*mm*xi^4*yi^4 + 4*er*mm*xi^2*yi^6 + er*mm*yi^8)/(mm*(xi^2 + yi^2)^4);
    dvxdymat = @fastpow (3*A*gr*rc^4*xi^4 - 24*A*er*rc^4*xi^3*yi - 18*A*gr*rc^4*xi^2*yi^2 + 24*A*er*rc^4*xi*yi^3 + 3*A*gr*rc^4*yi^4 - 4*A*gr*rc^2*xi^6 + 24*A*er*rc^2*xi^5*yi + 8*A*gr*rc^2*xi^4*yi^2 + 16*A*er*rc^2*xi^3*yi^3 + 12*A*gr*rc^2*xi^2*yi^4 - 8*A*er*rc^2*xi*yi^5 + 2*gr*mm*xi^8 + 8*gr*mm*xi^6*yi^2 + 12*gr*mm*xi^4*yi^4 + 8*gr*mm*xi^2*yi^6 + 2*gr*mm*yi^8)/(2*mm*(xi^2 + yi^2)^4);
    dvydxmat = @fastpow (A*rc^2*(3*gr*rc^2*xi^4 - 24*er*rc^2*xi^3*yi - 18*gr*rc^2*xi^2*yi^2 + 24*er*rc^2*xi*yi^3 + 3*gr*rc^2*yi^4 + 8*er*xi^5*yi + 12*gr*xi^4*yi^2 - 16*er*xi^3*yi^3 + 8*gr*xi^2*yi^4 - 24*er*xi*yi^5 - 4*gr*yi^6))/(2*mm*(xi^2 + yi^2)^4);
    dvydymat = @fastpow -(- 3*A*er*rc^4*xi^4 - 6*A*gr*rc^4*xi^3*yi + 18*A*er*rc^4*xi^2*yi^2 + 6*A*gr*rc^4*xi*yi^3 - 3*A*er*rc^4*yi^4 + 2*A*er*rc^2*xi^6 + 4*A*gr*rc^2*xi^5*yi - 10*A*er*rc^2*xi^4*yi^2 - 10*A*er*rc^2*xi^2*yi^4 - 4*A*gr*rc^2*xi*yi^5 + 2*A*er*rc^2*yi^6 + er*mm*xi^8 + 4*er*mm*xi^6*yi^2 + 6*er*mm*xi^4*yi^4 + 4*er*mm*xi^2*yi^6 + er*mm*yi^8)/(mm*(xi^2 + yi^2)^4);
    # Inside
    dvxdxinc = (2*er*mm)/(mm + mc);
    dvxdyinc = (gr*(3*mm + mc))/(2*(mm + mc));
    dvydxinc = (gr*(mm - mc))/(2*(mm + mc));
    dvydyinc = -(2*er*mm)/(mm + mc);
    # Where it's needed
    dvxdx      = (isinside==0)*dvxdxmat + (isinside==1)*dvxdxinc  
    dvxdy      = (isinside==0)*dvxdymat + (isinside==1)*dvxdyinc 
    dvydx      = (isinside==0)*dvydxmat + (isinside==1)*dvydxinc  
    dvydy      = (isinside==0)*dvydymat + (isinside==1)*dvydyinc 
    return @SMatrix [dvxdx dvxdy; dvydx dvydy]
end

@doc raw"""
    sol = Stokes2D_Schmid2003(x; params)  

Evaluates the analytical solution of [Schmid & Podladchikov (2003)](https://academic.oup.com/gji/article/155/1/269/713923):

    x      : is the coordinate vector or tuple
    params : optional parameter array, default (mm = 1.0, mc = 100, rc = 0.2, gr = 0.0, er =-1.0)
and returns:

    sol    : tuple containing the solution fields p (pressure), V (velocity vector), L (velocity gratdient tensor), ε̇ (deviatoric strain rate tensor) and τ (deviatoric stress tensor)

# Examples
```julia-repl
julia> Stokes2D_Schmid2003( [0, 0] )
(p = 0.0, V = [0.0, 0.0], L = [-0.019801980198019802 0.0; 0.0 0.019801980198019802], ε̇ = [-0.019801980198019802 0.0; 0.0 0.019801980198019802], τ = [-3.9603960396039604 0.0; 0.0 3.9603960396039604])
```
```julia-repl
julia> Stokes2D_Schmid2003( (0, 0) )
(p = 0.0, V = (x = 0.0, y = 0.0), L = (xx = -0.019801980198019802, xy = 0.0, yx = 0.0, yy = 0.019801980198019802), ε̇ = (xx = -0.019801980198019802, xy = 0.0, yx = 0.0, yy = 0.019801980198019802), τ = (xx = -3.9603960396039604, xy = 0.0, yx = 0.0, yy = 3.9603960396039604), η = 100)
```
"""
function Stokes2D_Schmid2003(x;
    params = (mm = 1.0, mc = 100.0, rc = 0.2, gr = 0.0, er =-1.0) )
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

function Stokes2D_Schmid2003(coords::Union{Tuple, NamedTuple};
    params = (mm = 1.0, mc = 100.0, rc = 0.2, gr = 0.0, er =-1.0) )
    X = SVector(values(coords)...)
    sol = Stokes2D_Schmid2003(X; params)
    return (p=sol.p, 
    V=(x=sol.V[1], y=sol.V[2]),
    L=(xx=sol.L[1,1], xy=sol.L[1,2], yx=sol.L[2,1], yy=sol.L[2,2]), 
    ε̇=(xx=sol.ε̇[1,1], xy=sol.ε̇[1,2], yx=sol.ε̇[2,1], yy=sol.ε̇[2,2]), 
    τ=(xx=sol.τ[1,1], xy=sol.τ[1,2], yx=sol.τ[2,1], yy=sol.τ[2,2]),
    η=sol.η) 
end