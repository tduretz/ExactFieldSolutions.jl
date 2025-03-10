using UnPack, StaticArrays
@doc raw"""
    sol = Poroelasstic2D_Rozhko2008(x; params)  

Evaluates the solution for an elastic medium with a circular hole under far field stress (derived from Kirsch solution, Jaeger & Cook etc...):

    x      : is the coordinate vector or tuple
    params : optional parameter array, default = (r_in=1.0, r_out=20.0, P0=0.0, dPf=1.0, m=0.0, nu= 0.49, G=1.0)
and returns:

    sol    : tuple containing the solution fields u (displacement in Cartesian coordinates), u_pol (displacement in polar coordinates), pt, pf, σ, σ_pol

# Examples
```julia-repl
julia> Poroelasstic2D_Rozhko2008([1,1])
(u = [0.005436794860536682 0.005436794860536681], u_pol = [0.007688789027611311 -0.0], pt = 0.014066799864165765, pf = 0.8843108934201204, σ = [-0.019291833584202365 0.003918775290027451; 0.003918775290027451 -0.008841766144129164], σ_pol = [-0.007535507714120015 -0.0; -0.0 -0.020598092014211516])
```
```julia-repl
julia> Poroelasstic2D_Rozhko2008((1,1))
(pt = 0.014066799864165765, pf = 0.8843108934201204, u = (x = 0.005436794860536682, y = 0.005436794860536681), σ = (xx = -0.019291833584202365, xy = 0.003918775290027451, yx = 0.003918775290027451, yy = -0.008841766144129164))
```
"""
function Poroelasticity2D_Rozhko2008(X;  
    params = (r_in=1.0, r_out=20.0, P0=0.0, dPf=1.0, m=0.0, nu= 0.49, G=1.0) )
    @unpack r_in, r_out, P0, dPf, m, nu, G = params

    # This is for the plane strain case and alpha = 1 
    # see original paper of Rozhko (2008)
    eta   = (1-2*nu)/(1-nu)/2
    kappa = 3-4*nu 

    # coordinate transform
    rho = sqrt(X[1]^2 + X[2]^2)
    phi = atan(X[2], X[1])

    if rho < r_in
        Pf   = G
        Pt   = 0.
        Ux   = 0.
        Ut   = 0.
        Ur   = 0.
        Uy   = 0.
        Sxx  = 0.  
        Syy  = 0.  
        Sxy  = 0. 
        Srr  = 0.  
        Stt  = 0.  
        Srt  = 0. 
    else
        Srr = (eta*rho^2*dPf*m^3*cos(2*phi)*log(1/(rho^32))+eta*rho^2*dPf*m^4*log(1/(r_out^8))+eta*dPf*m^4*log(rho^8*r_out^8)+eta*rho^6*dPf*log(r_out^8)+eta*m^2*rho^4*P0*log(1/(rho^32))+eta*m^2*rho^4*dPf*log(rho^32)+eta*m*rho^6*P0*cos(2*phi)*log(rho^32)+eta*m^2*rho^2*dPf*log(1/(r_out^8))+eta*rho^2*P0*m^4*log(r_out^8)+eta*rho^8*dPf*log(1/r_out^8*rho^8)+eta*rho^8*P0*log(1/rho^8*r_out^8)+8*eta*P0*m^4+24*eta*m^2*rho^2*dPf+16*eta*m^2*rho^4*P0-16*eta*m^2*rho^4*dPf-8*eta*rho^6*dPf*m^2+8*eta*rho^6*P0*m^2-8*eta*dPf*m^4-24*eta*m^2*rho^2*P0+eta*rho^6*dPf*m^2*log(r_out^8)+8*eta*rho^2*dPf*m^4-8*eta*rho^2*P0*m^4+8*eta*rho^2*dPf*m^3*cos(2*phi)-8*eta*m^3*dPf*cos(2*phi)+8*eta*m^3*P0*cos(2*phi)+eta*P0*m^4*log(1/(rho^8*r_out^8))+eta*rho^6*P0*log(1/(r_out^8))+24*eta*m*rho^4*P0*cos(2*phi)-8*eta*m^2*rho^2*P0*cos(4*phi)+8*eta*m^2*rho^4*P0*cos(4*phi)-24*eta*m*rho^6*P0*cos(2*phi)+8*eta*m^2*rho^2*dPf*cos(4*phi)-8*eta*m^2*rho^4*dPf*cos(4*phi)-8*eta*rho^2*P0*m^3*cos(2*phi)+eta*m^2*rho^4*dPf*log(rho^16)*cos(4*phi)+eta*m^2*rho^4*P0*log(1/(rho^16))*cos(4*phi)-24*eta*m*rho^4*dPf*cos(2*phi)+eta*rho^2*P0*m^3*cos(2*phi)*log(rho^32)+24*eta*m*rho^6*dPf*cos(2*phi)+eta*m^2*rho^2*P0*log(r_out^8)+eta*rho^6*P0*m^2*log(1/(r_out^8))+eta*m*rho^6*dPf*cos(2*phi)*log(1/(rho^32)))/(m^2*rho^4*cos(4*phi)*log(r_out^16)+m^2*rho^4*log(r_out^32)+rho^6*m*cos(2*phi)*log(1/(r_out^32))+rho^2*m^3*cos(2*phi)*log(1/(r_out^32))+rho^8*log(r_out^8)+m^4*log(r_out^8));
        Stt = (eta*dPf*m^4*log(rho^8*r_out^8)+eta*m^2*rho^4*P0*log(1/(rho^32))+eta*m^2*rho^4*dPf*log(rho^32)+eta*m^2*rho^2*dPf*log(r_out^8)+eta*rho^6*P0*m^2*log(r_out^8)+eta*rho^2*dPf*m^4*log(r_out^8)+eta*m^2*rho^2*P0*log(1/(r_out^8))+eta*rho^8*dPf*log(1/r_out^8*rho^8)+eta*rho^8*P0*log(1/rho^8*r_out^8)+16*eta*P0*m^4-24*eta*m^2*rho^2*dPf+16*eta*m^2*rho^4*P0-16*eta*m^2*rho^4*dPf+8*eta*rho^6*dPf*m^2-8*eta*rho^6*P0*m^2-16*eta*dPf*m^4+24*eta*m^2*rho^2*P0-8*eta*rho^2*dPf*m^4+8*eta*rho^2*P0*m^4+56*eta*rho^2*dPf*m^3*cos(2*phi)+8*eta*m^3*dPf*cos(2*phi)-8*eta*m^3*P0*cos(2*phi)+eta*P0*m^4*log(1/(rho^8*r_out^8))-24*eta*m*rho^4*P0*cos(2*phi)+8*eta*m^2*rho^2*P0*cos(4*phi)+8*eta*m^2*rho^4*P0*cos(4*phi)+24*eta*m*rho^6*P0*cos(2*phi)-8*eta*m^2*rho^2*dPf*cos(4*phi)-8*eta*m^2*rho^4*dPf*cos(4*phi)-56*eta*rho^2*P0*m^3*cos(2*phi)+eta*m^2*rho^4*dPf*log(rho^16)*cos(4*phi)+eta*m^2*rho^4*P0*log(1/(rho^16))*cos(4*phi)+24*eta*m*rho^4*dPf*cos(2*phi)-24*eta*m*rho^6*dPf*cos(2*phi)+eta*rho^6*dPf*log(1/(r_out^8))+eta*rho^6*dPf*m^2*log(1/(r_out^8))+eta*rho^2*P0*m^3*cos(2*phi)*log(rho^32*r_out^32)+eta*rho^6*P0*log(r_out^8)+eta*rho^2*P0*m^4*log(1/(r_out^8))+eta*m*rho^6*P0*cos(2*phi)*log(1/r_out^32*rho^32)+8*eta*rho^8*dPf-8*eta*rho^8*P0+eta*m*rho^6*dPf*cos(2*phi)*log(1/rho^32*r_out^32)+eta*rho^2*dPf*m^3*cos(2*phi)*log(1/(rho^32*r_out^32)))/(m^2*rho^4*cos(4*phi)*log(r_out^16)+m^2*rho^4*log(r_out^32)+rho^6*m*cos(2*phi)*log(1/(r_out^32))+rho^2*m^3*cos(2*phi)*log(1/(r_out^32))+rho^8*log(r_out^8)+m^4*log(r_out^8));
        Srt = eta*m*sin(2*phi)*(-2*rho^6*dPf*log(r_out)+2*rho^2*log(r_out)*P0*m^2-2*rho^4*log(r_out)*P0*m^2+2*rho^4*dPf*log(r_out)*m^2+2*m*dPf*rho^2*cos(2*phi)-2*m*P0*rho^2*cos(2*phi)-2*rho^4*dPf*m^2+2*rho^6*log(r_out)*P0+2*rho^4*P0*m^2-2*m*rho^4*dPf*cos(2*phi)-m^2*dPf+m^2*P0-3*rho^2*P0*m^2+3*rho^4*P0-3*rho^4*dPf+3*rho^2*dPf*m^2-3*rho^6*P0+3*rho^6*dPf+2*m*rho^4*P0*cos(2*phi)+2*rho^4*dPf*log(r_out)-2*rho^2*dPf*log(r_out)*m^2-2*rho^4*log(r_out)*P0)/log(r_out)/(4*m^2*rho^4*cos(2*phi)^2+2*m^2*rho^4-4*rho^6*m*cos(2*phi)-4*rho^2*m^3*cos(2*phi)+rho^8+m^4);
                
        Ux  = -1/8*eta*r_in*cos(phi)*(11*m*rho^4*dPf-11*m*rho^4*P0+kappa*rho^6*dPf+4*m^3*log(rho)*P0-4*m^3*log(rho)*dPf+5*rho^2*P0*m^2-kappa*rho^6*P0+4*rho^2*P0*m^3-3*kappa*m^3*dPf+3*kappa*m^3*P0-4*rho^2*dPf*m^3+12*m*P0*rho^2-12*rho^2*m*dPf-4*rho^2*m*log(r_out)*P0+2*kappa*log(r_out)*rho^6*P0+4*rho^2*dPf*log(r_out)*m^3+4*kappa*log(r_out)*m^3*dPf-20*rho^4*m*dPf*cos(phi)^2+6*rho^4*m*log(r_out)*P0+20*m*P0*rho^4*cos(phi)^2+12*dPf*m^2*cos(phi)^2*rho^2-16*m*P0*cos(phi)^2*rho^2-12*P0*m^2*cos(phi)^2*rho^2+16*m*dPf*cos(phi)^2*rho^2+dPf*m^3+4*P0*m^2-rho^6*P0+rho^6*dPf-5*rho^2*dPf*m^2-8*kappa*log(r_out)*m*P0*rho^4*cos(phi)^2+8*kappa*log(r_out)*m^2*dPf*rho^2+2*kappa*log(r_out)*m*P0*rho^4-2*kappa*log(r_out)*m^2*P0*rho^2-4*kappa*rho^4*dPf*m*cos(phi)^2+4*kappa*log(r_out)*m*dPf*rho^4+4*kappa*m*P0*rho^4*cos(phi)^2-8*rho^4*m*log(r_out)*P0*cos(phi)^2+16*cos(phi)^2*dPf*rho^2*log(rho)*m^2+4*rho^4*dPf*m^2-4*rho^4*P0*m^2-5*kappa*rho^2*dPf*m^2+5*kappa*m^2*P0*rho^2+16*cos(phi)^2*m*log(rho)*dPf*rho^4-16*cos(phi)^2*rho^2*P0*log(rho)*m^2-16*cos(phi)^2*m*log(rho)*P0*rho^4+8*kappa*log(r_out)*m^2*P0*cos(phi)^2*rho^2-16*kappa*log(r_out)*m^2*dPf*cos(phi)^2*rho^2-12*kappa*m^2*P0*cos(phi)^2*rho^2-16*log(r_out)*m^2*dPf*cos(phi)^2*rho^2+4*rho^6*dPf*log(r_out)+4*rho^4*log(r_out)*P0-4*rho^4*dPf*log(r_out)+8*log(r_out)*m^2*P0*cos(phi)^2*rho^2+12*kappa*m^2*dPf*cos(phi)^2*rho^2-2*rho^6*log(r_out)*P0-4*dPf*m^2-4*rho^4*dPf*log(r_out)*m^2+4*rho^4*log(r_out)*P0*m^2+12*rho^2*dPf*log(r_out)*m^2-6*rho^2*log(r_out)*P0*m^2-P0*m^3+kappa*m*P0*rho^4-kappa*rho^4*dPf*m-12*dPf*rho^2*log(rho)*m^2+12*rho^2*P0*log(rho)*m^2-4*rho^2*log(r_out)*P0*m^3-12*m*log(rho)*dPf*rho^4+12*m*log(rho)*P0*rho^4-2*kappa*log(r_out)*m^3*P0+4*rho^2*m*dPf*log(r_out)+4*rho^6*P0*log(rho)-4*dPf*rho^6*log(rho)+2*log(r_out)*m^3*P0)/rho/log(r_out)/G/(-m^2+4*m*rho^2*cos(phi)^2-2*m*rho^2-rho^4);       
        Uy  = -1/8*eta*r_in*sin(phi)*(-9*m*rho^4*dPf+9*m*rho^4*P0-kappa*rho^6*dPf+4*m^3*log(rho)*P0-4*m^3*log(rho)*dPf+7*rho^2*P0*m^2+kappa*rho^6*P0+4*rho^2*P0*m^3-3*kappa*m^3*dPf+3*kappa*m^3*P0-4*rho^2*dPf*m^3-4*m*P0*rho^2+4*rho^2*m*dPf-4*rho^2*m*log(r_out)*P0-2*kappa*log(r_out)*rho^6*P0+4*rho^2*dPf*log(r_out)*m^3+4*kappa*log(r_out)*m^3*dPf+20*rho^4*m*dPf*cos(phi)^2-2*rho^4*m*log(r_out)*P0-20*m*P0*rho^4*cos(phi)^2+12*dPf*m^2*cos(phi)^2*rho^2+16*m*P0*cos(phi)^2*rho^2-12*P0*m^2*cos(phi)^2*rho^2-16*m*dPf*cos(phi)^2*rho^2+dPf*m^3-4*P0*m^2+rho^6*P0-rho^6*dPf-7*rho^2*dPf*m^2+8*kappa*log(r_out)*m*P0*rho^4*cos(phi)^2+8*kappa*log(r_out)*m^2*dPf*rho^2-6*kappa*log(r_out)*m*P0*rho^4-6*kappa*log(r_out)*m^2*P0*rho^2+4*kappa*rho^4*dPf*m*cos(phi)^2+4*kappa*log(r_out)*m*dPf*rho^4-4*kappa*m*P0*rho^4*cos(phi)^2+8*rho^4*m*log(r_out)*P0*cos(phi)^2+16*cos(phi)^2*dPf*rho^2*log(rho)*m^2-4*rho^4*dPf*m^2+4*rho^4*P0*m^2-7*kappa*rho^2*dPf*m^2+7*kappa*m^2*P0*rho^2-16*cos(phi)^2*m*log(rho)*dPf*rho^4-16*cos(phi)^2*rho^2*P0*log(rho)*m^2+16*cos(phi)^2*m*log(rho)*P0*rho^4+8*kappa*log(r_out)*m^2*P0*cos(phi)^2*rho^2-16*kappa*log(r_out)*m^2*dPf*cos(phi)^2*rho^2-12*kappa*m^2*P0*cos(phi)^2*rho^2-16*log(r_out)*m^2*dPf*cos(phi)^2*rho^2-4*rho^6*dPf*log(r_out)-4*rho^4*log(r_out)*P0+4*rho^4*dPf*log(r_out)+8*log(r_out)*m^2*P0*cos(phi)^2*rho^2+12*kappa*m^2*dPf*cos(phi)^2*rho^2+2*rho^6*log(r_out)*P0+4*dPf*m^2+4*rho^4*dPf*log(r_out)*m^2-4*rho^4*log(r_out)*P0*m^2+4*rho^2*dPf*log(r_out)*m^2-2*rho^2*log(r_out)*P0*m^2-P0*m^3+5*kappa*m*P0*rho^4-5*kappa*rho^4*dPf*m-4*dPf*rho^2*log(rho)*m^2+4*rho^2*P0*log(rho)*m^2-4*rho^2*log(r_out)*P0*m^3+4*m*log(rho)*dPf*rho^4-4*m*log(rho)*P0*rho^4-2*kappa*log(r_out)*m^3*P0+4*rho^2*m*dPf*log(r_out)-4*rho^6*P0*log(rho)+4*dPf*rho^6*log(rho)+2*log(r_out)*m^3*P0)/rho/log(r_out)/G/(m^2-4*m*rho^2*cos(phi)^2+2*m*rho^2+rho^4);

        Ur  =  1/8*r_in*eta*(-4*rho^2*dPf*log(r_out)*m^2+4*m^2*log(rho)*dPf-4*m^2*log(rho)*P0-4*rho^2*log(r_out)*m*P0*cos(2*phi)-2*rho^4*log(r_out)*P0+4*rho^4*dPf*log(r_out)+4*rho^2*log(r_out)*P0*m^2+4*kappa*log(r_out)*m*dPf*rho^2*cos(2*phi)+4*rho^2*log(r_out)*m*dPf*cos(2*phi)+2*kappa*log(r_out)*m^2*P0-4*kappa*log(r_out)*m^2*dPf+2*kappa*log(r_out)*rho^4*P0-4*m*P0*cos(2*phi)-4*kappa*rho^2*dPf*m*cos(2*phi)+4*kappa*rho^2*P0*m*cos(2*phi)-2*log(r_out)*m^2*P0-3*kappa*P0*m^2+kappa*rho^4*dPf+4*rho^2*log(r_out)*P0+3*kappa*dPf*m^2-4*dPf*rho^2*log(r_out)-kappa*rho^4*P0+4*m*dPf*cos(2*phi)+P0*m^2-dPf*m^2-8*m*dPf*rho^2*cos(2*phi)+8*m*P0*rho^2*cos(2*phi)-4*rho^4*dPf*log(rho)+4*rho^4*P0*log(rho)-4*kappa*log(r_out)*rho^2*P0*m*cos(2*phi)-rho^4*P0+rho^4*dPf-4*rho^2*P0*m^2+4*rho^2*dPf*m^2)/(-2*m*rho^2*cos(2*phi)+rho^4+m^2)^(1/2)/rho/G/log(r_out);
        Ut  = -1/4*r_in*eta*m*sin(2*phi)*(2*dPf*rho^2*log(r_out)-kappa*rho^2*dPf+rho^2*kappa*P0+2*kappa*log(r_out)*dPf*rho^2-rho^2*P0+rho^2*dPf-2*dPf+2*P0-4*rho^2*dPf*log(rho)+4*rho^2*P0*log(rho))/(-2*m*rho^2*cos(2*phi)+rho^4+m^2)^(1/2)/rho/G/log(r_out);

        Pf  = P0 + dPf - dPf*log(rho)/log(r_out);
        Sxx =  1/2*(((-2*rho.^2+1+rho.^4).*Srr+(-2*rho.^2-1-rho.^4).*Stt).*cos(2*phi)+(-2*rho.^2+1+rho.^4).*Srr+(2*rho.^2+1+rho.^4).*Stt+(-2*rho.^4 .*sin(2*phi)+2*sin(2*phi)).*Srt)./(-2*rho.^2 .*cos(2*phi)+rho.^4+1);
        Syy = -1/2*(((2*rho.^2+1+rho.^4).*Srr+(-rho.^4+2*rho.^2-1).*Stt).*cos(2*phi)+(-2*rho.^2-1-rho.^4).*Srr+(-rho.^4+2*rho.^2-1).*Stt+(-2*rho.^4 .*sin(2*phi)+2*sin(2*phi)).*Srt)./(-2*rho.^2 .*cos(2*phi)+rho.^4+1);
        Sxy =  1/2*((2+2*rho.^4).*Srt.*cos(2*phi)+(-sin(2*phi)+rho.^4 .*sin(2*phi)).*Srr+(sin(2*phi)-rho.^4 .*sin(2*phi)).*Stt-4*Srt.*rho.^2)./(-2*rho.^2 .*cos(2*phi)+rho.^4+1);
        Pt  = -1/2*(Sxx + Syy) 
    end   
    return (u=[Ux Uy], u_pol=[Ur Ut], pt=Pt, pf=Pf, σ=[Sxx Sxy; Sxy Syy], σ_pol=[Srr Srt; Srt Stt])
end

function Poroelasticity2D_Rozhko2008(coords::Union{Tuple, NamedTuple};
    params = (r_in=1.0, r_out=20.0, P0=0.0, dPf=1.0, m=0.0, nu= 0.49, G=1.0) )
    X = SVector(values(coords)...)
    sol = Poroelasticity2D_Rozhko2008(X; params)
    return (pt=sol.pt, pf=sol.pf, u=(x=sol.u[1], y=sol.u[2]), σ=(xx=sol.σ[1,1], xy=sol.σ[1,2], yx=sol.σ[2,1], yy=sol.σ[2,2]) )
end
