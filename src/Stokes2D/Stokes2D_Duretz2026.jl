@doc raw"""
    sol = Stokes2D_Duretz2026(x; params)  

Provides full solution fields for a matrix-inclusion problem using a compressible viscous rheology and general far-field shear condition. 
Extends the analytical solution of [Schmid & Podladchikov (2003)](https://academic.oup.com/gji/article/155/1/269/713923) and [Jaeger, Cook & Zimmerman](https://www.wiley.com/en-us/Fundamentals+of+Rock+Mechanics%2C+4th+Edition-p-9780632057597):

    x      : is the coordinate vector or tuple
    params : optional parameter array, default (ηm=1.0, ηi=100.0, ξm=100.0, ξi=100.0, R=0.1, γ̇=0.0, ε̇=-1.0). 
and returns:

    sol    : tuple containing the solution fields p (pressure), V (velocity vector) and τ (deviatoric stress tensor)

# Examples
```julia-repl
julia> Stokes2D_Duretz2026( [0.2, 0.5] )
(V = [-0.2089090434970527, 0.47592372115980974], p = -0.09568944517573581, τ = [-2.00581033036802 0.12557580379582292; 0.12557580379582292 2.0064482600025246])
```
```julia-repl
julia> Stokes2D_Duretz2026( (0.2, 0.5) )
(p = -0.09568944517573581, V = (x = -0.2089090434970527, y = 0.47592372115980974), τ = (xx = -2.00581033036802, xy = 0.12557580379582292, yx = 0.12557580379582292, yy = 2.0064482600025246))
```
"""

function Stokes2D_Duretz2026(x; 
    params = (ηm=1.0, ηi=100.0, ξm=100.0, ξi=100.0, R=0.1, γ̇=0.0, ε̇=-1.0) )
    ηm, ηi, ξm, ξi, R, γ̇, ε̇ = params
    # Complex coordinate
    z  = x[1] + im*x[2]
    # Shear viscosity ratio
    β  = ηi/ηm
    # Poisson ratios
    νm = (3*ξm - 2*ηm) / (2*(3*ξm + ηm))
    νi = (3*ξi - 2*ηi) / (2*(3*ξi + ηi))
    # Kolosov numbers
    κm = 3 - 4*νm
    κi = 3 - 4*νi 
    # Constants
    A  = 2 *   2* (1 - β) / (β*κm + 1) # multiplied by 2 w.r.t Jaeger & Cook
    C  = 2 *      (β - 1) / (β*κm + 1) # multiplied by 2 w.r.t Jaeger & Cook
    Ci = 2 * (β*(κi + 1)) / (β*κi + 1) # multiplied by 2 w.r.t Jaeger & Cook
    # Far-field principal stress (Jaeger & Cook)
    # σ1  = 2 * ηm * ε̇  
    # Correct for compressible and incompressible limits
    if (x[1]^2 + x[2]^2) > R^2
        η, κ = ηm, κm
        # As in Schmid & Podladchikov (2004)
        ϕ    = -im/2*ηm*γ̇*z      -   (im*γ̇ + 2*ε̇) * (-A)/4 * R^2 / z^1 
        ϕ′   = -im/2*ηm*γ̇        +   (im*γ̇ + 2*ε̇) * (-A)/4 * R^2 / z^2
        ϕ′′  =                   - 2*(im*γ̇ + 2*ε̇) * (-A)/4 * R^2 / z^3
        ψ    = (im*γ̇ - 2*ε̇)*ηm*z -   (im*γ̇ + 2*ε̇) *    C/2 * R^4 / z^3 
        ψ′   = (im*γ̇ - 2*ε̇)*ηm   + 3*(im*γ̇ + 2*ε̇) *    C/2 * R^4 / z^4 
        # As in Jaeger & Cook
        # ϕ  = σ1/4 * ( A*R^2/z^1) 
        # ϕ′ = σ1/4 * (-A*R^2/z^2) 
        # ψ  =-σ1/2 * ( 2*z + C*R^4/z^3)  
    else
        # As in Schmid & Podladchikov (2004)
        η, κ = ηi, κi
        ϕ    = -im/2*ηi*γ̇*z
        ϕ′   = -im/2*ηi*γ̇
        ϕ′′  = 0.0
        ψ    = 2*(im*γ̇ - 2*ε̇) * Ci/4 * z 
        ψ′   = 2*(im*γ̇ - 2*ε̇) * Ci/4 
        # As in Jaeger & Cook
        # ϕ′ = 0.0
        # ϕ  = 0.0
        # ψ  = -σ1/2 * Ci * z
    end
    sxx_p_syy         = 4*real(ϕ′)
    sxx_m_syy_p_2isxy = 2*(conj(z)*ϕ′′ + ψ′)
    sxy               = 1/2*imag(sxx_m_syy_p_2isxy)
    syy               = 1/2*real((sxx_p_syy + sxx_m_syy_p_2isxy))
    sxx               = sxx_p_syy - syy
    szz               = νm * sxx_p_syy
    p                 = -1/3*(sxx_p_syy + szz)
    W                 = 1/(2*η) * (κ*ϕ - z*conj(ϕ′) - conj(ψ))
    return (V=@SVector([real(W), imag(W)]), p=p, τ=@SMatrix([sxx+p sxy; sxy syy+p]))
end

function Stokes2D_Duretz2026(coords::Union{Tuple, NamedTuple};
    params = (ηm=1.0, ηi=100.0, ξm=100.0, ξi=100.0, R=0.1, γ̇=0.0, ε̇=-1.0) )
    X = SVector(values(coords)...)
    sol = Stokes2D_Duretz2026(X; params)
    return (p=sol.p, 
    V=(x=sol.V[1], y=sol.V[2]),
    τ=(xx=sol.τ[1,1], xy=sol.τ[1,2], yx=sol.τ[2,1], yy=sol.τ[2,2]))
end