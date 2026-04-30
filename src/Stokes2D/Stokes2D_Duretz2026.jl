# ---- Matrix potentials ----
ϕ_mat(z, k, γ̇, ε̇, ηm, rc)    = -im/2*ηm*γ̇*z      -   (im*γ̇ + 2*ε̇) * k[1] * rc^2 / z^1 
ϕ′_mat(z, k, γ̇, ε̇, ηm, rc)   = -im/2*ηm*γ̇        +   (im*γ̇ + 2*ε̇) * k[1] * rc^2 / z^2
ϕ′′_mat(z, k, γ̇, ε̇, ηm, rc)  =                   - 2*(im*γ̇ + 2*ε̇) * k[1] * rc^2 / z^3
ψ_mat(z, k, γ̇, ε̇, ηm, rc)    = (im*γ̇ - 2*ε̇)*ηm*z -   (im*γ̇ + 2*ε̇) * k[2] * rc^4 / z^3 
ψ′_mat(z, k, γ̇, ε̇, ηm, rc)   = (im*γ̇ - 2*ε̇)*ηm   + 3*(im*γ̇ + 2*ε̇) * k[2] * rc^4 / z^4 
# ---- Inclusion potentials ----
ϕ_inc(z, k, γ̇, ε̇, ηi, rc)    = -im/2*ηi*γ̇*z 
ϕ′_inc(z, k, γ̇, ε̇, ηi, rc)   = -im/2*ηi*γ̇   
ϕ′′_inc(z, k, γ̇, ε̇, ηi, rc)  = 0.0
ψ_inc(z, k, γ̇, ε̇, ηi, rc)    = 2*(im*γ̇ - 2*ε̇) * k[3] * z  - (im*γ̇ + 2*ε̇)* k[4] *ηi*z
ψ′_inc(z, k, γ̇, ε̇, ηi, rc)   = 2*(im*γ̇ - 2*ε̇) * k[3]      - (im*γ̇ + 2*ε̇)* k[4] *ηi

function residual(k, rc, ηm, ηi, κm, κi, γ̇, ε̇)
    θ    = π/3
    z    = rc * exp(im*θ)
    # ---- Matrix at r=rc ----
    ϕm   = ϕ_mat(z, k, γ̇, ε̇, ηm, rc ) 
    ϕm′  = ϕ′_mat(z, k, γ̇, ε̇, ηm, rc) 
    ϕm′′ = ϕ′′_mat(z, k, γ̇, ε̇, ηm, rc) 
    ψm   = ψ_mat(z, k, γ̇, ε̇, ηm, rc ) 
    ψm′  = ψ′_mat(z, k, γ̇, ε̇, ηm, rc) 
    # ---- Inclusion at r=rc ----
    ϕi   = ϕ_inc(z, k, γ̇, ε̇, ηi, rc )
    ϕi′  = ϕ′_inc(z, k, γ̇, ε̇, ηi, rc)
    ϕi′′ = ϕ′′_inc(z, k, γ̇, ε̇, ηi, rc)
    ψi   = ψ_inc(z, k, γ̇, ε̇, ηi, rc )
    ψi′  = ψ′_inc(z, k, γ̇, ε̇, ηi, rc)
    # ---- Displacement at r=rc ----
    um   = 1/(2*ηm)*(κm*ϕm - z*conj(ϕm′) - conj(ψm))
    ui   = 1/(2*ηi)*(κi*ϕi - z*conj(ϕi′) - conj(ψi))
    # ---- Traction at r=rc ----
    Tm   = 4*real(ϕm′) - 2*(conj(z)* ϕm′′ + ψm′)*exp(2*im*θ)
    Ti   = 4*real(ϕi′) - 2*(conj(z)* ϕi′′ + ψi′)*exp(2*im*θ)
    return @SVector([real(um - ui); imag(um - ui); real(Tm - Ti); imag(Tm - Ti)])
end

function constants(rc, ηm, ηi, κm, κi, γ̇, ε̇)
    k        = @SVector ones(4)
    r        = residual(k, rc, ηm, ηi, κm, κi, γ̇, ε̇)
    J        = ForwardDiff.jacobian( x -> residual(x, rc, ηm, ηi, κm, κi, γ̇, ε̇), k)
    k       -= J \ r
    return k
end

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
(V = [-0.20890904349705283, 0.47592372115980935], p = -0.09568944517573741, τ = [-2.0058103303680195 0.12557580379582506; 0.12557580379582506 2.0064482600025246], τzz = -0.0006379296345049162)
```
```julia-repl
julia> Stokes2D_Duretz2026( (0.2, 0.5) )
(p = -0.09568944517573741, V = (x = -0.20890904349705283, y = 0.47592372115980935), τ = (xx = -2.0058103303680195, xy = 0.12557580379582506, yx = 0.12557580379582506, yy = 2.0064482600025246, zz = -0.0006379296345049162))
```
"""
function Stokes2D_Duretz2026(x; 
    params=(ηm=1.0, ηi=100.0, ξm=100.0, ξi=100.0, R=0.1, γ̇=0.0, ε̇=-1.0))
    ηm, ηi, ξm, ξi, rc, γ̇, ε̇= params
    # Complex coordinate
    z  = x[1] + im*x[2]
    # Poisson ratios
    νm = (3*ξm - 2*ηm) / (2*(3*ξm + ηm))
    νi = (3*ξi - 2*ηi) / (2*(3*ξi + ηi))
    # Kolosov numbers
    κm = 3 - 4*νm
    κi = 3 - 4*νi 
    # Constants
    k = constants(rc, ηm, ηi, κm, κi, γ̇, ε̇)
    # Correct for compressible and incompressible limits
    if (x[1]^2 + x[2]^2) > rc^2
        η, κ = ηm, κm
        # As in Schmid & Podladchikov (2004)
        ϕ    = ϕ_mat(z, k, γ̇, ε̇, ηm, rc )   
        ϕ′   = ϕ′_mat(z, k, γ̇, ε̇, ηm, rc)  
        ϕ′′  = ϕ′′_mat(z, k, γ̇, ε̇, ηm, rc) 
        ψ    = ψ_mat(z, k, γ̇, ε̇, ηm, rc )  
        ψ′   = ψ′_mat(z, k, γ̇, ε̇, ηm, rc)  
    else
        # As in Schmid & Podladchikov (2004)
        η, κ = ηi, κi
        ϕ    = ϕ_inc(z, k, γ̇, ε̇, ηi, rc ) 
        ϕ′   = ϕ′_inc(z, k, γ̇, ε̇, ηi, rc)    
        ϕ′′  = ϕ′′_inc(z, k, γ̇, ε̇, ηi, rc) 
        ψ    = ψ_inc(z, k, γ̇, ε̇, ηi, rc )
        ψ′   = ψ′_inc(z, k, γ̇, ε̇, ηi, rc)
    end
    sxx_p_syy         = 4*real(ϕ′)
    sxx_m_syy_p_2isxy = 2*(conj(z)*ϕ′′ + ψ′)
    sxy               = 1/2*imag(sxx_m_syy_p_2isxy)
    syy               = 1/2*real((sxx_p_syy + sxx_m_syy_p_2isxy))
    sxx               = sxx_p_syy - syy
    szz               = νm * sxx_p_syy
    p                 = -1/3*(sxx_p_syy + szz)
    W                 = 1/(2*η) * (κ*ϕ - z*conj(ϕ′) - conj(ψ))
    return (V=@SVector([real(W), imag(W)]), p=p, τ=@SMatrix([sxx+p sxy; sxy syy+p]), τzz=szz+p )
end

function Stokes2D_Duretz2026(coords::Union{Tuple, NamedTuple};
    params = (ηm=1.0, ηi=100.0, ξm=100.0, ξi=100.0, R=0.1, γ̇=0.0, ε̇=-1.0) )
    x = SVector(values(coords)...)
    sol = Stokes2D_Duretz2026(x; params)
    return (p=sol.p, 
    V=(x=sol.V[1], y=sol.V[2]),
    τ=(xx=sol.τ[1,1], xy=sol.τ[1,2], yx=sol.τ[2,1], yy=sol.τ[2,2], zz=sol.τzz))
end