@doc raw"""
    sol = Stokes2D_Moutzouris_circle(x; params)  

    or

    sol = Stokes2D_Moutzouris_ellipse(x; params)

Provide full solution fields for a host-inclusion problem using a compressible viscous or elastic rheology and general far-field shear or volumetric condition. 
Extends the analytical solution of [Schmid & Podladchikov (2003)](https://academic.oup.com/gji/article/155/1/269/713923) and [Jaeger, Cook & Zimmerman](https://www.wiley.com/en-us/Fundamentals+of+Rock+Mechanics%2C+4th+Edition-p-9780632057597):

    x      : is the coordinate vector or tuple
    params : optional parameter array, default (ηm=1.0, ηi=0.1, ξm=1.0, ξi=10.0, ri=0.1, t=2.0, α=0.0, ε̇=1.0, γ̇=0.0, ζ̇=0.0, ε̇zz=0.0). 
and returns:

    sol    : tuple containing the solution fields p (pressure), V (velocity vector) and τ (deviatoric stress tensor)

# Examples
```julia-repl
julia> Stokes2D_Moutzouris_circle( [0.2, 0.5] )
(V = [-0.20890904349705283, 0.47592372115980935], p = -0.09568944517573741, τ = [-2.0058103303680195 0.12557580379582506; 0.12557580379582506 2.0064482600025246], τzz = -0.0006379296345049162)
```
```julia-repl
julia> Stokes2D_Moutzouris_circle( (0.2, 0.5) )
(p = -0.09568944517573741, V = (x = -0.20890904349705283, y = 0.47592372115980935), τ = (xx = -2.0058103303680195, xy = 0.12557580379582506, yx = 0.12557580379582506, yy = 2.0064482600025246, zz = -0.0006379296345049162))
```
"""

# New analytical solution for circle (including OOP)
function Stokes2D_Moutzouris_circle(x; 
    params = (ηm=1.0, ηi=0.1, ξm=1.0, ξi=10.0, ri=0.1, t=2.0, α=0.0, ε̇=1.0, γ̇=0.0, ζ̇=0.0, ε̇zz=0.0))
    (; ηm, ηi, ξm, ξi, ri, t, α, ε̇, γ̇, ζ̇, ε̇zz) = params
    x, y = x[1], x[2]
    r1, r2, sc = ri, ri, 1.0
    νm = (3ξm - 2ηm) / (2*(3ξm + ηm))
    νi = (3ξi - 2ηi) / (2*(3ξi + ηi))
    κm = 3 - 4νm
    κi = 3 - 4νi
    Em = 2ηm*(1 + νm)
    Ei = 2ηi*(1 + νi)

    P_mR = 2*ηm*(ζ̇ + νm*ε̇zz) / (κm - 1)
    P_mI = -ηm*γ̇ / (κm + 1)
    P_m  = P_mR + im*P_mI
    Q_m  = (im*γ̇ - 2ε̇)*ηm

    D   = 2ηi + (κi - 1)*ηm
    B1R = ((κm + 1)*ηi*P_mR + 2*ηi*ηm*(νi - νm)*ε̇zz) / D
    B1I = (κm + 1)*ηi*P_mI / ((κi + 1)*ηm)
    B1  = B1R + im*B1I

    A1  = (ηi - ηm)*conj(Q_m)*ri^2 / (κm*ηi + ηm)
    A3  = ri^2 * A1
    B2  = conj(A1/ri^2 + conj(Q_m))
    A2  = 2*ri^2 * (B1R - P_mR)

    z   = x + im*y
    z̄   = conj(z)
    r   = abs(z)
    if r > ri
        ϕ   =  P_m*z  + A1/z
        ϕ′  =  P_m    - A1/z^2
        ϕ′′ =           2A1/z^3
        ψ   =  Q_m*z  + A2/z   + A3/z^3
        ψ′  =  Q_m    - A2/z^2 - 3A3/z^4
        η_loc, κ_loc, ν_loc, E_loc = ηm, κm, νm, Em
    else
        ϕ   =  B1*z
        ϕ′  =  B1  + 0im
        ϕ′′ =  0.0 + 0im
        ψ   =  B2*z
        ψ′  =  B2  + 0im
        η_loc, κ_loc, ν_loc, E_loc = ηi, κi, νi, Ei
    end
    Qval = z̄*ϕ′′ + ψ′
    sxx  = 2real(ϕ′) - real(Qval)
    syy  = 2real(ϕ′) + real(Qval)
    sxy  = imag(Qval)
    szz  = ν_loc*(sxx + syy) + E_loc*ε̇zz
    p    = -(sxx + syy + szz) / 3
    vel  = (κ_loc*ϕ - z*conj(ϕ′) - conj(ψ)) / (2*η_loc) - ν_loc*ε̇zz*z
    # Velocity gradient tensor
    ∂v∂z = (κ_loc*ϕ′ - conj(ϕ′)) / (2*η_loc) - ν_loc*ε̇zz
    ∂v∂z̄ = -conj(Qval) / (2*η_loc)
    ∂v∂x = ∂v∂z + ∂v∂z̄
    ∂v∂y = im*(∂v∂z - ∂v∂z̄)
    return (V   = @SVector([real(vel), imag(vel)]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p,
            L   = @SMatrix([real(∂v∂x)  real(∂v∂y);
                            imag(∂v∂x)  imag(∂v∂y)]))
end

# physical position from the conformal coordinate
joukowski(λ) = λ + 1.0/λ          # ω(λ) = z   (change here if your map differs)
t_to_ri(t) = sqrt((t - 1.0)*(t + 1.0)) / (t - 1.0)

# Inverse mapping of the Joukowsky transform
function inv_joukowski(z)
    disc = sqrt(z^2 - 4.0 + 0im)
    ζ1   = (z + disc)/2
    ζ2   = (z - disc)/2
    return abs(ζ1) >= abs(ζ2) ? ζ1 : ζ2
end
to_zeta(X) = (ζ = inv_joukowski(complex(X[1], X[2])); @SVector([real(ζ), imag(ζ)]))

# Function to acquire the axes of the inclusion based on t (
function ellipse_axes(t)
    ri = t_to_ri(t)
    a  = ri + 1.0/ri
    b  = ri - 1.0/ri
    return a, b
end

# New analytical solution for ellipse (including OOP)
function Stokes2D_Moutzouris_ellipse(x; 
    params= (ηm=1.0, ηi=0.1, ξm=1.0, ξi=10.0, ri=0.1, t=2.0, α=0.0, ε̇=1.0, γ̇=0.0, ζ̇=0.0, ε̇zz=0.0))
    (; ηm, ηi, ξm, ξi, ri, t, α, ε̇, γ̇, ζ̇, ε̇zz) = params
    r1, r2     = ellipse_axes(t) # true physical semi-axes (a >= 2 always)
    sc         = r2 / ri 
    ri      = t_to_ri(t)
    Ζ       = to_zeta(sc.*x)
    τ, σ    = Ζ[1], Ζ[2]
    # Kolosov
    νm      = (3.0ξm - 2.0ηm) / (2.0*(3.0ξm + ηm))
    νi      = (3.0ξi - 2.0ηi) / (2.0*(3.0ξi + ηi))
    κm      = 3 - 4νm
    κi      = 3 - 4νi
    Em      = 2.0*ηm*(1.0 + νm)
    Ei      = 2.0*ηi*(1.0 + νi)
    # Boundary terms with angle of rotation
    P_mR    = 2.0*ηm*(ζ̇ + νm*ε̇zz) / (κm - 1.0)
    P_mI    = -ηm*γ̇ / (κm + 1.0)
    P_m     = P_mR + im*P_mI
    Q_m     = (im*γ̇ - 2.0ε̇)*ηm*exp(-2.0*im*α)
    Q_mR    = real(Q_m); Q_mI = imag(Q_m)
    # Radius factor
    r4      = ri^4
    # Out-of-plane forcing of the velocity-continuity condition
    Δzz     = 2.0*ηi*ηm*ε̇zz*(νm - νi)
    # Solution coefficients
    K       = -2.0*ηi*κm*P_mR + 2.0*ηi*P_mR + 2.0*ηi*Q_mR*ri^2 + 2.0*ηm*κi*P_mR +
                ηm*κi*Q_mR*ri^2 - 2.0*ηm*P_mR - ηm*Q_mR*ri^2 + 2.0*Δzz
    L       = ηi*κm*P_mR*r4 - ηi*P_mR - ηi*Q_mR*ri^2 + ηm*P_mR*r4 + ηm*P_mR + ηm*Q_mR*ri^2
    M       = ηi*κm - ηi - ηm*κi + ηm
    DEN     = 2.0*ηi^2*κm*r4 - 2.0*ηi^2*κm + ηi*ηm*κi*κm*r4 + ηi*ηm*κi -
                ηi*ηm*κm*r4 + 2.0*ηi*ηm*κm + 2.0*ηi*ηm*r4 - ηi*ηm +
                ηm^2*κi*r4 - ηm^2*κi - ηm^2*r4 + ηm^2
    den1    = ηi*κm*r4 + ηi + ηm*r4 - ηm
    den2    = ηi*κm*r4 - ηi + ηm*r4 + ηm
    A1_R    = (ηi - ηm)*(r4 - 1.0)*K / DEN
    A1_I    = Q_mI*ri^2*(ηm - ηi)*(r4 - 1.0) / den1
    A1      = A1_R + im*A1_I
    A2_R    = 2.0*(r4 - 1.0)*(M*L - Δzz*den2) / (ri^2 * DEN)
    A2_I    = 0.0
    A2      = A2_R + im*A2_I
    A3_R    = (ηi - ηm)*(r4 - 1.0)*(r4 + 1.0)*K / (ri^2 * DEN)
    A3_I    = Q_mI*(ηm - ηi)*(ri^8 - 1.0) / den1
    A3      = A3_R + im*A3_I
    B1_R    = (ηi*(κm + 1.0)*L - Δzz*den1) / DEN
    B1_I    = ηi*(κm + 1.0)*(ηi*κm*P_mI*r4 + ηi*P_mI + ηi*Q_mI*ri^2 +
                ηm*P_mI*r4 - ηm*P_mI - ηm*Q_mI*ri^2) / (ηm*(κi + 1.0)*den1)
    B1      = B1_R + im*B1_I
    B2_R    = ηi*ri^2*(κm + 1.0)*K / DEN
    B2_I    = ηi*Q_mI*r4*(κm + 1.0) / den1
    B2      = B2_R + im*B2_I
    λ       = τ + im*σ
    r       = abs(λ)
    # Potentials for matrix and inclusion
    if r > ri
        ϕ   =  P_m*(λ + 1.0/λ) + A1/λ
        ϕ′  =  P_m - P_m*λ^-2.0 - A1*λ^-2.0
        ϕ′′ =  2.0*P_m*λ^-3.0 + 2.0*A1*λ^-3.0
        ψ   =  Q_m*(λ + 1.0/λ) + A2/λ + A3*(1.0/(λ^3.0 - λ))
        ψ′  =  Q_m - Q_m*λ^-2.0 - A2*λ^-2.0 - A3*(3.0*λ^2 - 1.0)/((λ^3.0 - λ)^2.0)
        η_loc, κ_loc, ν_loc, E_loc = ηm, κm, νm, Em
    else
        ϕ   =  B1*(λ + 1.0/λ)
        ϕ′  =  B1 - B1*λ^-2.0
        ϕ′′ =  2.0*B1*λ^-3.0
        ψ   =  B2*(λ + 1.0/λ)
        ψ′  =  B2 - B2*λ^-2.0
        η_loc, κ_loc, ν_loc, E_loc = ηi, κi, νi, Ei
    end
    # Field formulas
    ω        = λ + 1.0/λ                       # = z (scaled)
    ω′       = 1.0 - λ^-2.0
    ω′′      = 2.0*λ^-3.0
    zbar     = conj(λ) + 1.0/conj(λ)           # = conj(z)
    Φ        = ϕ′/ω′                           # = Φ(z)  = dϕ/dz
    Φp       = (ω′*ϕ′′ - ω′′*ϕ′) / ω′^3        # = Φ′(z) = d²ϕ/dz²
    S        = 4.0*real(Φ)                     # σxx+σyy
    D        = 2.0*( zbar*Φp + ψ′/ω′ )         # σyy-σxx+2iσxy
    sxx      = (S - real(D))/2
    syy      = (S + real(D))/2
    sxy      = imag(D)/2
    szz      = ν_loc*(sxx + syy) + E_loc*ε̇zz
    p        = -(sxx + syy + szz) / 3.0
    conj_ωpr = 1.0 - 1.0/conj(λ)^2
    vel      = (κ_loc*ϕ - ω/conj_ωpr*conj(ϕ′) - conj(ψ)) / (2*η_loc) - ν_loc*ε̇zz*ω

    # Velocity gradient
    ∂v∂z = (κ_loc*Φ - conj(Φ)) / (2*η_loc) - ν_loc*ε̇zz
    ∂v∂z̄ = -conj(D) / (4*η_loc)
    ∂v∂x = ∂v∂z + ∂v∂z̄
    ∂v∂y = im*(∂v∂z - ∂v∂z̄)

    return (V   = @SVector([real(vel)/sc, imag(vel)/sc]),
            p   = p,
            τ   = @SMatrix([sxx+p  sxy; sxy  syy+p]),
            τzz = szz + p,
            L   = @SMatrix([real(∂v∂x)  real(∂v∂y);
                            imag(∂v∂x)  imag(∂v∂y)]))
end
