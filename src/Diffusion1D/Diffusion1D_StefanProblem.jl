function Calc_λ1(L,c,Tm,T0)
    # guess
    λ1   = 0.5
    # function
    f(λ1,L,c,Tm,T0) = (exp(-λ1^ 2)) / (λ1 * erf(λ1)) - (L * sqrt(pi)) / (c * (Tm - T0))
    for iter=1:100 
        # Residuals
        𝑓        = f(λ1,L,c,Tm,T0)
        # closure
        f_cl     = λ1 -> f(λ1,L,c,Tm,T0)
        d𝑓dλ1    = ForwardDiff.derivative(f_cl, λ1)
        λ1      -= 𝑓/d𝑓dλ1
        abs(𝑓) < 1e-8 ? break : nothing
    end
    return λ1
end

@doc raw"""
    sol = Diffusion1D_StefanProblem(X; params)
    
Evaluates the analytical solution of the Stefan problem in 1D: 

    X      : is the coordinate and time vector or tuple, X[1] is space and X[2] is time 
    params : optional parameter array, default: (Tm=1050, T0=0, L=4e5, c=1e3 , κ=7e-6) 
and returns:

    sol    : tuple containing the solution fields T and the ym the depth of the solidification interface

# Examples
```julia-repl
julia> Diffusion1D_StefanProblem([1 1e5])
(T = 805.5540353379434, ym = 1.465974602531472)
```
```julia-repl
julia> Diffusion1D_StefanProblem( (coord_y=0.1, time=1e5 ) )
(T = 90.13122886071284, ym = 1.465974602531472)
```
"""
function Diffusion1D_StefanProblem(X;
    params = (Tm=1050, T0=0 , L=4e5, c=1e3 , κ=7e-6 ) )
    @unpack Tm, T0, L, c, κ = params
    y, t = X[1], X[2]
    # Compute lambda
    λ1 = Calc_λ1(L, c, Tm, T0)       
    # Calculation of the depth of the solid-liquid boundary
    ym = 2 * λ1 * sqrt(κ * t)         
    # Calculation of dimensionless coodinate η
    η  = y / (2 * sqrt(κ * t))        
    # Calculation of dimensionless temperature θ
    θ  = (erf(η)) / (erf(λ1))
    # Check
    y >= ym ? θ = 1.0 : nothing # ternary operator
    # Calculation of temperature with the dimentionneless temperature (θ)
    T  = (Tm - T0) * θ + T0
    return (T=T, ym=ym)
end

function Diffusion1D_StefanProblem(coords::Union{Tuple, NamedTuple};
    params = (Tm=1050, T0=0 , L=4e5, c=1e3 , κ=7e-6 ) )
    X = SVector(values(coords)...)
    sol = Diffusion1D_StefanProblem(X; params)
    return (T=sol.T, ym=sol.ym)
end

Diffusion1D_StefanProblem( (coord_y=0.1, time=1e5 ) )