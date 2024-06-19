# function Wave1D_Duru21_u_fwd(X, params)
#     @unpack c, k  = params
#     x, t = X[1], X[2]
#     c = c0 + ϵ*sin(π*x/L)
#     return cos(k*π*t)*sin(n*π*x/L + a0)
#     # return real(exp(im*k*(x[1] - c*x[2])))
# end

# @doc raw"""
#     sol = Diffusion1D_Gaussian(x; params)  

# Evaluates the manufactured solution of an initial 1D wave:

#     x      : is the coordinate and time vector, x[1] is space and x[2] is time 
#     params : optional parameter array, default: (c = 1.0, k = 8.0) 
# and returns:

#     sol    : tuple containing the solution fields u, ∇u and the source term s = -Δu 

# # Examples
# ```julia-repl
# julia> Wave1D_dAlembert([0 0])
# (u = 1.0, ∇u = [0.0 -0.0], s = 0.0)
# ```

# ```julia-repl
# julia> Wave1D_dAlembert([0 0]; params)
# (u = 1.0, ∇u = [0.0 0.0], s = 0.0)

# julia> Wave1D_dAlembert([0 2]; params)
# (u = 0.424179007336997, ∇u = [10.866940344079486 10.866940344079486], s = 0.0)
# ```
# """
# function Wave1D_dAlembert(x;
#     params = (c=1.0, k=8.0) )
#     u      = Wave1D_dAlembert_u_fwd(x, params)
#     f_cl   = x -> Wave1D_dAlembert_u_fwd(x, params)
#     gradu  = ForwardDiff.gradient(f_cl, x)
#     hessu  = ForwardDiff.hessian(f_cl, x)
#     s      = hessu[2,2] - params.c^2*(hessu[1,1])
#     return (u=u, ∇u=gradu, s=s)
# end

using ExactFieldSolutions, Plots, Printf
import LinearAlgebra: norm
using UnPack, ForwardDiff

E0(x, t, rho, sigma, k) = real(-rho .* (sqrt(-(4 * k .^ 2 .* sigma .^ 4 .* x .^ 2 - 4 * im * k .* sigma .^ 4 .* x + 4 * im * k .* sigma .^ 2 .* x .^ 3 - sigma .^ 4 + 6 * sigma .^ 2 .* x .^ 2 - 9 * x .^ 4) .* exp(2 * x .^ 2 ./ sigma .^ 2)) .* (2 * im * k .* sigma .^ 2 .* x + sigma .^ 2 - 3 * x .^ 2) + (4 * k .^ 2 .* sigma .^ 4 .* x .^ 2 - 4 * im * k .* sigma .^ 4 .* x + 8 * im * k .* sigma .^ 2 .* x .^ 3 - sigma .^ 4 + 6 * sigma .^ 2 .* x .^ 2 - 9 * x .^ 4) .* exp(x .^ 2 ./ sigma .^ 2)) ./ (2 * k .^ 2 .* t .^ 2 .* x .^ 4))

function Wave1D_X_u_fwd(X, params)
    x, t = X[1], X[2]
    @unpack ρ, σ, k = params
    E = E0(x, t, params...)
    c = sqrt(E/ρ)
    return real(exp(im*k*(x - c*t)))
end

function Wave1D_X(X;
    params = (ρ=1.0, σ=0.1, k=8.0) )
    x, t   = X[1], X[2]
    u      = Wave1D_X_u_fwd(X, params)
    f_cl   = X -> Wave1D_X_u_fwd(X, params)
    gradu  = ForwardDiff.gradient(f_cl, X)
    hessu  = ForwardDiff.hessian(f_cl, X)
    E      = E0(x, t, params...) * exp(-x^2/params.σ^2)
    σ      = E*gradu[1]
    ∂σ∂x   = ForwardDiff.gradient(f_cl, X)
    s      = hessu[2,2] - ∂σ∂x[1]
    @show s
    return (u=u, ∇u=gradu, s=s)
end

function main()
    # Define domain
    Nx   = 500
    Nt   = 5
    x    = LinRange(-1/2, 1/2, Nx)
    xc   = 0.5*(x[2:end] .+ x[1:end-1])
    c    = 1.0
    k    = 1
    t    = LinRange(0.1, 0.2, Nt)
    u    = zeros(Nt,Nx)

    ρ    = 1.0
    σ    = 0.1
    Eref = E0.(x, t[1], ρ, σ, k)
    E    = Eref .* exp.(-x.^2/σ.^2)
    u    = similar(x)

    # # Evaluate solution for a several times
    params = (ρ=1.0, σ=0.1, k=8.0) 
    # p = plot(xlabel = "x", ylabel = "u" )
    # for it in eachindex(t)
        # for i in eachindex(x)
        #     sol     = Wave1D_X([x[i],t[1]]; params)
        #     u[i] = sol.u
        # end
    #     plot!(x, u[it,:], label = "t = $(@sprintf("%1.2f", t[it]))")
    # end
    @show Eref
    p = plot(x, E)
    display(p)
end

main()