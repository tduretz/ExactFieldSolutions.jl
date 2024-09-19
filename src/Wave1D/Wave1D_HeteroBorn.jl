using Plots, ForwardDiff, UnPack, SpecialFunctions, StaticArrays
import LinearAlgebra: norm

# ElasticModulus(x, ρ, k, σ, Gbg, G0) = Gbg + G0*(1 - 0.5*erfc(x/σ))
ElasticModulus(x, ρ, k, σ, Gbg, G0) = Gbg + G0*0.5*erfc(x/σ)

# ElasticModulus(x, ρ, k, σ, Gbg, G0) = Gbg + G0*sin(2*π*x)/10

params = (ρ=1.0, k=6.0, σ=0.1, Gbg=1, G0=1.0)

function Wave1D_HeteroBorn_u_fwd(X; 
    params  = params)
    x, t, k̄ = X[1], X[2], X[3]
    @unpack ρ, k, σ, Gbg, G0 = params
    G = ElasticModulus(x, ρ, k, σ, Gbg, G0)
    c = sqrt(G/ρ)
    return real(exp(im*k*k̄*(x - c*t)))
end

function Stress(X, ρ, k, σ, Gbg, G0)
    G      = ElasticModulus(X[1], ρ, k, σ, Gbg, G0)
    params = (ρ=ρ, k=k, σ=σ, Gbg=Gbg, G0=G0) 
    f_cl   = X -> Wave1D_HeteroBorn_u_fwd(X; params)
    gradu  = ForwardDiff.gradient(f_cl, X)
    return G*gradu[1]
end

function Force(X;
    params  = params )
    x, t, k̄ = X[1], X[2], X[3]
    u       = Wave1D_HeteroBorn_u_fwd(X; params)
    f_cl    = X -> Wave1D_HeteroBorn_u_fwd(X; params)
    gradu   = ForwardDiff.gradient(f_cl, X)
    hessu   = ForwardDiff.hessian(f_cl, X)

    σ      = Stress(X, params...)
    σ_cl   = X -> Stress(X, params...)
    ∂σ∂x   = ForwardDiff.gradient(σ_cl , X)

    # Residual to minimise
    return  params.ρ*hessu[2,2] - ∂σ∂x[1]
end

function Wave1D_HeteroBorn(X;
    params  = params )
    x, t, k̄ = X[1], X[2], X[3]
    u       = Wave1D_HeteroBorn_u_fwd(X; params)
    f_cl    = X -> Wave1D_HeteroBorn_u_fwd(X; params)
    gradu   = ForwardDiff.gradient(f_cl, X)
    hessu   = ForwardDiff.hessian(f_cl, X)

    σ      = Stress(X, params...)
    σ_cl   = X -> Stress(X, params...)
    ∂σ∂x   = ForwardDiff.gradient(σ_cl , X)

    # Residual to minimise
    s      = Force(X; params)#

    # dsdk̄
    s_cl   = X -> Force(X; params)
    ∂s∂x   = ForwardDiff.gradient(s_cl, X)
    
    G      =  ElasticModulus(x, params...)
    return (u=u, ∇u=gradu, s=s, G=G, ρ=params.ρ, ∂s∂x=∂s∂x )
end

function kfun(x, t, rho, k, sigma, Gbg, G0) 
    G = ElasticModulus(x, rho, k, sigma, Gbg, G0)
    return imag(-im * rho .* (4 * pi ^ 3 * G0 .* sigma .^ 2 .* exp(x .^ 2 ./ sigma .^ 2) .* erfc(x ./ sigma) + sqrt(2) * pi ^ (5 // 2) * G0 .* sigma .* t .* sqrt((G0 .* erfc(x ./ sigma) + 2 * Gbg) ./ rho) + 2 * sqrt(2) * pi ^ 3 * G0 .* t .* x .* sqrt((G0 .* erfc(x ./ sigma) + 2 * Gbg) ./ rho) .* exp(x .^ 2 ./ sigma .^ 2) .* erfc(x ./ sigma) + 8 * pi ^ 3 * Gbg .* sigma .^ 2 .* exp(x .^ 2 ./ sigma .^ 2) + 4 * sqrt(2) * pi ^ 3 * Gbg .* t .* x .* sqrt((G0 .* erfc(x ./ sigma) + 2 * Gbg) ./ rho) .* exp(x .^ 2 ./ sigma .^ 2)) ./ (k .* sigma .* t .* (pi ^ (5 // 2) * G0 .^ 2 .* t .* erfc(x ./ sigma) + 2 * pi ^ (5 // 2) * G0 .* Gbg .* t + 2 * sqrt(2) * pi ^ 3 * G0 .* rho .* sigma .* sqrt((G0 .* erfc(x ./ sigma) + 2 * Gbg) ./ rho) .* exp(x .^ 2 ./ sigma .^ 2) .* erfc(x ./ sigma) + 4 * sqrt(2) * pi ^ 3 * Gbg .* rho .* sigma .* sqrt((G0 .* erfc(x ./ sigma) + 2 * Gbg) ./ rho) .* exp(x .^ 2 ./ sigma .^ 2))))
end
function main()

    Nx   = 600
    Nt   = 5
    x    = LinRange(-1/2, 1/2, Nx)
    t    = LinRange(0.1, 0.2, Nt)
    u    = zeros(Nx)
    F    = zeros(Nx)
    G    = zeros(Nx)
    k̄_num = zeros(Nx)
    k̄_ana = zeros(Nx)

    # u = d'Alembert
    k̄ = 1.0 # d'Alembert

    p1 = plot()
    p2 = plot()
    p3 = plot()

    for it =1:Nt
        k̄ = 1.0 # d'Alembert
        for i in eachindex(F)
            sol  = Wave1D_HeteroBorn( @SVector([x[i], t[it], k̄]) )
            F[i] = sol.s
            iter = 0
            F1 = F[i]
            while abs(F[i])>1e-9 && abs(F[i])/abs(F1)>1e-13 && iter<100000
                δk̄   = -F[i]/sol.∂s∂x[3]
                k̄   += δk̄
                sol  = Wave1D_HeteroBorn( @SVector([x[i], t[it], k̄]) )
                F[i] = sol.s
                iter += 1 
            end
            k̄_num[i] = k̄
            k̄_ana[i] = kfun(x[i], t[it], params...)
            u[i] = Wave1D_HeteroBorn_u_fwd(@SVector([x[i], t[it], k̄]) )
            G[i] = sol.G

            # sol  = Wave1D_HeteroBorn( @SVector([x[i], t[it], k̄_ana[i]]) )
            # @show sol.s

        end
        @show norm(F)
        p2 = plot!(x, u, label=:none, title="u")
    end
    p1 = plot(x, F, label=:none, title="F (should be 0)")

    p3 = plot(x, G, label=:none, title="G")

    p4 = plot(x, k̄_num, label=:none, title="k̄_num")
    # p4 = plot!(x, k̄_ana, label=:none, title="k̄_ana")


    display(plot(p1, p2, p3, p4))

end

main()


