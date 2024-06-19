ElasticModulus(x, ρ, k, σ, Gbg, G0) = Gbg + G0*(1 - 0.5*erfc(x/σ))

function Wave1D_HeteroPlusSource_u_fwd(X, params)
    x, t = X[1], X[2]
    @unpack ρ, k, σ, Gbg, G0 = params
    G = ElasticModulus(x, ρ, k, σ, Gbg, G0)
    c = sqrt(G/ρ)
    return real(exp(im*k*(x - c*t)))
end

function Stress(X, ρ, k, σ, Gbg, G0)
    G      = ElasticModulus(X[1], ρ, k, σ, Gbg, G0)
    params = (ρ=ρ, k=k, σ=σ, Gbg=Gbg, G0=G0) 
    f_cl   = X -> Wave1D_HeteroPlusSource_u_fwd(X, params)
    gradu  = ForwardDiff.gradient(f_cl, X)
    return G*gradu[1]
end

function Wave1D_HeteroPlusSource(X;
    params = (ρ=2.0, k=8.0, σ=0.02, Gbg=1, G0=9.0) )
    x, t   = X[1], X[2]
    u      = Wave1D_HeteroPlusSource_u_fwd(X, params)
    f_cl   = X -> Wave1D_HeteroPlusSource_u_fwd(X, params)
    gradu  = ForwardDiff.gradient(f_cl, X)
    hessu  = ForwardDiff.hessian(f_cl, X)

    σ      = Stress(X, params...)
    S_cl   = X -> Stress(X, params...)
    ∂σ∂x   = ForwardDiff.gradient(S_cl , X)
    s      = params.ρ*hessu[2,2] - ∂σ∂x[1]

    G      =  ElasticModulus(x, params...)
    # hessu  = ForwardDiff.hessian(f_cl, X)
    # s      = hessu[2,2] - params.Gbg/params.ρ*(hessu[1,1])
    # @show s
    return (u=u, ∇u=gradu, s=s, G=G)
end