function Wave1D_dAlembert_u_fwd(x, params)
    @unpack c, k  = params
    return real(exp(im*k*(x[1] + c*x[2])))
end

function Wave1D_dAlembert(x;
    params = (c=1.0, k=8.0) )
    u      = Wave1D_dAlembert_u_fwd(x, params)
    f_cl   = x -> Wave1D_dAlembert_u_fwd(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = hessu[2,2] - params.c^2*(hessu[1,1])
    @show s
    return (u=u, ∇u=gradu, s=s)
end

# TO BE MOVED IN VISUALISATION EXAMPLE IF OK
using ForwardDiff

function Visualize_Wave1D_dAlembert()
    # Define domain
    Nx   = 500
    Nt   = 5
    x    = LinRange(-1/2, 1/2, Nx)
    xc   = 0.5*(x[2:end] .+ x[1:end-1])
    c    = 1.0
    k    = 8
    t    = LinRange(0, 0.1, Nt)
    u    = zeros(Nt,Nx)

    # Evaluate solution for a several times
    params = (c=c, k=k)
    p = plot()
    for it in eachindex(t)
        for i in eachindex(x)
            sol     = Wave1D_dAlembert([x[i],t[it]]; params)
            u[it,i] = sol.u
        end
        plot!(x, u[it,:], label = "t = $(@sprintf("%1.2f", t[it]))")
    end
    display(p)
end

Visualize_Wave1D_dAlembert()