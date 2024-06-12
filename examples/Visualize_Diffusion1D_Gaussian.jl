using ExactFieldSolutions, Plots, Printf
import LinearAlgebra: norm

function main()

    # Define domain
    Nx   = 50
    x    = LinRange(-1/2, 1/2, Nx)
    xc   = 0.5*(x[2:end] .+ x[1:end-1]) 

    # Allocate arrays
    t    = [0.0 1e3 1e4 1e5]
    u    = zeros(length(t), Nx-1)
    s    = zeros(length(t), Nx-1)
    ∂u∂x = zeros(length(t), Nx-1)

    # Evaluate solution
    for it in eachindex(t)
        for i=1:Nx-1
            sol        = Diffusion1D_Gaussian( [xc[i]; t[it]] )
            u[it,i]    = sol.u
            s[it,i]    = sol.s
            ∂u∂x[it,i] = sol.∇u[1]
        end
        @printf("Net source term should be 0.0: norm(s) = %1.2e\n", norm(s[it,:,:]))
    end
    
    # Visualise
    p1 = plot(xlabel = "x", ylabel= "u")
    p1 = plot!(xc,    u[1,:], title=@sprintf("u @ t = %1.2e", t[1]), label="u₀")
    p1 = plot!(xc,    u[1,:], title=@sprintf("u @ t = %1.2e", t[1]), label="u")
    p2 = plot(xlabel = "x", ylabel= "u")
    p2 = plot!(xc,    u[1,:], title=@sprintf("u @ t = %1.2e", t[1]), label="u₀")
    p2 = plot!(xc,    u[2,:], title=@sprintf("u @ t = %1.2e", t[2]), label="u")
    p3 = plot(xlabel = "x", ylabel= "u")
    p3 = plot!(xc,    u[1,:], title=@sprintf("u @ t = %1.2e", t[1]), label="u₀")
    p3 = plot!(xc,    u[3,:], title=@sprintf("u @ t = %1.2e", t[3]), label="u")
    p4 = plot(xlabel = "x", ylabel= "u")
    p4 = plot!(xc,    u[1,:], title=@sprintf("u @ t = %1.2e", t[1]), label="u₀")
    p4 = plot!(xc,    u[4,:], title=@sprintf("u @ t = %1.2e", t[4]), label="u")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 

end

main()