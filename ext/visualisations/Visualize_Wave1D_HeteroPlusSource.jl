using ExactFieldSolutions, Plots, Printf

function main()
    # Define domain
    Nx   = 500
    Nt   = 5
    x    = LinRange(-1/2, 1/2, Nx)
    t    = LinRange(0.1, 0.2, Nt)
    u    = zeros(Nt, Nx)
    s    = zeros(Nt, Nx)

    # Solution parameters
    params = (ρ=1., k=8., σ=0.02, Gbg=1.0, G0=9.0) 

    # Evaluate solution for several times
    p1 = plot(xlabel = "x", ylabel = "u" )
    for it in eachindex(t)
        for i in eachindex(x)
            sol     = Wave1D_HeteroPlusSource([x[i],t[it]]; params)
            u[it,i] = sol.u
            s[it,i] = sol.s
        end
        p1 = plot!(x, u[it,:], label = "t = $(@sprintf("%1.2f", t[it]))")
    end

    # Evaluate source for several times
    p2 = plot(xlabel = "x", ylabel = "f" )
    for it in eachindex(t)
        for i in eachindex(x)
            sol     = Wave1D_HeteroPlusSource([x[i],t[it]]; params)
            u[it,i] = sol.u
            s[it,i] = sol.s
        end
        p2 = plot!(x, s[it,:], label = "t = $(@sprintf("%1.2f", t[it]))")
    end

    # Evaluate coefficient
    G = zeros(Nx)
    for i in eachindex(x)
        sol  = Wave1D_HeteroPlusSource([x[i],t[1]]; params)
        G[i] = sol.G
    end
    p3 = plot(xlabel = "x", ylabel = "G")
    p3 = plot!(x, G, label=:none)
    display(plot(p1, p2, p3))
end

main()