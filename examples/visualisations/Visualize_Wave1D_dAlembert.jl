using ExactFieldSolutions, Plots, Printf
import LinearAlgebra: norm

function main()
    # Define domain
    Nx   = 500
    Nt   = 5
    x    = LinRange(-1/2, 1/2, Nx)
    xc   = 0.5*(x[2:end] .+ x[1:end-1])
    c    = 1.0
    k    = 8
    t    = LinRange(0, 0.2, Nt)
    u    = zeros(Nt,Nx)

    # Evaluate solution for a several times
    params = (c=c, k=k)
    p = plot(xlabel = "x", ylabel = "u" )
    for it in eachindex(t)
        for i in eachindex(x)
            sol     = Wave1D_dAlembert([x[i],t[it]]; params)
            u[it,i] = sol.u
        end
        plot!(x, u[it,:], label = "t = $(@sprintf("%1.2f", t[it]))")
    end
    display(p)
end

main()