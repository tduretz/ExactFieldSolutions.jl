using ExactFieldSolutions, Plots, Printf
import LinearAlgebra: norm

function main()
    # Define domain
    Nx   = 500
    Nt   = 5
    x    = LinRange(-1/2, 1/2, Nx)
    G    = 1.
    ρ    = 1.
    k    = 8.
    t    = LinRange(0, 0.2, Nt)
    u    = zeros(Nt,Nx)

    # Evaluate solution for a several times
    params = (ρ=ρ, k=k, G=G)
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