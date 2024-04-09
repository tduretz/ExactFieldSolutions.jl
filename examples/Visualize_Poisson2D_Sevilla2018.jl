using ExactSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 9, 9
    x    = LinRange(0, 1, Nx)
    y    = LinRange(0, 1, Ny)
    xc   = 0.5*(x[2:end] .+ x[1:end-1]) 
    yc   = 0.5*(y[2:end] .+ y[1:end-1]) 

    # Allocate arrays
    u    = zeros(Nx-1, Ny-1)
    s    = zeros(Nx-1, Ny-1)
    ∂u∂x = zeros(Nx-1, Ny-1)
    ∂u∂y = zeros(Nx-1, Ny-1)

    # Evaluate solution
    for i=1:Nx-1, j=1:Ny-1
        sol       = Poisson2D_Sevilla2018( [xc[i]; yc[j]] )
        u[i,j]    = sol.u
        s[i,j]    = sol.s
        ∂u∂x[i,j] = sol.∇u[1]
        ∂u∂y[i,j] = sol.∇u[2]
    end
    @show s[1,1], s[end,end]
    
    # Visualise
    p1 = heatmap(x, y, u',    aspect_ratio=1, xlims=(0,1), color=:jet, title="u")
    p2 = heatmap(x, y, s',    aspect_ratio=1, xlims=(0,1), color=:jet, title="s")
    p3 = heatmap(x, y, ∂u∂x', aspect_ratio=1, xlims=(0,1), color=:jet, title="∂u∂x")
    p4 = heatmap(x, y, ∂u∂y', aspect_ratio=1, xlims=(0,1), color=:jet, title="∂u∂y")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()