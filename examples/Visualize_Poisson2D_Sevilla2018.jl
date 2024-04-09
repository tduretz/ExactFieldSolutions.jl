using ExactSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 8, 8
    x    = LinRange(0, 1, Nx)
    y    = LinRange(0, 1, Ny)

    # Allocate arrays
    u    = zeros(Nx, Ny)
    s    = zeros(Nx, Ny)
    ∂u∂x = zeros(Nx, Ny)
    ∂u∂y = zeros(Nx, Ny)

    # Evaluate solution
    for i=1:Nx, j=1:Ny
        sol       = Poisson2D_Sevilla2018( [x[i]; y[j]] )
        u[i,j]    = sol.u
        s[i,j]    = sol.s
        ∂u∂x[i,j] = sol.∇u[1]
        ∂u∂y[i,j] = sol.∇u[2]
    end
    
    # Visualise
    p1 = heatmap(x, y, u',    aspect_ratio=1, xlims=(0,1), color=:jet, title="u")
    p2 = heatmap(x, y, s',    aspect_ratio=1, xlims=(0,1), color=:jet, title="s")
    p3 = heatmap(x, y, ∂u∂x', aspect_ratio=1, xlims=(0,1), color=:jet, title="∂u∂x")
    p4 = heatmap(x, y, ∂u∂y', aspect_ratio=1, xlims=(0,1), color=:jet, title="∂u∂y")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()