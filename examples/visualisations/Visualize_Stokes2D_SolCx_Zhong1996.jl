using ExactFieldSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 100, 100
    x    = LinRange(0, 1.0, Nx)
    y    = LinRange(0, 1.0, Ny)

    # Allocate arrays
    p    = zeros(Nx, Ny)
    Vx   = zeros(Nx, Ny)
    Vy   = zeros(Nx, Ny)
   
    # Evaluate solution
    for j=1:Ny, i=1:Nx
        sol       = Stokes2D_SolCx_Zhong1996( [x[i], y[j]] )
        p[i,j]    = sol.p
        Vx[i,j]   = sol.V[1]
        Vy[i,j]   = sol.V[2]
    end
    
    # Visualise
    p1 = heatmap(x, y, p',  aspect_ratio=1, xlims=(0,1.0), color=:jet, title="p")
    p2 = heatmap(x, y, Vx', aspect_ratio=1, xlims=(0,1.0), color=:jet, title="Vx")
    p3 = heatmap(x, y, Vy', aspect_ratio=1, xlims=(0,1.0), color=:jet, title="Vy")
    display( plot(p1,p2,p3, layout=(2,2)) ) 
end

main()