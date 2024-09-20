using ExactFieldSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 100, 100
    x    = LinRange(-1., 1., Nx)
    y    = LinRange(-1., 1., Ny)
    xc   = 0.5*(x[2:end] .+ x[1:end-1]) 
    yc   = 0.5*(y[2:end] .+ y[1:end-1]) 

    # Allocate arrays
    u    = zeros(Nx-1, Ny-1)
    s    = zeros(Nx-1, Ny-1)
    qx   = zeros(Nx-1, Ny-1)
    qy   = zeros(Nx-1, Ny-1)

    # Evaluate solution
    for i=1:Nx-1, j=1:Ny-1
        sol       = Poisson2D_VarCoeff( [xc[i]; yc[j]] )
        u[i,j]    = sol.u
        s[i,j]    = sol.s
        qx[i,j]   = sol.q[1]
        qy[i,j]   = sol.q[2]
    end
    
    # A = -diff( diff(u[:,2:end-1],dims=1), dims=1)./(x[2]-x[1]).^2 .- diff( diff(u[2:end-1,:],dims=2), dims=2)./(x[2]-x[1]).^2
    # A = 2cos.(Array(x[2:end-1])).*cos.(Array(y[2:end-1])')

    # Visualise
    p1 = heatmap(x, y, u',    aspect_ratio=1, xlims=(-1,1), color=:jet, title="u")
    p2 = heatmap(x, y, s',    aspect_ratio=1, xlims=(-1,1), color=:jet, title="s")
    p3 = heatmap(x, y, qx', aspect_ratio=1, xlims=(-1,1), color=:jet, title="∂u∂x")
    # p4 = heatmap(0.5*(x[2:end]+x[1:end-1]), y, (diff(u,dims=1)./(x[2]-x[1]))', aspect_ratio=1, xlims=(-1,1), color=:jet, title="∂u∂x_num", clim=(-.75,.75))
    # p4 = heatmap((x[2:end-1]), y[2:end-1], (A)', aspect_ratio=1, xlims=(-1,1), color=:jet, title="∂u∂x_num", clim=(.5,2))
    p4 = heatmap(x, y, qy', aspect_ratio=1, xlims=(-1,1), color=:jet, title="∂u∂y")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()