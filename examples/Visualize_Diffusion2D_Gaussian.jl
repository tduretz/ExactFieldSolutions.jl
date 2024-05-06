using ExactFieldSolutions, Plots, Printf

function main()

    # Define domain
    Nx, Ny = 50, 50
    x    = LinRange(-1/2, 1/2, Nx)
    y    = LinRange(-1/2, 1/2, Ny)
    xc   = 0.5*(x[2:end] .+ x[1:end-1]) 
    yc   = 0.5*(y[2:end] .+ y[1:end-1]) 

    # Allocate arrays
    t    = [0.0 1e3 1e4 1e5]
    u    = zeros(length(t), Nx-1, Ny-1)
    s    = zeros(length(t), Nx-1, Ny-1)
    ∂u∂x = zeros(length(t), Nx-1, Ny-1)
    ∂u∂y = zeros(length(t), Nx-1, Ny-1)

    # Evaluate solution
    for it in eachindex(t)
        for i=1:Nx-1, j=1:Ny-1
            sol          = Diffusion2D_Gaussian( [xc[i]; yc[j]; t[it]] )
            u[it,i,j]    = sol.u
            s[it,i,j]    = sol.s
            ∂u∂x[it,i,j] = sol.∇u[1]
            ∂u∂y[it,i,j] = sol.∇u[2]
        end
    end
    
    # Visualise
    p1 = heatmap(x, y, u[1,:,:]',    aspect_ratio=1, xlims=(-1/2,1/2), color=cgrad(:roma, rev=true), title=@sprintf("u @ t = %1.2e", t[1]))
    p2 = heatmap(x, y, u[2,:,:]',    aspect_ratio=1, xlims=(-1/2,1/2), color=cgrad(:roma, rev=true), title=@sprintf("u @ t = %1.2e", t[2]))
    p3 = heatmap(x, y, u[3,:,:]',    aspect_ratio=1, xlims=(-1/2,1/2), color=cgrad(:roma, rev=true), title=@sprintf("u @ t = %1.2e", t[3]))
    p4 = heatmap(x, y, u[4,:,:]',    aspect_ratio=1, xlims=(-1/2,1/2), color=cgrad(:roma, rev=true), title=@sprintf("u @ t = %1.2e", t[4]))
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 

end

main()