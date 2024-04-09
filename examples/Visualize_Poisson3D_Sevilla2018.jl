using ExactSolutions, Plots

function main()

    # Define domain
    Nx, Ny, Nz = 8, 8, 8
    x    = LinRange(0, 1, Nx)
    y    = LinRange(0, 1, Ny)
    z    = LinRange(0, 1, Ny)

    # Allocate arrays
    u    = zeros(Nx, Ny, Nz)
    s    = zeros(Nx, Ny, Nz)
    ∂u∂x = zeros(Nx, Ny, Nz)
    ∂u∂y = zeros(Nx, Ny, Nz)
    ∂u∂z = zeros(Nx, Ny, Nz)

    # Evaluate solution
    for i=1:Nx, j=1:Ny, k=1:Ny
        sol         = Poisson3D_Sevilla2018( [x[i]; y[j]; z[j]] )
        u[i,j,k]    = sol.u
        s[i,j,k]    = sol.s
        ∂u∂x[i,j,k] = sol.∇u[1]
        ∂u∂y[i,j,k] = sol.∇u[2]
    end
    
    # Visualise
    p1 = heatmap(x, y, u[:,:,Int64(ceil(Nz/2))]',    aspect_ratio=1, xlims=(0,1), color=:jet, title="u - xy middle cut")
    p2 = heatmap(x, y, u[:,Int64(ceil(Ny/2)),:]',    aspect_ratio=1, xlims=(0,1), color=:jet, title="u - xz middle cut")
    p3 = heatmap(x, y, u[Int64(ceil(Nx/2)),:,:]',    aspect_ratio=1, xlims=(0,1), color=:jet, title="u - yz middle cut")
    display( plot(p1,p2,p3, layout=(2,2)) ) 
 
end

main()