using ExactSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 100, 100
    x    = LinRange(-0.5, 0.5, Nx)
    y    = LinRange(-0.5, 0.5, Ny)

    # Allocate arrays
    p    = zeros(Nx, Ny)
    τxx  = zeros(Nx, Ny)
    τyy  = zeros(Nx, Ny)
    τxy  = zeros(Nx, Ny)
   
    # Evaluate solution
    for i=1:Nx, j=1:Ny
        sol       = Elasticity2D_Hole( [x[i]; y[j]] )
        p[i,j]    = sol.p
        τxx[i,j]  = sol.τ[1,1]
        τyy[i,j]  = sol.τ[2,2]
        τxy[i,j]  = sol.τ[1,2]
    end
    
    # Visualise
    p1 = heatmap(x, y, p',   aspect_ratio=1, xlims=(-0.5,0.5), color=:jet, title="p")
    p2 = heatmap(x, y, τxy', aspect_ratio=1, xlims=(-0.5,0.5), color=:jet, title="τxy")
    p3 = heatmap(x, y, τxx', aspect_ratio=1, xlims=(-0.5,0.5), color=:jet, title="τxx")
    p4 = heatmap(x, y, τyy', aspect_ratio=1, xlims=(-0.5,0.5), color=:jet, title="τyy")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()