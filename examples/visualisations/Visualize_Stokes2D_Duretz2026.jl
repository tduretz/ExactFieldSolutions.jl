using ExactFieldSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 100, 100
    x    = LinRange(-0.5, 0.5, Nx)
    y    = LinRange(-0.5, 0.5, Ny)

    # Allocate arrays
    p    = zeros(Nx, Ny)
    Vx   = zeros(Nx, Ny)
    Vy   = zeros(Nx, Ny)
    τxx  = zeros(Nx, Ny)
    τyy  = zeros(Nx, Ny)
    τxy  = zeros(Nx, Ny)
   
    # Evaluate solution
    for i=1:Nx, j=1:Ny
        sol       = Stokes2D_Duretz2026( [x[i]; y[j]] )
        p[i,j]    = sol.p
        Vx[i,j]   = sol.V[1]
        Vy[i,j]   = sol.V[2]
        τxx[i,j]  = sol.τ[1,1]
        τyy[i,j]  = sol.τ[2,2]
        τxy[i,j]  = sol.τ[1,2]
    end
    
    # Visualise
    p1 = heatmap(x, y,   p',  aspect_ratio=1, xlims=(-0.5,0.5), color=:matter, reverse=true, title=  "p")
    p2 = heatmap(x, y,  Vx',  aspect_ratio=1, xlims=(-0.5,0.5), color=:matter, reverse=true, title= "Vx")
    p3 = heatmap(x, y,  Vy',  aspect_ratio=1, xlims=(-0.5,0.5), color=:matter, reverse=true, title= "Vy")
    p4 = heatmap(x, y, τxx',  aspect_ratio=1, xlims=(-0.5,0.5), color=:matter, reverse=true, title="τxx")
    p5 = heatmap(x, y, τyy',  aspect_ratio=1, xlims=(-0.5,0.5), color=:matter, reverse=true, title="τyy")
    p6 = heatmap(x, y, τxy',  aspect_ratio=1, xlims=(-0.5,0.5), color=:matter, reverse=true, title="τxy")
    display( plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(600,700)) ) 
 
end

main()