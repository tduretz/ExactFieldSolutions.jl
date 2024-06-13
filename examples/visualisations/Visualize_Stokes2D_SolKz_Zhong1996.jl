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
    ρ    = zeros(Nx, Ny)
   
    # Evaluate solution
    for i=1:Nx, j=1:Ny
        sol       = Stokes2D_SolKz_Zhong1996( [x[i]; y[j]] )
        p[i,j]    = sol.p
        Vx[i,j]   = sol.V[1]
        Vy[i,j]   = sol.V[2]
        ρ[i,j]    = sol.ρ
    end
    
    # Visualise
    p1 = heatmap(x, y, p',    aspect_ratio=1, xlims=(0,1.0), color=:jet, title="p")
    p2 = heatmap(x, y, ρ',  aspect_ratio=1,   xlims=(0,1.0), color=:jet, title="ρ")
    p3 = heatmap(x, y, Vx', aspect_ratio=1,   xlims=(0,1.0), color=:jet, title="Vx")
    p4 = heatmap(x, y, Vy', aspect_ratio=1,   xlims=(0,1.0), color=:jet, title="Vy")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()