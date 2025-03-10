using ExactFieldSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 101, 101
    x    = LinRange(-5, 5, Nx)
    y    = LinRange(-5, 5, Ny)

    # Allocate arrays
    pt   = zeros(Nx, Ny)
    pf   = zeros(Nx, Ny)
    ur   = zeros(Nx, Ny)
    ut   = zeros(Nx, Ny)

    # Evaluate solution
    for i=1:Nx, j=1:Ny
        sol       = Poroelasticity2D_Rozhko2008( [x[i]; y[j]] )
        pt[i,j]   = sol.pt
        pf[i,j]   = sol.pf
        ur[i,j]   = sol.u_pol[1]
        ut[i,j]   = sol.u_pol[2]
    end
    
    # Visualise
    p1 = heatmap(x, y, pf', aspect_ratio=1, xlims=extrema(x), color=(:matter), title="pf")
    p2 = heatmap(x, y, pt', aspect_ratio=1, xlims=extrema(x), color=(:matter), title="pt")
    p3 = heatmap(x, y, ur', aspect_ratio=1, xlims=extrema(x), color=(:matter), title="ur")
    p4 = heatmap(x, y, ut', aspect_ratio=1, xlims=extrema(x), color=(:matter), title="ut")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()