using ExactFieldSolutions, Plots

function main()

    params = (
        θ1  =  30π/180, 
        θ2  = 180π/180, 
        η1  = 1e0, 
        η2  = 1e3, 
        Vr0 = -1., 
        Vt0 =  0., 
        Vr2 =  0., 
        Vt2 =  0.,
    )

    # Define domain
    Lx  = 10.
    Ly  = 5.
    Nx  = 400
    Ny  = Int64(round(abs(Ly)/(Lx)*Nx))
    x   = LinRange(-Lx/2, Lx/2, Nx)
    y   = LinRange( 0, Ly, Ny)

    # Allocate arrays
    p    = zeros(Nx, Ny)
    Vx   = zeros(Nx, Ny)
    Vy   = zeros(Nx, Ny)
    V    = zeros(Nx, Ny)
    ε̇II  = zeros(Nx, Ny)
   
    # Evaluate solution
    for i=1:Nx, j=1:Ny
        sol       = Stokes2D_Moulas2021( [x[i]; y[j]]; params )
        V[i,j]    = sqrt(sol.V[1]^2 + sol.V[2]^2)
        p[i,j]    = sol.p
        Vx[i,j]   = sol.V[1]
        Vy[i,j]   = sol.V[2]
        ε̇II[i,j]  = sqrt( 1/2*(sol.ε̇[1,1]^2 + sol.ε̇[2,2]^2) + sol.ε̇[1,2]^2  )
    end
    
    # Visualise
    p1 = heatmap(x, y, p',   ylims=(0,Ly), aspect_ratio=1, color=:vik, clims=(-20,20), title="p - Fig. 4 Moulas et al. (2021)", titlefontsize=10)
    p2 = heatmap(x, y, ε̇II', ylims=(0,Ly), aspect_ratio=1, color=:vik, clims=(0,1), title="ε̇II", titlefontsize=10)
    p3 = heatmap(x, y, Vx',  ylims=(0,Ly), aspect_ratio=1 , color=:vik, title="Vx", titlefontsize=10)
    p4 = heatmap(x, y, Vy',  ylims=(0,Ly), aspect_ratio=1 , color=:vik, title="Vy", titlefontsize=10)
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()