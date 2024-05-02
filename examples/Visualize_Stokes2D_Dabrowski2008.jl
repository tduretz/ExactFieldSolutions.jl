using ExactFieldSolutions, Plots

function main()

    # Define domain
    Nx, Ny = 100, 100
    x    = LinRange(-0.0, 1.0, Nx)
    y    = LinRange(-0.0, 1.0, Ny)

    # Allocate arrays
    p    = zeros(Nx, Ny)
    Vx   = zeros(Nx, Ny)
    Vy   = zeros(Nx, Ny)
    Vm   = zeros(Nx, Ny)
    ε̇II  = zeros(Nx, Ny)

    δ  = 2
    γ  = sqrt( (sqrt(δ) - 1) / (sqrt(δ) + 1))
    ε̇0 = 1.0
    r  = 0.1
   
    # Evaluate solution
    for i=1:Nx, j=1:Ny

        X  = x[i] + im*y[j]
        X̄  = conj(X)
        Xm = X - γ*X̄
        Xp = X + γ*X̄
        X̄m = conj(Xm)
        X̄p = conj(Xp)
        V  = im*ε̇0*r/(4*γ^2) * (γ*(sqrt(Xp - 4*γ) - sqrt(Xm + 4*γ) + Xm - Xp) + sqrt(X̄p - 4*γ) + sqrt(X̄m + 4*γ) + X̄m - X̄p)
        
        Vx[i,j]   = real(V) #- ε̇0/2*y[j]
        Vy[i,j]   = imag(V) #- ε̇0/2*x[j]
        Vm[i,j]   = sqrt( Vx[i,j]^2 + Vy[i,j]^2 )


        if (x[i]^2 + y[j]^2)<r
            Vx[i,j]   = NaN
            Vy[i,j]   = NaN
            Vm[i,j]   = NaN
        end    
        
        # sol       = Stokes2D_Schmid2003( [x[i]; y[j]] )
        # p[i,j]    = sol.p
        # Vx[i,j]   = sol.V[1]
        # Vy[i,j]   = sol.V[2]
        # ε̇II[i,j]  = sqrt( 1/2*(sol.ε̇[1,1]^2 + sol.ε̇[2,2]^2) + sol.ε̇[1,2]^2  )
    end
    
    # Visualise
    p1 = heatmap(x, y, Vm',   aspect_ratio=1, xlims=(-0.0,1.0), color=:jet, title="|V|")
    p2 = heatmap(x, y, ε̇II',  aspect_ratio=1, xlims=(-0.0,1.0), color=:jet, title="ε̇II")
    p3 = heatmap(x, y, Vx',   aspect_ratio=1, xlims=(-0.0,1.0), color=:jet, title="Vx")
    p4 = heatmap(x, y, Vy',   aspect_ratio=1, xlims=(-0.0,1.0), color=:jet, title="Vy")
    display( plot(p1,p2,p3,p4, layout=(2,2)) ) 
 
end

main()