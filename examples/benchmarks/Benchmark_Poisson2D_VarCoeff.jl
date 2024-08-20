using Plots, SparseArrays, Symbolics, SparseDiffTools, Printf, StaticArrays, ExactFieldSolutions
import LinearAlgebra: norm
import Statistics: mean

ExactSolution(x) = Poisson2D_VarCoeff(x)

function Residual_Poisson2D!(F, T, f, k, Δ, Xc)
    ncx, ncy = size(T,1), size(T,2)
    for j = 1:ncy, i=1:ncx

        # Interpolate coefficient (second order) 
        kW  = 0.5*(k[i+1,j+1] + k[i,  j+1])
        kE  = 0.5*(k[i+1,j+1] + k[i+2,j+1])
        kWS = 0.5*(k[i+1,j+0] + k[i,  j+0])
        kES = 0.5*(k[i+1,j+0] + k[i+2,j+0])
        kWN = 0.5*(k[i+1,j+2] + k[i,  j+2])
        kEN = 0.5*(k[i+1,j+2] + k[i+2,j+2])
        
        kS  = 0.5*(k[i+1,j+1] + k[i+1,  j])
        kN  = 0.5*(k[i+1,j+1] + k[i+1,j+2])
        kSW = 0.5*(k[i+0,j+1] + k[i+0,  j])
        kNW = 0.5*(k[i+0,j+1] + k[i+0,j+2])
        kSE = 0.5*(k[i+2,j+1] + k[i+2,  j])
        kNE = 0.5*(k[i+2,j+1] + k[i+2,j+2])
    
        ############################

        if (i==1)
            # West boundary
            sol   = ExactSolution( @SVector([Xc.x[i]-Δ.x; Xc.y[j]]) )
            fW    = sol.s
            TW    = sol.u
            sol   = ExactSolution( @SVector([Xc.x[i]-Δ.x; Xc.y[j]-Δ.y]) )
            TSW   = sol.u
            sol   = ExactSolution( @SVector([Xc.x[i]-Δ.x; Xc.y[j]+Δ.y]) )
            TNW   = sol.u
        else
            fW    = f[i-1,j]
            TW    = T[i-1,j]
            if j==1
                sol   = ExactSolution( @SVector([Xc.x[i]-Δ.x; Xc.y[j]-Δ.y]) )
                TSW   = sol.u
            else
                TSW = T[i-1,j-1]
            end
            if j==ncy
                sol   = ExactSolution( @SVector([Xc.x[i]-Δ.x; Xc.y[j]+Δ.y]) )
                TNW   = sol.u
            else
                TNW   = T[i-1,j+1]
            end
        end

        if i==ncx
            # East boundary
            sol   = ExactSolution( @SVector([Xc.x[i]+Δ.x; Xc.y[j]]) )
            fE    = sol.s
            TE    = sol.u
            sol   = ExactSolution( @SVector([Xc.x[i]+Δ.x; Xc.y[j]-Δ.y]) )
            TSE   = sol.u
            sol   = ExactSolution( @SVector([Xc.x[i]+Δ.x; Xc.y[j]+Δ.y]) )
            TNE   = sol.u
        else
            fE    = f[i+1,j]
            TE    = T[i+1,j]
            if j==1
                sol   = ExactSolution( @SVector([Xc.x[i]+Δ.x; Xc.y[j]-Δ.y]) )
                TSE   = sol.u
            else
                TSE = T[i+1,j-1]
            end
            if j==ncy
                sol   = ExactSolution( @SVector([Xc.x[i]+Δ.x; Xc.y[j]+Δ.y]) )
                TNE   = sol.u
            else
                TNE   = T[i+1,j+1]
            end
        end

        if j==1
            # South boundary
            sol   = ExactSolution( @SVector([Xc.x[i]; Xc.y[j]-Δ.y]) )
            fS    = sol.s
            TS    = sol.u
        else
            fS    = f[i,j-1]
            TS    = T[i,j-1]
        end
        
        if j==ncy
            # North boundary
            sol   = ExactSolution( @SVector([Xc.x[i]; Xc.y[j]+Δ.y]) )
            fN    = sol.s
            TN    = sol.u
        else
            fN    = f[i,j+1]
            TN    = T[i,j+1]
        end

        # Central point
        fC = f[i,j]
        TC = T[i,j]

        # Balance equation 
        ii = i + (j-1)*ncx

        # Stencil connections
        u  = @SVector([TSW; TS; TSE; TW; TC; TE; TNW; TN; TNE])

        # # 5-point
        # Mx = @SVector([0.;     0.;          0.;    -kW;   (kW+kE);           -kE;     0.;    0.;           0.]) 
        # My = @SVector([0.;    -kS;          0.;     0.;   (kS+kN);            0.;     0.;   -kN;           0.]) 

        # Skewed 9-point
        Mx = @SVector([-kWS/12;   (kWS+kES)/12;   -kES/12;       -10*kW/12;   10*(kW+kE)/12;     -10*kE/12;    -kWN/12;   (kWN+kEN)/12;    -kEN/12]) 
        My = @SVector([-kSW/12;      -10*kS/12;   -kSE/12;    (kSW+kNW)/12;   10*(kS+kN)/12;  (kSE+kNE)/12;    -kNW/12;      -10*kN/12;    -kNE/12]) 

        # Right-hand side with Laplacian term !
        𝑓 = fC - ( 4*fC - fW - fE - fS - fN )/12

        # Residual
        F[ii] = Mx'*u./Δ.x^2 + My'*u./Δ.y^2 - 𝑓
    end
    return nothing
end

Residual!(F, T, f, k, Δ, Xc) = Residual_Poisson2D!(F, T, f, k, Δ, Xc)

function main(Δ, nc, L)

    # Parameters
    x   = (min=-L.x/2, max=L.x/2)
    xv  = LinRange(x.min,        x.max,        nc.x+1)
    xc  = LinRange(x.min+Δ.x/2., x.max-Δ.x/2., nc.x  )
    xce = LinRange(x.min-Δ.x/2., x.max+Δ.x/2., nc.x+2)
    y   = (min=-L.y/2, max=L.y/2)
    yv  = LinRange(y.min,        y.max,        nc.y+1)
    yc  = LinRange(y.min+Δ.y/2., y.max-Δ.y/2., nc.y  )
    yce = LinRange(y.min-Δ.y/2., y.max+Δ.y/2., nc.y+2)

    Xc  = (x=xc, y=yc)    

    # Allocations
    T   = zeros(nc.x, nc.y) .+ rand(nc.x, nc.y)
    F   = zeros(nc.x*nc.y)
    f   = zeros(nc.x, nc.y)
    δT  = zeros(nc.x*nc.y)
    Ta  = zeros(nc.x, nc.y)
    kc_ex = zeros(nc.x+2, nc.y+2) 

    # Initial condition: Evaluate exact initial solution
    for j=1:nc.y, i=1:nc.x 
        sol     = ExactSolution( @SVector([xc[i]; yc[j]]) )
        Ta[i,j] = sol.u
        f[i,j]  = sol.s 
    end

    # Variable coefficient on vertices
    for j=1:nc.y+2, i=1:nc.x+2 
        sol     = ExactSolution( @SVector([xce[i]; yce[j]]) )
        kc_ex[i,j] = sol.β
    end

    # Sparsity pattern
    input       = rand(nc.x, nc.y)
    output      = similar(input)
    Res_closed! = (F, T) -> Residual!(F, T, f, kc_ex, Δ, Xc)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))

    # Makes coloring
    colors      = matrix_colors(J)

    @show  r = norm(F)/nc.x
    
    # Time loop
    for iter=1:10

        # Residual evaluation: T is found if F = 0
        Residual!(F, T, f, kc_ex, Δ, Xc)
        Res_closed! = (F, T) -> Residual!(F, T, f, kc_ex, Δ, Xc)
        r = norm(F)/length(F)
        @printf("## Iteration %06d: r = %1.2e ##\n", iter, r)
        if r < 1e-10 break end
            
        # Jacobian assembly
        forwarddiff_color_jacobian!(J, Res_closed!, T, colorvec = colors)

        # Solve
        δT   .= .-J\F

        # update
        T    .+= reshape(δT, nc...) 
    end

    p1 = heatmap(xc, yc, (Ta)')
    p2 = heatmap(xc, yc, (T)')
    display(plot(p1, p2))
    display(heatmap(xce, yce, log10.(kc_ex)'))

    # Error
    return mean(abs.(T .- Ta))
end

function ConvergenceAnalysis()

    # Space
    L   = (x=2., y=2.)
    Ncx = [40, 80, 160]  
    Δxv = L.x ./ Ncx
    ϵx  = zero(Δxv)
    for i in eachindex(Ncx)
        Δ   = (x=Δxv[i], y=Δxv[i]) 
        nc  = (x=Ncx[i], y=Ncx[i])
        @show (Δ, nc, L)
        ϵx[i] = main(Δ, nc, L)
    end

    p1 = plot(xlabel="log10(1/Δx)", ylabel="log10(ϵx)")
    p1 = scatter!( log10.( 1.0./Δxv ), log10.(ϵx), label="ϵ")
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 1.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O1"  ) 
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 2.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O2"  ) 
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 4.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O4"  ) 
    display(plot(p1))
end

@time ConvergenceAnalysis()

# nc = (x = 50, y = 50) 
# L  = (x = 2.0, y = 2.0)
# Δ  = (x = L.x/nc.x, y = L.y/nc.y)
# main(Δ, nc, L)

