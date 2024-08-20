using Plots, SparseArrays, Symbolics, SparseDiffTools, Printf, StaticArrays, ExactFieldSolutions
import LinearAlgebra: norm
import Statistics: mean

function Residual_Poisson2D!(F, T, f, k, Δ, Xc)
    ncx, ncy = size(T,1), size(T,2)
    for j = 1:ncy, i=1:ncx

        # Interpolate coefficient (second order) 
        kW = 0.5*(k[i,  j] + k[i,  j+1])
        kE = 0.5*(k[i+1,j] + k[i+1,j+1])
        kS = 0.5*(k[i,  j] + k[i+1,  j])
        kN = 0.5*(k[i,j+1] + k[i+1,j+1])

        ############################

        if (i==1)
            # West boundary
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]-Δ.x; Xc.y[j]]) )
            fW    = sol.s
            TW    = sol.u
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]-Δ.x; Xc.y[j]-Δ.y]) )
            TSW   = sol.u
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]-Δ.x; Xc.y[j]+Δ.y]) )
            TNW   = sol.u
        else
            fW    = f[i-1,j]
            TW    = T[i-1,j]
            if j==1
                sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]-Δ.x; Xc.y[j]-Δ.y]) )
                TSW   = sol.u
            else
                TSW = T[i-1,j-1]
            end
            if j==ncy
                sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]-Δ.x; Xc.y[j]+Δ.y]) )
                TNW   = sol.u
            else
                TNW   = T[i-1,j+1]
            end
        end

        if i==ncx
            # East boundary
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]+Δ.x; Xc.y[j]]) )
            fE    = sol.s
            TE    = sol.u
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]+Δ.x; Xc.y[j]-Δ.y]) )
            TSE   = sol.u
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]+Δ.x; Xc.y[j]+Δ.y]) )
            TNE   = sol.u
        else
            fE    = f[i+1,j]
            TE    = T[i+1,j]
            if j==1
                sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]+Δ.x; Xc.y[j]-Δ.y]) )
                TSE   = sol.u
            else
                TSE = T[i+1,j-1]
            end
            if j==ncy
                sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]+Δ.x; Xc.y[j]+Δ.y]) )
                TNE   = sol.u
            else
                TNE   = T[i+1,j+1]
            end
        end

        if j==1
            # South boundary
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]; Xc.y[j]-Δ.y]) )
            fS    = sol.s
            TS    = sol.u
        else
            fS    = f[i,j-1]
            TS    = T[i,j-1]
        end
        
        if j==ncy
            # North boundary
            sol   = Poisson2D_Sevilla2018( @SVector([Xc.x[i]; Xc.y[j]+Δ.y]) )
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

        # # 5-point classic with variable coefficient
        # qxW = -kW*(TC - TW)/Δ.x # West flux
        # qxE = -kE*(TE - TC)/Δ.x # East flux
        # qyS = -kS*(TC - TS)/Δ.y # South flux
        # qyN = -kN*(TN - TC)/Δ.y # North flux
        # F[ii] = (qxE - qxW)/Δ.x + (qyN - qyS)/Δ.y - f[i,j]
       
        # # 5-point
        # Mx = @SVector([0.;     1.;          0.;    1.;   -4.;           1.;     0.;    1.;           0.]) 

        # Skewed 9-point
        Mx = @SVector([-1/12;     2/12;       -1/12;    -10/12;   20/12;    -10/12;    -1/12;     2/12;      -1/12]) 
        My = @SVector([-1/12;   -10/12;       -1/12;      2/12;   20/12;      2/12;    -1/12;   -10/12;      -1/12]) 

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
    x   = (min=0., max=L.x)
    xc  = LinRange(x.min+Δ.x/2., x.max-Δ.x/2., nc.x)
    y   = (min=0., max=L.y)
    yc  = LinRange(y.min+Δ.y/2., y.max-Δ.y/2., nc.y)
    Xc  = (x=xc, y=yc)    

    # Allocations
    T   = zeros(nc.x, nc.y) .+ rand(nc.x, nc.y)
    F   = zeros(nc.x*nc.y)
    f   = zeros(nc.x, nc.y)
    δT  = zeros(nc.x*nc.y)
    Ta  = zeros(nc.x, nc.y)
    k   = zeros(nc.x+1, nc.y+1) .+ ones(nc.x+1, nc.y+1)

    # Initial condition: Evaluate exact initial solution
    for j=1:nc.y, i=1:nc.x 
        sol     = Poisson2D_Sevilla2018( @SVector([xc[i]; yc[j]]) )
        Ta[i,j] = sol.u
        f[i,j]  = sol.s 
    end

    # Sparsity pattern
    input       = rand(nc.x, nc.y)
    output      = similar(input)
    Res_closed! = (F, T) -> Residual!(F, T, f, k, Δ, Xc)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))

    # Makes coloring
    colors      = matrix_colors(J)

    @show  r = norm(F)/nc.x
    
    # Time loop
    for iter=1:10

        # Residual evaluation: T is found if F = 0
        Residual!(F, T, f, k, Δ, Xc)
        Res_closed! = (F, T) -> Residual!(F, T, f, k, Δ, Xc)
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

    display(heatmap(xc, yc, (T.-0*Ta)'))

    # Error
    return mean(abs.(T .- Ta))
end

function ConvergenceAnalysis()

    # Space
    L   = (x=1., y=1.)
    Ncx = [20, 40, 80]  
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

# nc = (x = 40, y = 40) 
# L  = (x = 1.0, y = 1.0)
# Δ  = (x = L.x/nc.x, y = L.y/nc.y)
# main(Δ, nc, L)

