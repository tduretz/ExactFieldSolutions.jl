using Plots, SparseArrays, Symbolics, SparseDiffTools, Printf, ExactFieldSolutions
import LinearAlgebra: norm
import Statistics: mean

function Residual_BackwardEuler!(F, T, b, k, Δx, x)
    ncx = length(F)
    bC, bW, bE = 0., 0., 0.
    for i in eachindex(T)
        bC = b[i]
        # West boundary
        if i==1
            sol = Poisson1D_VarCoeff([x.min-Δx/2.0] )
            qxW = -k[i]*(T[i] - sol.u)/Δx
            qfW = -k[i]*(b[i] - sol.s)/Δx
        end
        # East boundary
        if i==ncx
            sol = Poisson1D_VarCoeff([x.max+Δx/2.0] )
            qxE = -k[i+1]*(sol.u - T[i])/Δx
            qfE = -k[i+1]*(sol.s - b[i])/Δx
        end
        # West flux
        if i>1 # Left flux 
            qxW = -k[i]*(T[i] - T[i-1])/Δx
            qfW = -k[i]*(b[i] - b[i-1])/Δx

        end
        # East
        if i<ncx # Right flux
            qxE = -k[i+1]*(T[i+1] - T[i])/Δx
            qfE = -k[i+1]*(b[i+1] - b[i])/Δx
        end
        # Balance
        F[i] = (qxE - qxW)/Δx - b[i] + (qfE - qfW)/Δx * Δx^2/12
    end
    return nothing
end

Residual!(F, T, b, k, Δx, x)    = Residual_BackwardEuler!(F, T, b, k, Δx, x)

function main(Δx, ncx, L)

    # Parameters
    x   = (min=-L/2, max=L/2)
    xc  = LinRange(x.min+Δx/2., x.max-Δx/2., ncx)
    xv  = LinRange(x.min, x.max, ncx+1)

    # Allocations
    T   = zeros(ncx)
    b   = zeros(ncx)
    k   =  ones(ncx+1)
    F   = zeros(ncx)
    δT  = zeros(ncx)
    Ta  = zeros(ncx)

    # Sparsity pattern
    input       = rand(ncx)
    output      = similar(input)
    Res_closed! = (F, T) -> Residual!(F, T, b, k, Δx, x)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))
    colors      = matrix_colors(J)

    # Initial condition: Evaluate exact initial solution
    for i in eachindex(T)
        sol   = Poisson1D_VarCoeff([xc[i]] )
        T[i] = sol.u
        b[i] =-sol.s
    end
    for i in eachindex(k)
        sol   = Poisson1D_VarCoeff([xv[i]] )
        k[i] = sol.β
    end
    
    # Newton iterations
    for iter=1:10

        # Residual evaluation: T is found if F = 0
        Residual!(F, T, b, k, Δx, x)
        Res_closed! = (F, T) -> Residual!(F, T, b, k, Δx, x)
        r = norm(F)/ncx
        @printf("## Iteration %06d: r = %1.2e ##\n", iter, r)
        if r < 1e-10 break end
            
        # Jacobian assembly
        forwarddiff_color_jacobian!(J, Res_closed!, T, colorvec = colors)

        # Solve
        δT   .= .-J\F

        # update
        T    .+= δT
    end

    # Evaluate exact solution
    for i in eachindex(T)
        sol   = Poisson1D_VarCoeff([xc[i]] )
        Ta[i] = sol.u
    end

    p1=plot(xc, T)
    p1=plot!(xc, Ta)
    display(p1)

    # Error
    return mean(abs.(T .- Ta))
end

function ConvergenceAnalysis()

    # Space
    L   = 2.
    # Ncx = [20, 40, 80].*3  
    Ncx = [10, 20, 30].*1 
    Δxv = L ./ Ncx
    ϵx  = zero(Δxv)
    for i in eachindex(Ncx)
       ϵx[i] = main(Δxv[i], Ncx[i], L)
    end

    @show ϵx

    p1 = plot(xlabel="log10(1/Δx)", ylabel="log10(ϵx)")
    p1 = scatter!( log10.( 1.0./Δxv ), log10.(ϵx), label="ϵ")
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 1.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O1"  ) 
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 2.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O2"  ) 
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 4.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O4"  ) 

    display(plot(p1))
end

ConvergenceAnalysis()