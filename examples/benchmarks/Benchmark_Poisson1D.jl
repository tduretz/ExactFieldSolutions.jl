using Plots, SparseArrays, Symbolics, SparseDiffTools, Printf, ExactFieldSolutions
import LinearAlgebra: norm
import Statistics: mean

function main_test(Δx, ncx, L)

    # Parameters
    x   = (min=-L/2, max=L/2)
    xc  = LinRange(x.min+Δx/2., x.max-Δx/2., ncx)
    xce = LinRange(x.min-Δx/2., x.max+Δx/2., ncx+2)
    xv  = LinRange(x.min, x.max, ncx+1)
    xve = LinRange(x.min-Δx, x.max+Δx, ncx+3)

    # Allocations
    T      = zeros(ncx+2)
    b      = zeros(ncx+2)
    k      =  ones(ncx+3)
    dTdx   =  ones(ncx+3)
    q      =  ones(ncx+1)
    qa     =  ones(ncx+3)
    dqdx   =  ones(ncx+2)
    F_dTdx = zeros(ncx+1)
    F_dqdx = zeros(ncx-0)
    F      = zeros(ncx)
    Ta     = zeros(ncx)

    # Initial condition: Evaluate exact initial solution
    for i in eachindex(T)
        sol     = Poisson1D_VarCoeff([xce[i]] )
        T[i]    = sol.u
        b[i]    = sol.s
        dqdx[i] = sol.s
    end
    for i in eachindex(k)
        sol     = Poisson1D_VarCoeff([xve[i]] )
        k[i]    = sol.β
        dTdx[i] = sol.∇u[1]
        qa[i]   = sol.q[1]
    end

    for i in eachindex(q)
        sol    = Poisson1D_VarCoeff([xv[i]] )
        q[i]   = sol.q[1]
    end
    
    α = 1/22
    a = 12/11
    α = 0.
    a = 1.0

    # T[2:end-1] .= 0.

    for iter=1:300000
        @. F_dTdx = α*dTdx[1:end-2] + dTdx[2:end-1] + α*dTdx[3:end-0] - a*(T[2:end] - T[1:end-1])/Δx
        @. F_dqdx = α*dqdx[1:end-2] + dqdx[2:end-1] + α*dqdx[3:end-0] - a*(q[2:end] - q[1:end-1])/Δx
        @. F      = -dqdx[2:end-1] - b[2:end-1]
        # @. q[2:end-1]      = -k[3:end-2]*dTdx[3:end-2]
        @. q[1:end-0]     = -k[2:end-1]*dTdx[2:end-1]
        @. dTdx[2:end-1] -= F_dTdx ./10
        @. dqdx[2:end-1] -= F_dqdx ./10
        @. T[2:end-1]    += F/10000
        norm(F)<1e-11 ? break : nothing
        mod(iter, 10000)==0 ?  println(norm(F), ' ', norm(F_dTdx), ' ', norm(F_dqdx)) : nothing
    end

    # Evaluate exact solution
    for i in eachindex(Ta)
        sol   = Poisson1D_VarCoeff([xc[i]] )
        Ta[i] = sol.u
    end

    @show T[1]
    @show Ta[1]

    # p1=plot(xc, abs.(T[2:end-1].-Ta))
    # p1=plot(xc, T[2:end-1])
    # p1=plot!(xc, Ta)
    # p1=plot(xv, qa[2:end-1])
    # p1=plot!(xv, q)
    p1=plot(xv, abs.(qa[2:end-1] .-q))
    display(p1)

    # Error
    return mean(abs.(T[2:end-1] .- Ta))
end

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

    # Makes coloring
    colors      = matrix_colors(J)

    # Initial condition: Evaluate exact initial solution
    for i in eachindex(T)
        sol   = Poisson1D_VarCoeff([xc[i]] )
        T[i] = sol.u
        b[i] = sol.s
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
       ϵx[i] = main_test(Δxv[i], Ncx[i], L)
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


# L   = 2.
# Ncx = 10 
# Δx = L ./ Ncx
# main_test(Δx, Ncx, L)