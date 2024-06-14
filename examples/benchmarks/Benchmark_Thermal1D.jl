using Plots, SparseArrays, Symbolics, SparseDiffTools, Printf, ExactFieldSolutions
import LinearAlgebra: norm
import Statistics: mean

function Residual_BackwardEuler!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t)
    ncx = length(F)
    params = (T0=1.0, K=1.0, σ=0.1)
    for i in eachindex(T)
        # West boundary
        if i==1
            sol   = Diffusion1D_Gaussian([x.min-Δx/2.0; t]; params )
            Ta = sol.u
            qxW = -k*(T[i] - Ta)/Δx
        end
        # East boundary
        if i==ncx
            sol   = Diffusion1D_Gaussian([x.min+Δx/2.0; t]; params )
            Ta = sol.u
            qxE = -k*(Ta - T[i])/Δx
        end
        # West flux
        if i>1 # Left flux 
            qxW = -k*(T[i] - T[i-1])/Δx
        end
        # East
        if i<ncx # Right flux
            qxE = -k*(T[i+1] - T[i])/Δx
        end
        # Balance
        F[i] = ρ*Cp*(T[i] - T0[i])/Δt + (qxE - qxW)/Δx 
    end
    return nothing
end

function Residual_CrankNicolson!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t, θ)
    ncx = length(F)
    params = (T0=1.0, K=1.0, σ=0.1)
    for i in eachindex(T)
        # West boundary
        if i==1
            # Previous step
            sol   = Diffusion1D_Gaussian([x.min-Δx/2.0; t-Δt]; params )
            Ta    = sol.u
            qxW0  = -k*(T0[i] - Ta)/Δx
            # Next step
            sol   = Diffusion1D_Gaussian([x.min-Δx/2.0; t]; params )
            Ta    = sol.u
            qxW   = -k*(T[i] - Ta)/Δx
        end
        # East boundary
        if i==ncx
            # Previous step
            sol   = Diffusion1D_Gaussian([x.min+Δx/2.0; t-Δt]; params )
            Ta    = sol.u
            qxE0  = -k*(Ta - T0[i])/Δx
            # Next step
            sol   = Diffusion1D_Gaussian([x.min+Δx/2.0; t]; params )
            Ta    = sol.u
            qxE   = -k*(Ta - T[i])/Δx
        end
        # West flux
        if i>1 # Left flux 
            qxW0 = -k*(T0[i] - T0[i-1])/Δx
            qxW  = -k*( T[i] -  T[i-1])/Δx
        end
        # East
        if i<ncx # Right flux
            qxE0 = -k*(T0[i+1] - T0[i])/Δx
            qxE  = -k*( T[i+1] -  T[i])/Δx

        end
        # Balance
        F[i] = ρ*Cp*(T[i] - T0[i])/Δt + θ*(qxE - qxW)/Δx + (1.0-θ)*(qxE0 - qxW0)/Δx
    end
    return nothing
end

Residual!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t)    = Residual_BackwardEuler!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t)
Residual!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t, θ) = Residual_CrankNicolson!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t, θ)

function main(Δx, Δt, ncx, nt, L)

    # Parameters
    x   = (min=-L/2, max=L/2)
    xc  = LinRange(x.min+Δx/2., x.max-Δx/2., ncx)
    ρ   = 1.0
    Cp  = 1.0
    k   = 1.0
    t   = 0.
    θ   = 1.0 # Crank-Nicolson scheme

    # Allocations
    T   = zeros(ncx)
    T0  = zeros(ncx)
    F   = zeros(ncx)
    δT  = zeros(ncx)
    Ta  = zeros(ncx)

    # Sparsity pattern
    input       = rand(ncx)
    output      = similar(input)
    Res_closed! = (F, T) -> Residual!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t, θ)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))

    # Makes coloring
    colors      = matrix_colors(J)

    # Initial condition: Evaluate exact initial solution
    params = (T0=1.0, K=1.0, σ=0.1)
    for i in eachindex(T)
        sol   = Diffusion1D_Gaussian([xc[i]; t]; params )
        T[i] = sol.u
    end
    
    # Time loop
    for it=1:nt
        T0 .= T
        t  += Δt 
        @printf("########### Step %06d ###########\n", it)

        for iter=1:10

            # Residual evaluation: T is found if F = 0
            Residual!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t, θ)
            Res_closed! = (F, T) -> Residual!(F, T, T0, k, ρ, Cp, Δx, Δt, x, t, θ)
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
    end

    # Evaluate exact solution
    for i in eachindex(T)
        sol   = Diffusion1D_Gaussian([xc[i]; t]; params )
        Ta[i] = sol.u
    end

    # Error
    return mean(abs.(T .- Ta))
end

function ConvergenceAnalysis()

    L   = 2.
    tt  = 0.01
    K   = 1.

    # Time
    ncx = 400
    Δx  = L/ncx
    Nt  = [40, 80, 160, 320]  
    Δtv = 0.01 ./ Nt
    ϵt  = zero(Δtv)
    for i in eachindex(Nt)
        Δx    = sqrt(Δtv[i]*sqrt(K)) / 20
        ncx   = Int64(floor(L/Δx))
        L     = ncx*Δx
        ϵt[i] = main(Δx, Δtv[i], ncx, Nt[i], L)
    end

     # Time
     L   = 2.
     nt  = 100
     Ncx = [40, 80, 160, 320]  
     Δxv = 2.0 ./ Ncx
     ϵx  = zero(Δtv)
     for i in eachindex(Ncx)
        Δt    = Δxv[i]^2/sqrt(K)
        nt    = Int64(floor(tt/Δt))
        ϵx[i] = main(Δxv[i], Δt, Ncx[i], nt, L)
     end

    p1 = plot(xlabel="log10(1/Δx)", ylabel="log10(ϵx)")
    p1 = scatter!( log10.( 1.0./Δxv ), log10.(ϵx), label="ϵ")
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 1.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O1"  ) 
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 2.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O2"  ) 

    p2 = plot(xlabel="log10(1/Δt)", ylabel="log10(ϵt)")
    p2 = scatter!( log10.( 1.0./Δtv ), log10.(ϵt), label="ϵ")
    p2 = plot!( log10.( 1.0./Δtv ), log10.(ϵt[1]) .- 1.0* ( log10.( 1.0./Δtv ) .- log10.( 1.0./Δtv[1] ) ), label="O1" ) 
    p2 = plot!( log10.( 1.0./Δtv ), log10.(ϵt[1]) .- 2.0* ( log10.( 1.0./Δtv ) .- log10.( 1.0./Δtv[1] ) ), label="O2" ) 
    display(plot(p1, p2))
end

ConvergenceAnalysis()