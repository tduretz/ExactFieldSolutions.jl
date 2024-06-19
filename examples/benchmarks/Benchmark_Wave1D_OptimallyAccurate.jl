using Plots, SparseArrays, Symbolics, SparseDiffTools, Printf, ExactFieldSolutions, StaticArrays
import LinearAlgebra: norm
import Statistics: mean

function Residual_OSA!(F, U, U0, U00, E, ρ, Δx, Δt, x, t)
    ncx = length(F)
    params = (c=sqrt(E/ρ), k=8.0)
    for i in eachindex(U)

        # Center point at 3 time levels
        UC00  = U00[i]
        UC0   = U0[i] 
        UC    = U[i] 
        
        # West point
        if i==1 # West boundary
            # Previous previous step
            sol   = Wave1D_dAlembert([x.min-Δx/2.0; t-2*Δt]; params)
            UW00  = sol.u
            # Previous step
            sol   = Wave1D_dAlembert([x.min-Δx/2.0; t-1*Δt]; params)
            UW0   = sol.u
            # Current step
            sol   = Wave1D_dAlembert([x.min-Δx/2.0; t-0*Δt]; params)
            UW    = sol.u
        else # Inside
            UW00 = U00[i-1]
            UW0  = U0[i-1]
            UW   = U[i-1]
        end
       
        # East point
        if i==ncx  # East boundary
            # Previous previous step
            sol   = Wave1D_dAlembert([x.max+Δx/2.0; t-2*Δt]; params)
            UE00  = sol.u
            # Previous step
            sol   = Wave1D_dAlembert([x.max+Δx/2.0; t-1*Δt]; params)
            UE0   = sol.u
            # Current step
            sol   = Wave1D_dAlembert([x.max+Δx/2.0; t-0*Δt]; params)
            UE    = sol.u
        else # Inside
            UE00 = U00[i+1]
            UE0  = U0[i+1]
            UE   = U[i+1]
        end
    
        # Balance

        # # Conventional - style 1:
        # F[i] = ρ*(1*UC -2*UC0 + 1*UC00)/Δt^2 - E*(1*UW0  -2*UC0 + 1*UE0)/Δx^2

        # # Conventional - style 2:
        # Mt = @SVector([0.;   1.;   0.;   0.;  -2.; 0.;  0.; 1.; 0.]) 
        # Mx = @SVector([0.;   0.;   0.;   1.;  -2.; 1.;  0.; 0.; 0.])
        # u  = @SVector([UW00; UC00; UE00; UW0; UC0; UE0; UW; UC; UE])
        # F[i] = ρ/Δt^2*Mt'*u  - E/Δx^2*Mx'u

        # Conventional - style 2:
        Mt = @SVector([1/12; 10/12; 1/12; -2/12; -20/12; -2/12;  1/12; 10/12; 1/12]) 
        Mx = @SVector([1/12; -2/12; 1/12; 10/12; -20/12; 10/12;  1/12; -2/12; 1/12])
        u  = @SVector([UW00; UC00; UE00; UW0; UC0; UE0; UW; UC; UE])
        F[i] = ρ/Δt^2*Mt'*u  - E/Δx^2*Mx'u
    end
    return nothing
end

Residual!(F, U, U0, U00, E, ρ, Δx, Δt, x, t)  = Residual_OSA!(F, U, U0, U00, E, ρ, Δx, Δt, x, t)

function main_Wave1D_OSA(Δx, Δt, ncx, nt, L)

    # Parameters
    x   = (min=-L/2, max=L/2)
    xc  = LinRange(x.min+Δx/2., x.max-Δx/2., ncx)
    ρ   = 1.0
    E   = 1.0
    t   = 0.
    k   = 8.
    c   = sqrt(E/ρ)

    # Allocations
    U   = zeros(ncx)
    U0  = zeros(ncx)
    U00 = zeros(ncx)
    δU  = zeros(ncx)
    Ua  = zeros(ncx)
    F   = zeros(ncx)

    # Sparsity pattern
    input       = rand(ncx)
    output      = similar(input)
    Res_closed! = (F, U) -> Residual!(F, U, U0, U00, E, ρ, Δx, Δt, x, t)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))

    # Makes coloring
    colors      = matrix_colors(J)

    # Initial condition: Evaluate exact initial solution
    params = (c=c, k=k)
    for i in eachindex(U)
        sol   = Wave1D_dAlembert([xc[i]; t   ]; params)
        U[i]  = sol.u
        sol   = Wave1D_dAlembert([xc[i]; t-Δt]; params)
        U0[i] = sol.u
    end
    
    # Time loop
    for it=1:nt
        U00 .= U0
        U0  .= U
        t   += Δt 
        @printf("########### Step %06d ###########\n", it)

        for iter=1:10

            # Residual evaluation: T is found if F = 0
            Residual!(F, U, U0, U00, E, ρ, Δx, Δt, x, t)
            Res_closed! = (F, U) -> Residual!(F, U, U0, U00, E, ρ, Δx, Δt, x, t)
            r = norm(F)/ncx
            @printf("## Iteration %06d: r = %1.2e ##\n", iter, r)
            if r < 1e-9 break end
                
            # Jacobian assembly
            forwarddiff_color_jacobian!(J, Res_closed!, U, colorvec = colors)

            # Solve
            δU   .= .-J\F

            # update
            U    .+= δU
        end
    end

    # Evaluate exact solution
    for i in eachindex(U)
        sol   = Wave1D_dAlembert([xc[i]; t]; params)
        Ua[i] = sol.u
    end

    # p = plot()
    # p = plot!(xc, U,  label = "u" )
    # p = plot!(xc, Ua, label = "ua")
    # display(p)
    # sleep(0.1)

    # Error
    return mean(abs.(U .- Ua))
end

function ConvergenceAnalysis()

    L   = 1.0
    tt  = 0.2

    ρ   = 1.0
    E   = 1.0
    c   = sqrt(E/ρ)

    # Time
    ncx = 400
    Δx  = L/ncx
    Nt  = [20, 20, 40, 80, 160]  
    Δtv = 0.01 ./ Nt
    ϵt  = zero(Δtv)
    for i in eachindex(Nt)
        Δx    = 2.1 * (Δtv[i]*sqrt(E/ρ)) *20 # Somehow dx has to be large for proper dt convergence!
        ncx   = Int64(floor(L/Δx))
        L     = ncx*Δx
        ϵt[i] = main_Wave1D_Conventional(Δx, Δtv[i], ncx, Nt[i], L)
    end

     # Space
     L   = 2.
     nt  = 100
     Ncx = [20, 40, 80, 160]  
     Δxv = 2.0 ./ Ncx
     ϵx  = zero(Δtv)
     for i in eachindex(Ncx)
        Δt    = Δxv[i]/sqrt(E/ρ)/2.1
        nt    = Int64(floor(tt/Δt))
        ϵx[i] = main_Wave1D_Conventional(Δxv[i], Δt, Ncx[i], nt, L)
     end

    p1 = plot(xlabel="log10(1/Δx)", ylabel="log10(ϵx)")
    p1 = scatter!( log10.( 1.0./Δxv ), log10.(ϵx), label="ϵ")
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 1.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O1"  ) 
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 2.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O2"  ) 
    p1 = plot!( log10.( 1.0./Δxv ), log10.(ϵx[1]) .- 4.0* ( log10.( 1.0./Δxv ) .- log10.( 1.0./Δxv[1] ) ), label="O4"  ) 

    p2 = plot(xlabel="log10(1/Δt)", ylabel="log10(ϵt)")
    p2 = scatter!( log10.( 1.0./Δtv ), log10.(ϵt), label="ϵ")
    p2 = plot!( log10.( 1.0./Δtv ), log10.(ϵt[1]) .- 1.0* ( log10.( 1.0./Δtv ) .- log10.( 1.0./Δtv[1] ) ), label="O1" ) 
    p2 = plot!( log10.( 1.0./Δtv ), log10.(ϵt[1]) .- 2.0* ( log10.( 1.0./Δtv ) .- log10.( 1.0./Δtv[1] ) ), label="O2" ) 
    p2 = plot!( log10.( 1.0./Δtv ), log10.(ϵt[1]) .- 4.0* ( log10.( 1.0./Δtv ) .- log10.( 1.0./Δtv[1] ) ), label="O4" ) 

    display(plot(p1, p2))
end


ConvergenceAnalysis()

# begin
#     ρ   = 1.0
#     E   = 1.0
#     t   = 0.
#     k   = 8.
#     c   = sqrt(E/ρ)
#     L   = 1.0
#     t   = 0.2
#     ncx = 500
#     Δx  = L/ncx
#     nt  = 1000
#     Δt  = Δx/sqrt(E/ρ)/2.1
#     nt  = Int64(floor(t/Δt))
#     main(Δx, Δt, ncx, nt, L)
# end