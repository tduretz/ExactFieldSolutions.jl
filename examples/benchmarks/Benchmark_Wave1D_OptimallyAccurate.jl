using Plots, SparseArrays, Symbolics, SparseDiffTools, Printf, ExactFieldSolutions, StaticArrays
import LinearAlgebra: norm
import Statistics: mean

# Dynamic dispatch based on a rule
_type(rule) = Val{rule}()
Analytics(rule, x)                 = _func(_type(rule), x) 
_func(::Val{:dAlembert}       , x) = Wave1D_dAlembert(x)
_func(::Val{:HeteroPlusSource}, x) = Wave1D_HeteroPlusSource(x)

# Problem type
problem = :dAlembert
problem = :HeteroPlusSource

noisy =  false

function Residual_OSOA!(F, U, U0, U00, f, G, ρ, Δx, Δt, x, t)
    ncx = length(F)
    for i in eachindex(U)

        # Center point at 3 time levels
        UC00  = U00[i]
        UC0   = U0[i] 
        UC    = U[i] 

        GW    = G[i]
        GE    = G[i+1]
        
        # West point
        if i==1 # West boundary
            # Previous previous step
            sol   = Analytics(problem, @SVector([x.min-Δx/2.0; t-2*Δt]))
            UW00  = sol.u
            # Previous step
            sol   = Analytics(problem, @SVector([x.min-Δx/2.0; t-1*Δt]))
            UW0   = sol.u
            # Current step
            sol   = Analytics(problem, @SVector([x.min-Δx/2.0; t-0*Δt]))
            UW    = sol.u
        else # Inside
            UW00 = U00[i-1]
            UW0  = U0[i-1]
            UW   = U[i-1]
        end
       
        # East point
        if i==ncx  # East boundary
            # Previous previous step
            sol   = Analytics(problem, @SVector([x.max+Δx/2.0; t-2*Δt]))
            UE00  = sol.u
            # Previous step
            sol   = Analytics(problem, @SVector([x.max+Δx/2.0; t-1*Δt]))
            UE0   = sol.u
            # Current step
            sol   = Analytics(problem, @SVector([x.max+Δx/2.0; t-0*Δt]))
            UE    = sol.u
        else # Inside
            UE00 = U00[i+1]
            UE0  = U0[i+1]
            UE   = U[i+1]
        end
    
        # Balance equation

        # # Conventional - style 1:
        # F[i] = ρ*(1*UC -2*UC0 + 1*UC00)/Δt^2 - E*(1*UW0  -2*UC0 + 1*UE0)/Δx^2

        # # Conventional - style 2:
        # Mt = @SVector([0.;   1.;   0.;   0.;  -2.; 0.;  0.; 1.; 0.]) 
        # Mx = @SVector([0.;   0.;   0.;   1.;  -2.; 1.;  0.; 0.; 0.])
        # u  = @SVector([UW00; UC00; UE00; UW0; UC0; UE0; UW; UC; UE])
        # F[i] = ρ/Δt^2*Mt'*u  - E/Δx^2*Mx'u

        # Conventional - style 2:
        Mt = @SVector([1/12;     10/12;          1/12;    -2/12;    -20/12;           -2/12;     1/12;    10/12;           1/12]) 
        Mx = @SVector([1/12*GW; -2/12*(GW+GE)/2; 1/12*GE; 10/12*GW; -20/12*(GW+GE)/2; 10/12*GE;  1/12*GW; -2/12*(GW+GE)/2; 1/12*GE])
        u  = @SVector([UW00; UC00; UE00; UW0; UC0; UE0; UW; UC; UE])
        F[i] = ρ[i]/Δt^2*Mt'*u  - 1/Δx^2*Mx'u - f[i]
    end
    return nothing
end

Residual!(F, U, U0, U00, f, G, ρ, Δx, Δt, x, t)  = Residual_OSOA!(F, U, U0, U00, f, G, ρ, Δx, Δt, x, t)

function main_Wave1D_OSOA(Δx, Δt, ncx, nt, L)

    # Parameters
    x   = (min=-L/2, max=L/2)
    xv  = LinRange(x.min,       x.max,       ncx+1)
    xc  = LinRange(x.min+Δx/2., x.max-Δx/2., ncx  )
    t   = 0.

    # Allocations
    U   = zeros(ncx)
    U0  = zeros(ncx)
    U00 = zeros(ncx)
    f   = zeros(ncx)
    δU  = zeros(ncx)
    Ua  = zeros(ncx)
    F   = zeros(ncx)
    ρ   = zeros(ncx)
    G   = zeros(ncx+1)
    Gc  = zeros(ncx)

    for i in eachindex(G)
        sol   = Analytics(problem, @SVector([xv[i]; t]))
        G[i]  = sol.G
    end

    for i in eachindex(ρ)
        sol   = Analytics(problem, @SVector([xc[i]; t]))
        Gc[i] = sol.G
        ρ[i]  = sol.ρ
    end

    # If we assume tat coefficient is defined on the U points, then one has to introduce averaging
    G[2:end-1] .= 0.5*(Gc[1:end-1] .+ Gc[2:end])

    # Sparsity pattern
    input       = rand(ncx)
    output      = similar(input)
    Res_closed! = (F, U) -> Residual!(F, U, U0, U00, f, G, ρ, Δx, Δt, x, t)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))

    # Makes coloring
    colors      = matrix_colors(J)

    # Initial condition: Evaluate exact initial solution
    for i in eachindex(U)
        sol   = Analytics(problem, @SVector([xc[i]; t   ]))
        U[i]  = sol.u
        sol   = Analytics(problem, @SVector([xc[i]; t-Δt]))
        U0[i] = sol.u
    end
    
    # Time loop
    for it=1:nt
        U00 .= U0
        U0  .= U
        t   += Δt 
        noisy==true ? @printf("########### Step %06d ###########\n", it) : nothing

        for i in eachindex(f)
            sol   = Analytics(problem, @SVector([xc[i]; t-Δt])) # ???
            f[i]  = sol.s
        end

        r1 = 1.0
        for iter=1:10

            # Residual evaluation: T is found if F = 0
            Residual!(F, U, U0, U00, f, G, ρ, Δx, Δt, x, t)
            Res_closed! = (F, U) -> Residual!(F, U, U0, U00, f, G, ρ, Δx, Δt, x, t)
            r = norm(F)/ncx
            if iter==1 r1 = r; end
            noisy==true ? @printf("## Iteration %06d: r/r1 = %1.2e ##\n", iter, r/r1) : nothing
            if r/r1 < 1e-8 break end
                
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
        sol   = Analytics(problem, @SVector([xc[i]; t]))
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

    ρ   = 1.0 # Indicative values
    G   = 1.0 # Indicative values

    # Time
    ncx = 400
    Δx  = L/ncx
    Nt  = [80, 160, 320, 640]    
    Δtv = 0.01 ./ Nt
    ϵt  = zero(Δtv)
    for i in eachindex(Nt)
        Δx    = 2.1 * (Δtv[i]*sqrt(G/ρ)) * 100 # Somehow dx has to be large for proper dt
        ncx   = Int64(floor(L/Δx))
        L     = ncx*Δx
        ϵt[i] = main_Wave1D_OSOA(Δx, Δtv[i], ncx, Nt[i], L)
    end

    # Space
    L   = 2.
    nt  = 100
    Ncx = [160, 320, 640]  
    Δxv = 2.0 ./ Ncx
    ϵx  = zero(Δtv)
    for i in eachindex(Ncx)
        Δt    = Δxv[i]/sqrt(G/ρ)/2.1 / 20
        nt    = Int64(floor(tt/Δt))
        ϵx[i] = main_Wave1D_OSOA(Δxv[i], Δt, Ncx[i], nt, L)
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