using Plots, ExactFieldSolutions
import Statistics: mean

function main_VelStress(Δx, Δt, ncx, nt, L)

    # Parameters
    x   = (min=-L/2, max=L/2)
    xc  = LinRange(x.min-Δx/2., x.max+Δx/2., ncx+2) # with ghosts nodes
    ρ   = 1.0
    E   = 1.0
    t   = 0.
    k   = 8.
    c   = sqrt(E/ρ)

    # Allocations
    V   = zeros(ncx+2)
    Va  = zeros(ncx+2)
    σ   = zeros(ncx+1)

    # Initial condition: Evaluate exact initial solution
    params = (c=c, k=k)
    for i in eachindex(V)
        sol   = Wave1D_dAlembert([xc[i]; t]; params)
        V[i] = sol.∇u[2]
    end

    for i in eachindex(σ)
        sol   = Wave1D_dAlembert([xc[i]; t]; params)
        σ[i] = E*sol.∇u[1]
    end

    # Time loop
    for it=1:nt
        t += Δt 
        # Set BCs
        sol          = Wave1D_dAlembert([xc[1];   t]; params)
        V[1]         = sol.∇u[2]
        sol          = Wave1D_dAlembert([xc[end]; t]; params)
        V[end]       = sol.∇u[2]
        σ          .+= Δt*E*diff(V, dims=1)/Δx
        V[2:end-1] .+= Δt/ρ*diff(σ, dims=1)/Δx 
    end

    # Evaluate exact solution
    for i in eachindex(V)
        sol   = Wave1D_dAlembert([xc[i]; t]; params)
        Va[i] = sol.∇u[2]
    end

    # p = plot()
    # p = plot!(xc, V,  label = "u" )
    # p = plot!(xc, Ua, label = "ua")
    # display(p)
    # sleep(0.1)

    # Error
    return mean(abs.(V .- Va))
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
        Δx    = 2.1 * (Δtv[i]*sqrt(E/ρ))
        ncx   = Int64(floor(L/Δx))
        L     = ncx*Δx
        ϵt[i] = main_VelStress(Δx, Δtv[i], ncx, Nt[i], L)
    end

     # Time
     L   = 2.
     nt  = 100
     Ncx = [20, 40, 80, 160]  
     Δxv = 2.0 ./ Ncx
     ϵx  = zero(Δtv)
     for i in eachindex(Ncx)
        Δt    = Δxv[i]/sqrt(E/ρ)/2.1
        nt    = Int64(floor(tt/Δt))
        ϵx[i] = main_VelStress(Δxv[i], Δt, Ncx[i], nt, L)
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

# begin

    # ρ   = 1.0
    # E   = 1.0
    # t   = 0.
    # k   = 8.
    # c   = sqrt(E/ρ)

#     L   = 1.0
#     t   = 0.2
#     ncx = 500
#     Δx  = L/ncx
#     nt  = 1000
#     Δt  = Δx/sqrt(E/ρ)/2.1
#     nt  = Int64(floor(t/Δt))
#     main(Δx, Δt, ncx, nt, L)
# end