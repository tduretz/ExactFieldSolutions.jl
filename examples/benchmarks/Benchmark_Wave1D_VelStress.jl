using Plots, ExactFieldSolutions, StaticArrays
import Statistics: mean

# Dynamic dispatch based on a rule
_type(rule) = Val{rule}()
Analytics(rule, x)                 = _func(_type(rule), x) 
_func(::Val{:dAlembert}       , x) = Wave1D_dAlembert(x)
_func(::Val{:HeteroPlusSource}, x) = Wave1D_HeteroPlusSource(x)

# Problem type
problem = :dAlembert
problem = :HeteroPlusSource

function main_VelStress(Δx, Δt, ncx, nt, L)

    # Parameters
    x   = (min=-L/2, max=L/2)
    xv  = LinRange(x.min,       x.max,       ncx+1)
    xc  = LinRange(x.min-Δx/2., x.max+Δx/2., ncx+2) # with ghosts nodes
    t   = 0.0

    # Allocations
    V   = zeros(ncx+2)
    Va  = zeros(ncx+2)
    σ   = zeros(ncx+1)
    G   = zeros(ncx+1)
    f   = zeros(ncx+2)
    ρ   = zeros(ncx+2)

    # Initial condition: Evaluate exact initial solution    
    for i in eachindex(V)
        sol  = Analytics(problem, [xc[i]; t])
        V[i] = sol.∇u[2]
        ρ[i] = sol.ρ
    end

    for i in eachindex(σ)
        sol  = Analytics(problem,[xv[i]; t-1*Δt/2]) # ⚠ time staggering
        G[i] = sol.G
        σ[i] = G[i]*sol.∇u[1]
    end

    # Time loop
    for it=1:nt

        # Update force
        for i in eachindex(V)
            sol  = Analytics(problem, @SVector([xc[i]; t]))
            f[i] = sol.s
        end

        # Set BCs
        sol          = Analytics(problem,[xc[1];   t])
        V[1]         = sol.∇u[2]
        sol          = Analytics(problem,[xc[end]; t])
        V[end]       = sol.∇u[2]
        σ          .+= Δt *G         .* diff(V, dims=1)/Δx
        V[2:end-1] .+= Δt./ρ[2:end-1].*(diff(σ, dims=1)/Δx + f[2:end-1]) 
        # Time update
        t += Δt 
    end

    # Evaluate exact solution
    for i in eachindex(V)
        sol   = Analytics(problem,[xc[i]; t])
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

    ρ    = 2.0
    maxG = 10.0

    # Time
    ncx = 400
    Δx  = L/ncx
    Nt  = [160, 320, 640]  
    Δtv = 0.01 ./ Nt
    ϵt  = zero(Δtv)
    for i in eachindex(Nt)
        Δx    = 2.1 * (Δtv[i]*sqrt(maxG/ρ)) * 80
        ncx   = Int64(floor(L/Δx))
        L     = ncx*Δx
        ϵt[i] = main_VelStress(Δx, Δtv[i], ncx, Nt[i], L)
    end

     # Space
     L   = 2.
     nt  = 100
     Ncx = [160, 320, 640]  
     Δxv = 2.0 ./ Ncx
     ϵx  = zero(Δtv)
     for i in eachindex(Ncx)
        Δt    = Δxv[i]/sqrt(maxG/ρ)/2.1 /20
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