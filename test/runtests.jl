using ExactFieldSolutions, Test

@testset "Schmid2003" verbose=true begin
    p   = 0.0
    V   = (x = -0.15099009900990099, y = 0.15099009900990104)
    L   = (xx = -0.7549504950495054, xy = -0.9801980198019802, yx = 0.9801980198019797, yy = 0.7549504950495045)
    ε̇   = (xx = -0.7549504950495054, xy = -2.220446049250313e-16, yx = -2.220446049250313e-16, yy = 0.7549504950495045)
    τ   = (xx = -1.5099009900990108, xy = -4.440892098500626e-16, yx = -4.440892098500626e-16, yy = 1.509900990099009)
    η   = 1.0   
    sol = Stokes2D_Schmid2003_circle( (0.2, 0.2) )
    @test sol.η              ≈ η 
    @test sol.p              ≈ p
    @test all(values(sol.V) .≈ values(V))
    @test all(values(sol.L) .≈ values(L))
    @test all(values(sol.ε̇) .≈ values(ε̇))
    @test all(values(sol.τ) .≈ values(τ))
end


@testset "Moutzouris2026" verbose=true begin
    p    = -0.05393579072532709
    V    = [0.2006760424781664, -0.5364073967772356]
    τ    = [2.0133874937966665  0.09406162472906811;
            0.09406162472906811 -1.9774302999797817]
    τzz  = -0.035957193816884725
    sol = Stokes2D_Moutzouris_circle( (0.2, 0.5) )
    @test sol.p    ≈ p
    @test all(sol.V .≈ V)
    @test all(sol.τ .≈ τ)
    @test sol.τzz  ≈ τzz
end

@testset "far-field" begin
    for prm in ((ηm=1.0, ηi=0.1,  ξm=1.0, ξi=10.0, ri=0.1, t=2.0, α=0.0, ε̇=1.0,  γ̇=0.0, ζ̇=0.0,  ε̇zz=0.0),
                (ηm=1.0, ηi=10.0, ξm=1.0, ξi=0.5,  ri=0.1, t=2.0, α=0.0, ε̇=0.3,  γ̇=1.5, ζ̇=0.2,  ε̇zz=0.0),
                (ηm=2.0, ηi=0.05, ξm=3.0, ξi=1.0,  ri=0.1, t=2.0, α=0.0, ε̇=-0.5, γ̇=0.7, ζ̇=-0.1, ε̇zz=0.4))
        L∞ = [prm.ζ̇+prm.ε̇  prm.γ̇; 0.0  prm.ζ̇-prm.ε̇]
        for x in ([100.0, 0.0], [0.0, 100.0], [70.0, -70.0])
            sol = Stokes2D_Moutzouris_circle(x; params=prm)
            @test sol.L ≈ L∞          rtol=1e-5
            @test sol.V ≈ L∞*x        rtol=1e-5
        end
    end
end

@testset "incompressible limit" begin
    prm = (ηm=1.0, ηi=100.0, ξm=1e8, ξi=1e8, ri=0.2, t=2.0, α=0.0, ε̇=-1.0, γ̇=0.0, ζ̇=0.0, ε̇zz=0.0)
    for x in ([0.1,0.05], [0.5,0.3], [1.0,-0.7])
        m = Stokes2D_Moutzouris_ellipse(x; params=prm)
        s = Stokes2D_Schmid2003_ellipse(x; params=prm)
        @test m.V ≈ s.V  rtol=1e-6
        @test m.p ≈ s.p  rtol=1e-6
        @test m.τ ≈ s.τ  rtol=1e-6
    end
end