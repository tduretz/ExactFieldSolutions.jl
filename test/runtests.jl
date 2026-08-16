using ExactFieldSolutions, Test

@testset "Schmid2003" verbose=true begin
    p   = 0.0
    V   = (x = 0.0, y = 0.0)
    L   = (xx = -0.019801980198019802, xy = 0.0, yx = 0.0, yy = 0.019801980198019802)
    ε̇   = (xx = -0.019801980198019802, xy = 0.0, yx = 0.0, yy = 0.019801980198019802)
    τ   = (xx = -3.9603960396039604, xy = 0.0, yx = 0.0, yy = 3.9603960396039604)
    η   = 100   
    sol = Stokes2D_Schmid2003( (0, 0) )
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