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


@testset "Duretz2026" verbose=true begin
    p   = -0.09568944517573587
    V   = (x = -0.2089090434970527, y = 0.47592372115980974)
    τ = (xx = -2.00581033036802, xy = 0.12557580379582292, yx = 0.12557580379582292, yy = 2.0064482600025246, zz = -0.0006379296345049024)
    sol = Stokes2D_Duretz2026( (0.2, 0.5) )
    @test sol.p              ≈ p
    @test all(values(sol.V) .≈ values(V))
    @test all(values(sol.τ) .≈ values(τ))
end