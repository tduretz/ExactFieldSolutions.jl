using ExactFieldSolutions, Test

@testset "Schmid2003" verbose=true begin
    p = 0.0
    V = (x = 0.0, y = 0.0)
    L = (xx = -0.019801980198019802, xy = 0.0, yx = 0.0, yy = 0.019801980198019802)
    ε̇ = (xx = -0.019801980198019802, xy = 0.0, yx = 0.0, yy = 0.019801980198019802)
    τ = (xx = -3.9603960396039604, xy = 0.0, yx = 0.0, yy = 3.9603960396039604)
    η = 100   
    sol = Stokes2D_Schmid2003( (0, 0) )
    @test sol.η    - η    ≈ 0.0
    @test sol.p    - p    ≈ 0.0
    @test sol.V.x  - V.x  ≈ 0.0
    @test sol.V.y  - V.y  ≈ 0.0
    @test sol.L.xx - L.xx ≈ 0.0
    @test sol.L.xy - L.xy ≈ 0.0
    @test sol.L.yx - L.yx ≈ 0.0
    @test sol.L.yy - L.yy ≈ 0.0
    @test sol.ε̇.xx - ε̇.xx ≈ 0.0
    @test sol.ε̇.xy - ε̇.xy ≈ 0.0
    @test sol.ε̇.yx - ε̇.yx ≈ 0.0
    @test sol.ε̇.yy - ε̇.yy ≈ 0.0
    @test sol.τ.xx - τ.xx ≈ 0.0
    @test sol.τ.xy - τ.xy ≈ 0.0
    @test sol.τ.yx - τ.yx ≈ 0.0
    @test sol.τ.yy - τ.yy ≈ 0.0
end