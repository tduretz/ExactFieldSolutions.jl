using ExactFieldSolutions, Test

@testset "Schmid2003" verbose=true begin
    sol = Stokes2D_Schmid2003( (0, 0) )
    L = (xx = -0.019801980198019802, xy = 0.0, yx = 0.0, yy = 0.019801980198019802)
    @test sol.L.xx - L.xx ≈ 0.0
    @test sol.L.xy - L.xy ≈ 0.0
    @test sol.L.yx - L.yx ≈ 0.0
    @test sol.L.yy - L.yy ≈ 0.0
end