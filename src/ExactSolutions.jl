module ExactSolutions

using Plots, SpecialFunctions, ForwardDiff, UnPack

include("Poisson2D/Poisson2D_Sevilla2018.jl")
# export Poisson2D_Sevilla2018_enz
export Poisson2D_Sevilla2018

include("Poisson3D/Poisson3D_Sevilla2018.jl")
export Poisson3D_Sevilla2018

include("Stokes2D/Stokes2D_Schmid2003.jl")
export Stokes2D_Schmid2003

include("Stokes2D/Stokes2D_SolKz_Zhong1996.jl")
export Stokes2D_SolKz_Zhong1996

include("Stokes2D/Stokes2D_SolCx_Zhong1996.jl")
export Stokes2D_SolCx_Zhong1996

include("Elasticity2D/Elasticity2D_Hole.jl")
export Elasticity2D_Hole

end # module ExactSolutions
