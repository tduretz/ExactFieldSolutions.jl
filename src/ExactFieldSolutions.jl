module ExactFieldSolutions

using Plots, SpecialFunctions, ForwardDiff, UnPack, StaticArrays, Printf

include("Wave1D/Wave1D_dAlembert.jl")
export Wave1D_dAlembert

include("Wave1D/Wave1D_HeteroPlusSource.jl")
export Wave1D_HeteroPlusSource

include("Diffusion1D/Diffusion1D_Gaussian.jl")
export Diffusion1D_Gaussian

include("Diffusion2D/Diffusion2D_Gaussian.jl")
export Diffusion2D_Gaussian

include("Poisson2D/Poisson2D_Sevilla2018.jl")
# export Poisson2D_Sevilla2018_enz
export Poisson2D_Sevilla2018

include("Poisson3D/Poisson3D_Sevilla2018.jl")
export Poisson3D_Sevilla2018

include("Stokes2D/Stokes2D_Donea2003.jl")
export Stokes2D_Donea2003

include("Stokes2D/Stokes2D_Schmid2003.jl")
export Stokes2D_Schmid2003

include("Stokes2D/Stokes2D_Moulas2021.jl")
export Stokes2D_Moulas2021

include("Stokes2D/Stokes2D_SolKz_Zhong1996.jl")
export Stokes2D_SolKz_Zhong1996

include("Stokes2D/Stokes2D_SolCx_Zhong1996.jl")
export Stokes2D_SolCx_Zhong1996

include("Elasticity2D/Elasticity2D_Hole.jl")
export Elasticity2D_Hole

end # module ExactFieldSolutions
