module ExactFieldSolutions 

using SpecialFunctions, ForwardDiff, UnPack, StaticArrays, Printf, FastPow

include("Wave1D/Wave1D_dAlembert.jl")
export Wave1D_dAlembert

include("Wave1D/Wave1D_HeteroPlusSource.jl")
export Wave1D_HeteroPlusSource

include("Wave1D/Wave1D_Source.jl")
export Wave1D_Source

include("Diffusion1D/Diffusion1D_Gaussian.jl") 
export Diffusion1D_Gaussian

include("Diffusion1D/Diffusion1D_StefanProblem.jl") 
export Diffusion1D_StefanProblem

include("Diffusion2D/Diffusion2D_Gaussian.jl")
export Diffusion2D_Gaussian

include("Poisson1D/Poisson1D_VarCoeff.jl")
export Poisson1D_VarCoeff

include("Poisson2D/Poisson2D_Sevilla2018.jl")
export Poisson2D_Sevilla2018

include("Poisson2D/Poisson2D_VarCoeff.jl")
export Poisson2D_VarCoeff

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

include("Poroelasticity2D/Poroelasticity2D_Rozhko2008.jl")
export Poroelasticity2D_Rozhko2008

end # module ExactFieldSolutions
