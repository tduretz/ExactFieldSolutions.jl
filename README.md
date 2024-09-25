## ExactFieldSolutions.jl

<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliageodynamics.github.io/GeoParams.jl/dev/) -->
[![CI](https://github.com/tduretz/ExactFieldSolutions.jl//actions/workflows/blank.yml/badge.svg)](https://github.com/tduretz/ExactFieldSolutions.jl//actions/workflows/blank.yml)
<!-- [![DOI](https://zenodo.org/badge/369433137.svg)](https://zenodo.org/doi/10.5281/zenodo.8089230) -->
[![codecov](https://codecov.io/gh/tduretz/ExactFieldSolutions.jl)](https://codecov.io/gh/tduretz/ExactFieldSolutions.jl)

Full field solutions are essential for the verification of numerical codes that are based on the solution of Partial Differential Equations (PDE).
They allow for checking if numerical solutions are meaningful in eyeball norm but, more importantly, they allow for the quantification of discretisation errors.
This quantification further enables to check whether numerical solutions converge to exact solutions at expected theoretical rates.
`ExactFieldSolutions` compiles full field solutions for 1D, 2D and 3D PDE problems including Poisson-type and mechanical problems (Stokes, elasticity).
Contributions are welcome and full field solutions to other problems (electric, magnetic, MHD) are more than welcome. Feel free to make a PR.

`ExactFieldSolutions` benefits from automatic differentiation tools available within the Julia ecosystem (e.g., ForwardDiff). These allow to evaluate fluxes and sources terms in a simplified way. [See this manufactured solution for the 2D Poisson problem](src/Poisson2D_Sevilla2018.jl).

Please note that `ExactFieldSolutions` is a registered package, so you can install it simply by typing `add ExactFieldSolutions` in package mode.

## Visualisation examples

It is necessary to activate the example environment in order to reproduce the visualisations and benchmarks, one can use the package mode for this purpose (`]`):

```julia-repl
(ExactFieldSolutions) pkg> activate examples/
  Activating project at `~/REPO/ExactFieldSolutions/examples`

(examples) pkg>
```

### Poisson 2D
[Sevilla et al. (2018)](ext/visualisations/Visualize_Poisson2D_Sevilla2018.jl)

![alt text](img/Poisson2D_Sevilla2018.svg "Sevilla et al. (2018)")

### Poisson 3D
[Sevilla et al. (2018)](ext/visualisations/Visualize_Poisson3D_Sevilla2018.jl)

![alt text](img/Poisson3D_Sevilla2018.svg "Sevilla et al. (2018)")

### Diffusion 1D
[Diffusion of a 1D Gaussian](ext/visualisations/Visualize_Diffusion1D_Gaussian.jl)

![alt text](img/Diffusion1D_Gaussian.svg)

### Diffusion 2D
[Diffusion of a 2D Gaussian](ext/visualisations/Visualize_Diffusion2D_Gaussian.jl)

![alt text](img/Diffusion2D_Gaussian.svg)

### Wave 1D 
[Propagation of a 1D wave based on d'Alembert's solution](ext/visualisations/Visualize_Wave1D_dAlembert.jl)

![alt text](img/Wave1D_dAlembert.svg)

[Propagation of a 1D wave with a time-space dependent source and variable coefficient](ext/visualisations/Visualize_Wave1D_HeteroPlusSource.jl)

![alt text](img/Wave1D_HeteroPlusSource.svg)

### Stokes 2D
[Viscous inclusion - Schmid & Podladchikov (2003)](ext/visualisations/Visualize_Stokes2D_Schmid2003.jl)

![alt text](img/Stokes2D_Schmid2003.svg "Schmid & Podladchikov (2003)")

[Double corner flow - Moulas et al., (2021)](ext/visualisations/Visualize_Stokes2D_Moulas2021.jl)

![alt text](img/Stokes2D_Moulas2021.svg "Moulas et al. (2021)")

[Donea & Huerta (2003)](ext/visualisations/Visualize_Stokes2D_Donea2003.jl)

![alt text](img/Stokes2D_Donea2003.svg "Donea & Huerta (2003)")

[SolKz - Zhong et al. (1996)](ext/visualisations/Visualize_Stokes2D_SolKz_Zhong1996.jl)

![alt text](img/Stokes2D_SolKz_Zhong1996.svg "Zhong et al. (1996)")


[SolCx - Zhong et al. (1996)](ext/visualisations/Visualize_Stokes2D_SolCx_Zhong1996.jl)

![alt text](img/Stokes2D_SolCx_Zhong1996.svg "Zhong et al. (1996)")

### Elasticity 2D
[Elastic plate with a hole](ext/visualisations/Visualize_Elasticity2D_Hole.jl)

![alt text](img/Elasticity2D_Hole.svg "Elastic plate with a hole")

## Benchmarking examples

### Diffusion 1D

[1D diffusion: Finite Difference Method (FDM) with backward-Euler integration and spatial staggering (θ = 1.0)](ext/benchmarks/Benchmark_Diffusion1D.jl)

![alt text](img/Benchmark_Diffusion1D_BackwardEuler_FDM.svg "Diffusion in 1D using the Finite Difference Method (FDM): backward-Euler and spatial staggering") 

[1D diffusion: Finite Difference Method (FDM) with Crank-Nicolson integration and spatial staggering (θ = 0.5)](ext/benchmarks/Benchmark_Diffusion1D.jl)

![alt text](img/Benchmark_Diffusion1D_CrankNicolson_FDM.svg "Diffusion in 1D using the Finite Difference Method (FDM): Crank-Nicolson and spatial staggering")

### Poisson 2D

[2D Poisson: Finite Difference Method (FDM) using an O(4) compact discretisation with constant coefficient](ext/benchmarks/Benchmark_Poisson2D.jl)

![alt text](img/Benchmark_Poisson2D_O4_FDM.svg "2D Poisson problem using O(4) compact finite difference scheme")

### Wave 1D

[1D wave: Finite Difference Method (FDM) with velocity-stress discretisation](ext/benchmarks/Benchmark_Wave1D_VelStress_FDM.jl)

![alt text](img/Benchmark_Wave1D_VelStress_FDM.svg "Wave in 1D using the Finite Difference Method (FDM): velocity-stress scheme") 

[1D wave: Finite Difference Method (FDM) with conventional O(2) discretisation](ext/benchmarks/Benchmark_Wave1D_Conventional_FDM.jl)

![alt text](img/Benchmark_Wave1D_Conventional_FDM.svg "Wave in 1D using the Finite Difference Method (FDM): conventional O(2) discretisation") 

[1D wave: Finite Difference Method (FDM) with optimised O(4) discretisation with constant coefficient](ext/benchmarks/Benchmark_Wave1D_OptimallyAccurate_FDM.jl)

![alt text](img/Benchmark_Wave1D_OptimallyAccurate_FDM.svg "Wave in 1D using the Finite Difference Method (FDM): optimised O(4) discretisation with constant coefficient") 
