## ExactFieldSolutions.jl

Full field solutions are essential for the verification of numerical codes that are based on the solution of Partial Differential Equations (PDE).
They allow to check if numerical solutions are meaningful in eyeball norm but, more importantly, they allow to quantify discretisation errors.
This quantification further enables to check whether numerical solutions converge to exact solutions at expected theoretical rates.
`ExactFieldSolutions` compiles full field solutions for 2D and 3D PDE problems including Poisson-type and mechanical problems (Stokes, elasticity).
Contributions are welcome and full field solutions to other problems (electric, magnetic, MHD) are more than welcome. Feel free to make a PR.

`ExactFieldSolutions` benefits from automatic differentiation tools available within the Julia ecosystem (e.g., ForwardDiff). These allow to evaluate fluxes and sources terms in a simplified way. [See this manufactured solution for the 2D Poisson problem](src/Poisson2D_Sevilla2018.jl).

Please note that `ExactFieldSolutions` is a registered package, so you can install it simply by typing `add ExactFieldSolutions` in package mode.

### Poisson 2D
[Sevilla et al. (2018)](examples/Visualize_Poisson2D_Sevilla2018.jl)

![alt text](img/Poisson2D_Sevilla2018.svg "Sevilla et al. (2018)")

### Poisson 3D
[Sevilla et al. (2018)](examples/Visualize_Poisson3D_Sevilla2018.jl)

![alt text](img/Poisson3D_Sevilla2018.svg "Sevilla et al. (2018)")

### Stokes 2D
[Viscous inclusion - Schmid & Podladchikov (2003)](examples/Visualize_Stokes2D_Schmid2003.jl)

![alt text](img/Stokes2D_Schmid2003.svg "Schmid & Podladchikov (2003)")

[Double corner flow - Moulas et al., (2021)](examples/Visualize_Stokes2D_Moulas2021.jl)

![alt text](img/Stokes2D_Moulas2021.svg "Moulas et al. (2021)")

[Donea & Huerta (2003)](examples/Visualize_Stokes2D_Donea2003.jl)

![alt text](img/Stokes2D_Donea2003.svg "Donea & Huerta (2003)")

[SolKz - Zhong et al. (1996)](examples/Visualize_Stokes2D_SolKz_Zhong1996.jl)

![alt text](img/Stokes2D_SolKz_Zhong1996.svg "Zhong et al. (1996)")


[SolCx - Zhong et al. (1996)](examples/Visualize_Stokes2D_SolCx_Zhong1996.jl)

![alt text](img/Stokes2D_SolCx_Zhong1996.svg "Zhong et al. (1996)")

### Elasticity 2D
[Elastic plate with a hole](examples/Visualize_Elasticity2D_Hole.jl)

![alt text](img/Elasticity2D_Hole.svg "Elastic plate with a hole")
