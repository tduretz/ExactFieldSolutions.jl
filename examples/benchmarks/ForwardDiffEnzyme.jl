# Hessian, jacobian and gradient with parameters - 2D Poisson problem
using UnPack, StaticArrays
import LinearAlgebra: tr, norm
import ForwardDiff, Enzyme 

function ScalarFunction(x, params)
    @unpack α, β, a, b, c, d = params
    return exp(α*sin(a*x[1] + c*x[2]) + β*cos(b*x[1] + d*x[2]))
end

function VectorFunction(x, params)
    @unpack α, β, a, b, c, d = params
    return [exp(α*sin(a*x[1] + c*x[2])),  exp(β*cos(b*x[1] + d*x[2]))]
end

function Poisson2D_Sevilla2018(x::SVector{2,Float64};
    params = (α = 0.1, β = 0.3, a = 5.1, b = 4.3, c = -6.2, d = 3.4) )
    
    # Doing things with ForwardDiff
    u      = ScalarFunction(x, params)
    f_cl   = x -> ScalarFunction(x, params)
    g_cl   = x -> VectorFunction(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    jacu   = ForwardDiff.jacobian(g_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = -tr(hessu)

    # Check source term
    @unpack α, β, a, b, c, d = params
    s_check = u*(-a*α*cos(a*x[1] + c*x[2]) + b*β*sin(b*x[1] + d*x[2]))*(a*α*cos(a*x[1] + c*x[2]) - b*β*sin(b*x[1] + d*x[2])) + u*(a^2*α*sin(a*x[1] + c*x[2]) + b^2*β*cos(b*x[1] + d*x[2])) + u*(-α*c*cos(a*x[1] + c*x[2]) + β*d*sin(b*x[1] + d*x[2]))*(α*c*cos(a*x[1] + c*x[2]) - β*d*sin(b*x[1] + d*x[2])) + u*(α*c^2*sin(a*x[1] + c*x[2]) + β*d^2*cos(b*x[1] + d*x[2]))
    @assert abs(s-s_check) ≤ 1e-13

    # Doing the same with Enzyme:
    # 1) Using gradient works fine! 
    gradu_enz = Enzyme.gradient(Enzyme.Reverse, f_cl, x)
    @assert norm(gradu_enz[1] - gradu) ≤ 1e-13

    # 2) Using jacobian returns an error
    jacu_enz = Enzyme.jacobian(Enzyme.Reverse, g_cl, x)
    @assert norm(jacu_enz[1] - jacu) ≤ 1e-13

    # 3) Is there a Hessian API?

    return (u=u, ∇u=gradu, s=s)
end

Poisson2D_Sevilla2018(@SVector([1.0; 2.0]))
