function Wave1D_dAlembert_u_fwd(x, params)
    @unpack c, k  = params
    return real(exp(im*k*(x[1] - c*x[2])))
end

function Wave1D_dAlembert(x;
    params = (c=1.0, k=8.0) )
    u      = Wave1D_dAlembert_u_fwd(x, params)
    f_cl   = x -> Wave1D_dAlembert_u_fwd(x, params)
    gradu  = ForwardDiff.gradient(f_cl, x)
    hessu  = ForwardDiff.hessian(f_cl, x)
    s      = hessu[2,2] - params.c^2*(hessu[1,1])
    return (u=u, ∇u=gradu, s=s)
end