@doc raw"""
    sol = Stokes2D_SolKz_Zhong1996(x; params)  

    Evaluates the manufactured solution of [Zhong et al., (1996)](https://booksite.elsevier.com/brochures/geophysics/PDFs/00118.pdf). 
    Analytical solution adapted from Underworld (https://github.com/underworldcode, LGPL-3 license).
    The solution requires the following expressions the density and viscosity fields: 
    ρ = -σ*sin(km*x[2])*cos(kn*x[1]) 
    η = exp(2*B*x[2])
    
    Input parameters are:
    x      : is the coordinate vector or tuple
    params : optional parameter array, default = (Δη=1e6, B = 6.9, km = 1.6*π, n = 3, σ = 1.0)
and returns:

    sol    : tuple containing the solution fields p (pressure), V (velocity vector), ρ (density) and η (viscosity)

# Examples
```julia-repl
julia> Stokes2D_SolKz_Zhong1996([0; 0])
(p = -0.017494712727797308, V = [0.0, 0.0], ρ = 0.0, η = 1.0)
```

```julia-repl
julia> Stokes2D_SolKz_Zhong1996((0,0))
(p = -0.017494712727797308, V = (x = 0.0, y = 0.0), η = 1.0, ρ = 0.0)
```
"""
function Stokes2D_SolKz_Zhong1996(x; 
    params = (Δη=1e6, B = 6.9, km = 1.6*π, n = 3, σ = 1.0  ) )
    @unpack Δη, B, km, n, σ = params
	#************************************************************************%
	#************************************************************************%
	kn = 3*π
	a  = B*B + kn*kn
	b  = 2.0*kn*B
	r  = sqrt(a*a + b*b)
	Rp = sqrt( (r+a)/2.0 )
	Rm = sqrt( (r-a)/2.0 )
	UU = Rp - B
	VV = Rp + B
	
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	sum6 = 0.0;
	sum7 = 0.0;
	
	#*******************************************%
	#          calculate the constants          %
	#*******************************************%
	
	t3 = kn * kn;
	t4 = km * km;
	t6 = B * B;
	t8 = 0.4e1 * t3 * t6;
	t10 = 0.4e1 * t4 * t6;
	t13 = 0.8e1 * kn * t6 * km;
	t14 = t4 * t4;
	t16 = 0.2e1 * t3 * t4;
	t17 = t3 * t3;
	aa = -0.4e1 * B * km * kn * (t3 + t4) / (t8 + t10 + t13 + t14 + t16 + t17) / (-t13 + t8 + t10 + t14 + t16 + t17);
	
	t1 = kn * kn;
	t2 = t1 * t1;
	t3 = B * B;
	t5 = 0.4e1 * t1 * t3;
	t6 = km * km;
	t7 = t6 * t6;
	t9 = 0.2e1 * t1 * t6;
	t11 = 0.4e1 * t3 * t6;
	t16 = 0.8e1 * kn * t3 * km;
	bb = kn * (t2 + t5 + t7 + t9 - t11) / (t5 + t11 + t16 + t7 + t9 + t2) / (-t16 + t5 + t11 + t7 + t9 + t2);
	
	AA = aa;
	BB = bb;
	
	t1 = B * B;
	t2 = t1 * Rp;
	t4 = Rm * Rm;
	t5 = t4 * Rp;
	t7 = t4 * B;
	t8 = km * km;
	t12 = Rp * Rp;
	t13 = B * t12;
	t21 = 0.8e1 * t1 * km * BB * Rp;
	t23 = 0.2e1 * Rm;
	t24 = cos(t23);
	t26 = Rm * Rp;
	t38 = sin(t23);
	t51 = exp(-0.2e1 * Rp);
	t53 = B + Rp;
	t54 = Rm * t53;
	t55 = Rm * B;
	t57 = 0.2e1 * B * km;
	t58 = t55 + t57 - t26;
	t62 = 0.3e1 * t1;
	t64 = 0.2e1 * Rp * B;
	t65 = t62 + t64 + t4 - t8 - t12;
	t67 = t54 * t65 * BB;
	t69 = Rm - km;
	t70 = cos(t69);
	t72 = -t57 + t55 - t26;
	t77 = Rm + km;
	t78 = cos(t77);
	t81 = t54 * t65 * AA;
	t86 = sin(t77);
	t92 = sin(t69);
	t96 = exp(-t53);
	t98 = B - Rp;
	t99 = Rm * t98;
	t100 = t55 + t57 + t26;
	t104 = t62 - t64 + t4 - t8 - t12;
	t106 = t99 * t104 * BB;
	t109 = -t57 + t55 + t26;
	t116 = t99 * t104 * AA;
	t130 = exp(-0.3e1 * Rp - B);
	t135 = exp(-0.4e1 * Rp);
	t144 = t4 * t1;
	t150 = t4 * t12;
	C1 = (((0.2e1 * Rp * (0.2e1 * t2 + 0.2e1 * t5 + t7 + B * t8 - 0.3e1 * t1 * B + t13) * AA + t21) * t24 + (-0.2e1 * t26 * (t4 - t8 - t12 + 0.5e1 * t1) * AA + 0.8e1 * B * BB * km * Rm * Rp) * t38 - 0.2e1 * B * (0.2e1 * t13 + t12 * Rp - 0.3e1 * t2 + t5 + 0.2e1 * t7 + t8 * Rp) * AA - t21) * t51 + ((0.2e1 * t54 * t58 * AA + t67) * t70 + (0.2e1 * t54 * t72 * AA - t67) * t78 + (t81 + 0.2e1 * t54 * t72 * BB) * t86 + (t81 - 0.2e1 * t54 * t58 * BB) * t92) * t96 + ((-0.2e1 * t99 * t100 * AA - t106) * t70 + (-0.2e1 * t99 * t109 * AA + t106) * t78 + (-t116 - 0.2e1 * t99 * t109 * BB) * t86 + (-t116 + 0.2e1 * t99 * t100 * BB) * t92) * t130 + 0.4e1 * t4 * t98 * t53 * AA * t135) / (((-0.8e1 * t4 - 0.8e1 * t1) * t12 * t24 + 0.8e1 * t144 + 0.8e1 * t12 * t1) * t51 + (0.4e1 * t150 - 0.4e1 * t144) * t135 + 0.4e1 * t150 - 0.4e1 * t144);
	
	t1 = Rm * Rp;
	t2 = Rm * Rm;
	t3 = km * km;
	t4 = Rp * Rp;
	t5 = B * B;
	t12 = km * Rm;
	t17 = 0.2e1 * Rm;
	t18 = cos(t17);
	t22 = t2 * Rp;
	t25 = B * t3;
	t26 = t5 * B;
	t33 = t5 * km;
	t38 = sin(t17);
	t40 = Rm * B;
	t41 = 0.3e1 * t5;
	t51 = exp(-0.2e1 * Rp);
	t53 = B + Rp;
	t54 = Rm * t53;
	t57 = t41 + 0.2e1 * Rp * B + t2 - t3 - t4;
	t59 = t54 * t57 * AA;
	t60 = B * km;
	t61 = 0.2e1 * t60;
	t62 = t40 + t61 - t1;
	t67 = Rm - km;
	t68 = cos(t67);
	t70 = -t61 + t40 - t1;
	t75 = Rm + km;
	t76 = cos(t75);
	t82 = t54 * t57 * BB;
	t84 = sin(t75);
	t90 = sin(t67);
	t94 = exp(-t53);
	t97 = 0.3e1 * Rm * t26;
	t98 = t2 * Rm;
	t99 = t98 * B;
	t100 = t3 * Rm;
	t101 = t100 * Rp;
	t103 = Rm * t4 * B;
	t104 = t4 * Rp;
	t105 = Rm * t104;
	t107 = 0.8e1 * t33 * Rp;
	t109 = 0.5e1 * t1 * t5;
	t110 = t98 * Rp;
	t111 = t100 * B;
	t112 = t97 + t99 - t101 + t103 - t105 + t107 + t109 + t110 - t111;
	t114 = t2 * t4;
	t116 = 0.2e1 * t60 * t1;
	t117 = t2 * t5;
	t119 = 0.3e1 * t26 * Rp;
	t120 = t104 * B;
	t121 = t4 * t5;
	t122 = 0.2e1 * t121;
	t123 = t22 * B;
	t125 = 0.2e1 * t33 * Rm;
	t126 = t25 * Rp;
	t127 = t114 + t116 + t117 - t119 + t120 + t122 + t123 + t125 + t126;
	t132 = -t107 + t103 - t105 - t101 + t97 - t111 + t110 + t109 + t99;
	t134 = t120 - t125 + t123 - t116 + t122 + t117 + t114 + t126 - t119;
	t152 = exp(-0.3e1 * Rp - B);
	t161 = exp(-0.4e1 * Rp);
	C2 = (((0.2e1 * t1 * (t2 - t3 - t4 + 0.5e1 * t5) * AA - 0.8e1 * B * BB * t12 * Rp) * t18 + (0.2e1 * Rp * (0.2e1 * t5 * Rp + 0.2e1 * t22 + t2 * B + t25 - 0.3e1 * t26 + B * t4) * AA + 0.8e1 * t33 * BB * Rp) * t38 + 0.2e1 * t40 * (t41 + t4 + t2 - t3) * AA - 0.8e1 * t5 * BB * t12) * t51 + ((-t59 + 0.2e1 * t54 * t62 * BB) * t68 + (-t59 - 0.2e1 * t54 * t70 * BB) * t76 + (0.2e1 * t54 * t70 * AA - t82) * t84 + (0.2e1 * t54 * t62 * AA + t82) * t90) * t94 + ((t112 * AA - 0.2e1 * t127 * BB) * t68 + (t132 * AA + 0.2e1 * t134 * BB) * t76 + (-0.2e1 * t134 * AA + t132 * BB) * t84 + (-0.2e1 * t127 * AA - t112 * BB) * t90) * t152 + (-0.2e1 * t59 + 0.8e1 * t40 * km * t53 * BB) * t161) / (((-0.8e1 * t2 - 0.8e1 * t5) * t4 * t18 + 0.8e1 * t117 + 0.8e1 * t121) * t51 + (0.4e1 * t114 - 0.4e1 * t117) * t161 + 0.4e1 * t114 - 0.4e1 * t117);
	
	t1 = B * B;
	t2 = t1 * Rp;
	t4 = Rm * Rm;
	t5 = t4 * Rp;
	t7 = Rp * Rp;
	t8 = B * t7;
	t11 = km * km;
	t13 = t4 * B;
	t21 = 0.8e1 * t1 * km * BB * Rp;
	t23 = 0.2e1 * Rm;
	t24 = cos(t23);
	t26 = Rm * Rp;
	t38 = sin(t23);
	t51 = exp(-0.2e1 * Rp);
	t53 = B + Rp;
	t54 = Rm * t53;
	t55 = Rm * B;
	t57 = 0.2e1 * B * km;
	t58 = t55 + t57 - t26;
	t62 = 0.3e1 * t1;
	t64 = 0.2e1 * Rp * B;
	t65 = t62 + t64 + t4 - t11 - t7;
	t67 = t54 * t65 * BB;
	t69 = Rm - km;
	t70 = cos(t69);
	t72 = -t57 + t55 - t26;
	t77 = Rm + km;
	t78 = cos(t77);
	t81 = t54 * t65 * AA;
	t86 = sin(t77);
	t92 = sin(t69);
	t96 = exp(-t53);
	t98 = B - Rp;
	t99 = Rm * t98;
	t100 = t55 + t57 + t26;
	t104 = t62 - t64 + t4 - t11 - t7;
	t106 = t99 * t104 * BB;
	t109 = -t57 + t55 + t26;
	t116 = t99 * t104 * AA;
	t130 = exp(-0.3e1 * Rp - B);
	t141 = t4 * t1;
	t147 = t4 * t7;
	t151 = exp(-0.4e1 * Rp);
	C3 = (((-0.2e1 * Rp * (-0.2e1 * t2 - 0.2e1 * t5 + t8 - 0.3e1 * t1 * B + B * t11 + t13) * AA - t21) * t24 + (0.2e1 * t26 * (t4 - t11 - t7 + 0.5e1 * t1) * AA - 0.8e1 * B * BB * km * Rm * Rp) * t38 - 0.2e1 * B * (0.2e1 * t8 + 0.2e1 * t13 + 0.3e1 * t2 - t7 * Rp - t5 - t11 * Rp) * AA + t21) * t51 + ((-0.2e1 * t54 * t58 * AA - t67) * t70 + (-0.2e1 * t54 * t72 * AA + t67) * t78 + (-t81 - 0.2e1 * t54 * t72 * BB) * t86 + (-t81 + 0.2e1 * t54 * t58 * BB) * t92) * t96 + ((0.2e1 * t99 * t100 * AA + t106) * t70 + (0.2e1 * t99 * t109 * AA - t106) * t78 + (t116 + 0.2e1 * t99 * t109 * BB) * t86 + (t116 - 0.2e1 * t99 * t100 * BB) * t92) * t130 + 0.4e1 * t4 * t98 * t53 * AA) / (((-0.8e1 * t4 - 0.8e1 * t1) * t7 * t24 + 0.8e1 * t141 + 0.8e1 * t7 * t1) * t51 + (0.4e1 * t147 - 0.4e1 * t141) * t151 + 0.4e1 * t147 - 0.4e1 * t141);
	
	t1 = Rm * Rp;
	t2 = Rm * Rm;
	t3 = km * km;
	t4 = Rp * Rp;
	t5 = B * B;
	t12 = km * Rm;
	t17 = 0.2e1 * Rm;
	t18 = cos(t17);
	t22 = t2 * Rp;
	t25 = t5 * B;
	t27 = B * t3;
	t33 = t5 * km;
	t38 = sin(t17);
	t40 = Rm * B;
	t41 = 0.3e1 * t5;
	t51 = exp(-0.2e1 * Rp);
	t53 = t2 * Rm;
	t54 = t53 * B;
	t56 = 0.5e1 * t1 * t5;
	t58 = Rm * t4 * B;
	t59 = t3 * Rm;
	t60 = t59 * Rp;
	t62 = 0.8e1 * t33 * Rp;
	t64 = 0.3e1 * Rm * t25;
	t65 = t53 * Rp;
	t66 = t59 * B;
	t67 = t4 * Rp;
	t68 = Rm * t67;
	t69 = t54 - t56 + t58 + t60 - t62 + t64 - t65 - t66 + t68;
	t71 = t2 * t4;
	t73 = 0.3e1 * t25 * Rp;
	t74 = t2 * t5;
	t75 = t27 * Rp;
	t76 = B * km;
	t78 = 0.2e1 * t76 * t1;
	t80 = 0.2e1 * t33 * Rm;
	t81 = t22 * B;
	t82 = t4 * t5;
	t83 = 0.2e1 * t82;
	t84 = t67 * B;
	t85 = t71 + t73 + t74 - t75 - t78 + t80 - t81 + t83 - t84;
	t89 = Rm - km;
	t90 = cos(t89);
	t92 = t60 - t66 - t65 + t58 + t54 - t56 + t62 + t68 + t64;
	t94 = t73 + t78 - t81 + t74 - t80 - t84 - t75 + t83 + t71;
	t98 = Rm + km;
	t99 = cos(t98);
	t105 = sin(t98);
	t111 = sin(t89);
	t115 = exp(-Rp - B);
	t117 = B - Rp;
	t118 = Rm * t117;
	t121 = t41 - 0.2e1 * Rp * B + t2 - t3 - t4;
	t123 = t118 * t121 * AA;
	t124 = 0.2e1 * t76;
	t125 = t40 + t124 + t1;
	t131 = -t124 + t40 + t1;
	t141 = t118 * t121 * BB;
	t152 = exp(-0.3e1 * Rp - B);
	t171 = exp(-0.4e1 * Rp);
	C4 = (((-0.2e1 * t1 * (t2 - t3 - t4 + 0.5e1 * t5) * AA + 0.8e1 * B * BB * t12 * Rp) * t18 + (-0.2e1 * Rp * (-0.2e1 * t5 * Rp - 0.2e1 * t22 + t4 * B - 0.3e1 * t25 + t27 + t2 * B) * AA - 0.8e1 * t33 * BB * Rp) * t38 + 0.2e1 * t40 * (t41 + t4 + t2 - t3) * AA - 0.8e1 * t5 * BB * t12) * t51 + ((t69 * AA - 0.2e1 * t85 * BB) * t90 + (t92 * AA + 0.2e1 * t94 * BB) * t99 + (-0.2e1 * t94 * AA + t92 * BB) * t105 + (-0.2e1 * t85 * AA - t69 * BB) * t111) * t115 + ((-t123 + 0.2e1 * t118 * t125 * BB) * t90 + (-t123 - 0.2e1 * t118 * t131 * BB) * t99 + (0.2e1 * t118 * t131 * AA - t141) * t105 + (0.2e1 * t118 * t125 * AA + t141) * t111) * t152 - 0.2e1 * t123 + 0.8e1 * t40 * km * t117 * BB) / (((-0.8e1 * t2 - 0.8e1 * t5) * t4 * t18 + 0.8e1 * t74 + 0.8e1 * t82) * t51 + (0.4e1 * t71 - 0.4e1 * t74) * t171 + 0.4e1 * t71 - 0.4e1 * t74);
	
	#******************************************************************%
	#******************************************************************%
	
	#*******************************************%
	#*       calculate the velocities etc      *%
	#*******************************************%
	
	t2 = exp(UU * x[2]);
	t3 = Rm * x[2];
	t4 = cos(t3);
	t6 = sin(t3);
	t11 = exp(-VV * x[2]);
	t18 = exp(-0.2e1 * x[2] * B);
	t19 = km * x[2];
	t20 = cos(t19);
	t22 = sin(t19);
	u1 = kn * (t2 * (C1 * t4 + C2 * t6) + t11 * (C3 * t4 + C4 * t6) + t18 * (AA * t20 + BB * t22));
	
	t1 = Rm * x[2];
	t2 = cos(t1);
	t4 = sin(t1);
	t14 = exp(UU * x[2]);
	t26 = exp(-VV * x[2]);
	t28 = km * x[2];
	t29 = cos(t28);
	t31 = sin(t28);
	t43 = exp(-0.2e1 * x[2] * B);
	u2 = (-UU * (C1 * t2 + C2 * t4) + C1 * t4 * Rm - C2 * t2 * Rm) * t14 + (VV * (C3 * t2 + C4 * t4) + C3 * t4 * Rm - C4 * t2 * Rm) * t26 + (0.2e1 * B * (AA * t29 + BB * t31) + AA * t31 * km - BB * t29 * km) * t43;
	
	t2 = 0.2e1 * x[2] * B;
	t3 = exp(t2);
	t4 = t3 * kn;
	t5 = Rm * x[2];
	t6 = cos(t5);
	t8 = sin(t5);
	t18 = exp(UU * x[2]);
	t31 = exp(-VV * x[2]);
	t34 = km * x[2];
	t35 = cos(t34);
	t37 = sin(t34);
	t47 = exp(-t2);
	u3 = 0.2e1 * t4 * (UU * (C1 * t6 + C2 * t8) - C1 * t8 * Rm + C2 * t6 * Rm) * t18 + 0.2e1 * t4 * (-VV * (C3 * t6 + C4 * t8) - C3 * t8 * Rm + C4 * t6 * Rm) * t31 + 0.2e1 * t4 * (-0.2e1 * B * (AA * t35 + BB * t37) - AA * t37 * km + BB * t35 * km) * t47;
	
	t1 = Rm * Rm;
	t3 = UU * UU;
	t8 = kn * kn;
	t11 = Rm * x[2];
	t12 = sin(t11);
	t14 = cos(t11);
	t20 = t14 * Rm;
	t27 = 0.2e1 * x[2] * B;
	t28 = exp(t27);
	t31 = exp(UU * x[2]);
	t38 = VV * VV;
	t54 = exp(-VV * x[2]);
	t56 = km * km;
	t59 = B * B;
	t66 = km * x[2];
	t67 = sin(t66);
	t69 = cos(t66);
	t83 = exp(-t27);
	u4 = ((C2 * t1 - t3 * C2 + 0.2e1 * UU * C1 * Rm - C2 * t8) * t12 + C1 * t14 * t1 - t3 * C1 * t14 - 0.2e1 * UU * C2 * t20 - t8 * C1 * t14) * t28 * t31 + ((-0.2e1 * VV * C3 * Rm + C4 * t1 - C4 * t8 - t38 * C4) * t12 + 0.2e1 * VV * C4 * t20 + C3 * t14 * t1 - t8 * C3 * t14 - t38 * C3 * t14) * t28 * t54 + ((BB * t56 - t8 * BB - 0.4e1 * t59 * BB - 0.4e1 * B * AA * km) * t67 + AA * t69 * t56 - t8 * AA * t69 - 0.4e1 * t59 * AA * t69 + 0.4e1 * B * BB * t69 * km) * t28 * t83;
	
	
	t1 = Rm * x[2];
	t2 = sin(t1);
	t3 = Rm * Rm;
	t4 = t3 * Rm;
	t5 = t2 * t4;
	t6 = UU * UU;
	t7 = t6 * UU;
	t8 = cos(t1);
	t15 = 0.2e1 * B * t8 * t3;
	t19 = B * UU;
	t20 = t2 * Rm;
	t23 = kn * kn;
	t24 = B * t23;
	t26 = 0.2e1 * t24 * t8;
	t27 = t23 * UU;
	t29 = B * t6;
	t33 = t23 * t2 * Rm;
	t35 = 0.1e1 / kn;
	t42 = 0.2e1 * B * t2 * t3;
	t43 = t8 * t4;
	t45 = 0.2e1 * t24 * t2;
	t52 = t23 * t8 * Rm;
	t53 = t8 * Rm;
	t64 = 0.2e1 * x[2] * B;
	t65 = exp(t64);
	t68 = exp(UU * x[2]);
	t70 = B * VV;
	t76 = t23 * VV;
	t78 = VV * VV;
	t79 = t78 * VV;
	t84 = B * t78;
	t108 = exp(-VV * x[2]);
	t111 = km * x[2];
	t112 = sin(t111);
	t113 = km * km;
	t118 = cos(t111);
	t119 = t118 * km;
	t121 = B * B;
	t123 = t112 * km;
	t130 = t113 * km;
	t148 = exp(-t64);
	u5 = (-(-t5 - t7 * t8 + 0.3e1 * UU * t8 * t3 + t15 + 0.3e1 * t6 * t2 * Rm + 0.4e1 * t19 * t20 - t26 + t27 * t8 - 0.2e1 * t29 * t8 - t33) * t35 * C1 - (-t7 * t2 + t27 * t2 + t42 + t43 - t45 + 0.3e1 * UU * t2 * t3 - 0.2e1 * t29 * t2 + t52 - 0.4e1 * t19 * t53 - 0.3e1 * t6 * t8 * Rm) * t35 * C2) * t65 * t68 + (-(t15 - 0.4e1 * t70 * t20 - t33 - 0.3e1 * VV * t8 * t3 - t76 * t8 + t79 * t8 + 0.3e1 * t78 * t2 * Rm - 0.2e1 * t84 * t8 - t26 - t5) * t35 * C3 - (t52 - 0.3e1 * VV * t2 * t3 + t79 * t2 + 0.4e1 * t70 * t53 - 0.3e1 * t78 * t8 * Rm - 0.2e1 * t84 * t2 + t43 - t76 * t2 + t42 - t45) * t35 * C4) * t65 * t108 - t65 * (-0.4e1 * B * BB * t112 * t113 + t23 * BB * t119 + 0.4e1 * t121 * AA * t123 - 0.4e1 * t121 * BB * t119 + BB * t118 * t130 - AA * t112 * t130 - 0.4e1 * B * AA * t118 * t113 - t23 * AA * t123 - 0.4e1 * t24 * AA * t118 - 0.4e1 * t24 * BB * t112) * t35 * t148;
	
	
	
	t2 = 0.2e1 * x[2] * B;
	t3 = exp(t2);
	t4 = t3 * kn;
	t5 = Rm * x[2];
	t6 = cos(t5);
	t8 = sin(t5);
	t18 = exp(UU * x[2]);
	t31 = exp(-VV * x[2]);
	t34 = km * x[2];
	t35 = cos(t34);
	t37 = sin(t34);
	t47 = exp(-t2);
	u6 = -0.2e1 * t4 * (UU * (C1 * t6 + C2 * t8) - C1 * t8 * Rm + C2 * t6 * Rm) * t18 - 0.2e1 * t4 * (-VV * (C3 * t6 + C4 * t8) - C3 * t8 * Rm + C4 * t6 * Rm) * t31 - 0.2e1 * t4 * (-0.2e1 * B * (AA * t35 + BB * t37) - AA * t37 * km + BB * t35 * km) * t47;
		
	#******************************************************************%
	#******************************************************************%
	u1   = u1 * cos(n*π*x[1]); # x[2] velocity 
	sum1 = sum1 + u1;
	u2   = u2 * sin(n*π*x[1]); # x[1] velocity 
	sum2 = sum2 + u2;

	sum5 = sum5 + u5*cos(n*π*x[1]);  # pressure 

	u6   = u6 - u5; # get total stress 
	sum6 = sum6 + u6*cos(n*π*x[1]);  # xx stress 
	u3   = u3 - u5; # get total stress 
	u3   = u3 * cos(n*π*x[1]); # zz stress 
    sum3 = sum3 + u3;
	u4   = u4 * sin(n*π*x[1]); # zx stress 
	sum4 = sum4 + u4;
	
	ρ    = -σ*sin(km*x[2])*cos(kn*x[1]); # density 
	sum7 = sum7 + ρ;
	
	SS  = exp(UU*x[2])*(C1*cos(Rm*x[2])+C2*sin(Rm*x[2])) +exp(-VV*x[2])*(C3*cos(Rm*x[2])+C4*sin(Rm*x[2])) + exp(-2*x[2]*B)*(AA*cos(km*x[2])+BB*sin(km*x[2]));
	SS  = SS * sin(kn*x[1]); # stream function 
	mag = sqrt(u1*u1+u2*u2);
	
	# Output 
    vx  = sum2;
    vz  = sum1;
    p   = sum5;
    ρ   = sum7;
    η   = exp( 2.0 * B * x[2] );
    
    return ( p=p, V=[vx, vz], ρ=ρ, η=η)
end

function Stokes2D_SolKz_Zhong1996(coords::Tuple;
    params = (Δη=1e6, B = 6.9, km = 1.6*π, n = 3, σ = 1.0  ) )
    X = SVector(values(coords)...)
    sol = Stokes2D_SolKz_Zhong1996(X; params)
    return (p=sol.p, V=(x=sol.V[1], y=sol.V[2]), η=sol.η, ρ=sol.ρ) 
end