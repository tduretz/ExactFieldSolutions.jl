function ComputeConstants( params )
    @unpack θ1, θ2, η1, η2, Vr0, Vt0, Vr2, Vt2 = params
    A1 = -Vt0

    sinθ2, cosθ2  = sincos(θ2)
    sinθ1, cosθ1  = sincos(θ1)
    sinθ2_square  = @fastpow sinθ2^2
    coseθ2_square = @fastpow cosθ2^2
    sinθ1_square  = @fastpow sinθ1^2
    cosθ1_square  = @fastpow cosθ1^2
    B1 = @fastpow ((sinθ1_square*(η1*η2*coseθ2_square+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1*η2*sinθ2_square+η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2))*Vr0+(
         sinθ1_square*(η1*η2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square)+θ1*(η1*η2*(coseθ2_square+sinθ2_square)+η2^2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1_square*
         (η1*η2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square)+θ1*(η1*η2*(coseθ2_square+sinθ2_square)+η2^2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1*sinθ1*
         (η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1*η2*(sinθ2_square+coseθ2_square)))*Vt0+cosθ1*sinθ1*(η2^2*(sinθ2*Vr2+cosθ2*Vt2)+η1*η2*(-1*sinθ2*Vr2-1*cosθ2*Vt2))+
         sinθ1_square*
         (θ1*(η1*η2*(-1*cosθ2*Vt2-1*sinθ2*Vr2)+η2^2*(cosθ2*Vt2+sinθ2*Vr2))+η2^2*(sinθ2*Vt2-1*cosθ2*Vr2)+η1*η2*((θ2*sinθ2+cosθ2)*Vr2+θ2*cosθ2*Vt2))
         +cosθ1_square*(θ1*(η1*η2*(-1*cosθ2*Vt2-1*sinθ2*Vr2)+η2^2*(cosθ2*Vt2+sinθ2*Vr2))+η1*η2*(θ2*sinθ2*Vr2+θ2*cosθ2*Vt2+sinθ2*Vt2)))/(sinθ1_square*(2*η1*
         η2*coseθ2_square+η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2));
    C1  = @fastpow  -(1*((sinθ1_square*(η2^2*(sinθ2_square+coseθ2_square)-1*η1*η2*coseθ2_square)+η1*η2*cosθ1_square*sinθ2_square)*Vr0+(sinθ1_square*
         (η1*η2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square)+θ1*(η1*η2*(coseθ2_square+sinθ2_square)+η2^2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1_square*
         (η1*η2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square)+θ1*(η1*η2*(coseθ2_square+sinθ2_square)+η2^2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1*sinθ1*
         (η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1*η2*(sinθ2_square+coseθ2_square)))*Vt0+cosθ1*sinθ1*(η2^2*(sinθ2*Vr2+cosθ2*Vt2)+η1*η2*(-1*sinθ2*Vr2-1*cosθ2*Vt2))+
         sinθ1_square*
         (θ1*(η1*η2*(-1*cosθ2*Vt2-1*sinθ2*Vr2)+η2^2*(cosθ2*Vt2+sinθ2*Vr2))+η2^2*(sinθ2*Vt2-1*cosθ2*Vr2)+η1*η2*((θ2*sinθ2+cosθ2)*Vr2+θ2*cosθ2*Vt2))
         +cosθ1_square*(θ1*(η1*η2*(-1*cosθ2*Vt2-1*sinθ2*Vr2)+η2^2*(cosθ2*Vt2+sinθ2*Vr2))+η1*η2*(θ2*sinθ2*Vr2+θ2*cosθ2*Vt2+sinθ2*Vt2))))/(sinθ1_square*(2*η1*
         η2*coseθ2_square+η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2));
    D1 = @fastpow -(1*((sinθ1_square*(η1*η2*(θ2*sinθ2_square-1*cosθ2*sinθ2+θ2*coseθ2_square)+θ1*(η2^2*(coseθ2_square+sinθ2_square)+η1*η2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1_square*
         (η1*η2*(θ2*sinθ2_square-1*cosθ2*sinθ2+θ2*coseθ2_square)+θ1*(η2^2*(coseθ2_square+sinθ2_square)+η1*η2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1*sinθ1*
         (η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1*η2*(sinθ2_square+coseθ2_square)))*Vr0+(sinθ1_square*(η1*η2*coseθ2_square+η2^2*(-1*sinθ2_square-1*coseθ2_square))-1*η1*η2*cosθ1_square*sinθ2_square)*
         Vt0+cosθ1*sinθ1*(η1*η2*(sinθ2*Vt2-1*cosθ2*Vr2)+η2^2*(cosθ2*Vr2-1*sinθ2*Vt2))+cosθ1_square*
         (θ1*(η1*η2*(cosθ2*Vr2-1*sinθ2*Vt2)+η2^2*(sinθ2*Vt2-1*cosθ2*Vr2))+η1*η2*((sinθ2-1*θ2*cosθ2)*Vr2+θ2*sinθ2*Vt2))+sinθ1_square*(θ1*
         (η1*η2*(cosθ2*Vr2-1*sinθ2*Vt2)+η2^2*(sinθ2*Vt2-1*cosθ2*Vr2))+η1*η2*(-1*θ2*cosθ2*Vr2-1*cosθ2*Vt2+θ2*sinθ2*Vt2)+η2^2*(sinθ2*Vr2+cosθ2*Vt2))
         ))/(sinθ1_square*(2*η1*η2*coseθ2_square+η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2));
    A2 = @fastpow -(1*((sinθ1_square*
         (θ1*(η1*η2*sinθ2_square-1*η1^2*sinθ2_square)+η1*η2*(-1*θ2*sinθ2_square+cosθ2*sinθ2-1*θ2*coseθ2_square)+η1^2*(θ2*sinθ2_square-1*cosθ2*sinθ2+θ2*coseθ2_square))+θ1*
         cosθ1_square*(η1*η2*sinθ2_square-1*η1^2*sinθ2_square)+cosθ1*sinθ1*(η1^2*sinθ2_square-1*η1*η2*sinθ2_square))*Vr0+(sinθ1_square*(-1*η1*η2*sinθ2_square+η1^2*
         (θ2^2*sinθ2_square+θ2^2*coseθ2_square)+θ1*(η1*η2*(θ2*coseθ2_square-1*cosθ2*sinθ2+θ2*sinθ2_square)+η1^2*(-1*θ2*coseθ2_square+cosθ2*sinθ2-1*θ2*sinθ2_square)))+cosθ1_square*
         (η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1*(η1*η2*(θ2*coseθ2_square-1*cosθ2*sinθ2+θ2*sinθ2_square)+η1^2*(-1*θ2*coseθ2_square+cosθ2*sinθ2-1*θ2*sinθ2_square)))+
         cosθ1*sinθ1*(η1^2*(-1*θ2*sinθ2_square+cosθ2*sinθ2-1*θ2*coseθ2_square)+η1*η2*(θ2*sinθ2_square-1*cosθ2*sinθ2+θ2*coseθ2_square)))*Vt0+cosθ1*sinθ1*
         (η1^2*(θ2*sinθ2*Vr2+θ2*cosθ2*Vt2+sinθ2*Vt2)+η1*η2*(-1*θ2*sinθ2*Vr2-1*θ2*cosθ2*Vt2-1*sinθ2*Vt2))+cosθ1_square*(θ1*
         (η1^2*(-1*sinθ2*Vt2-1*θ2*cosθ2*Vt2-1*θ2*sinθ2*Vr2)+η1*η2*(sinθ2*Vt2+θ2*cosθ2*Vt2+θ2*sinθ2*Vr2))+θ1^2*
         (η1*η2*(-2*cosθ2*Vt2-2*sinθ2*Vr2)+η1^2*(cosθ2*Vt2+sinθ2*Vr2)+η2^2*(cosθ2*Vt2+sinθ2*Vr2)))+sinθ1_square*(θ1*
         (η1^2*(-1*sinθ2*Vt2-1*θ2*cosθ2*Vt2-1*θ2*sinθ2*Vr2)+η1*η2*(sinθ2*Vt2+θ2*cosθ2*Vt2+θ2*sinθ2*Vr2))+θ1^2*
         (η1*η2*(-2*cosθ2*Vt2-2*sinθ2*Vr2)+η1^2*(cosθ2*Vt2+sinθ2*Vr2)+η2^2*(cosθ2*Vt2+sinθ2*Vr2))+η1^2*(-1*θ2*cosθ2*Vr2-1*cosθ2*Vt2+θ2*sinθ2*Vt2)+
         η1*η2*((sinθ2+θ2*cosθ2)*Vr2+2*cosθ2*Vt2-1*θ2*sinθ2*Vt2)+η2^2*(-1*sinθ2*Vr2-1*cosθ2*Vt2))))/(sinθ1_square*(2*η1*η2*coseθ2_square+η2^2*
         (-1*sinθ2_square-1*coseθ2_square)+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*(η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+
         θ1*(η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2));
    B2 = @fastpow ((cosθ1_square*
         (η1^2*(θ2^2*sinθ2_square+θ2^2*coseθ2_square)+θ1*(η1*η2*(θ2*coseθ2_square+cosθ2*sinθ2+θ2*sinθ2_square)+η1^2*(-1*θ2*coseθ2_square-1*cosθ2*sinθ2-1*θ2*sinθ2_square)))+sinθ1_square*(
         η1*η2*coseθ2_square+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1*
         (η1*η2*(θ2*coseθ2_square+cosθ2*sinθ2+θ2*sinθ2_square)+η1^2*(-1*θ2*coseθ2_square-1*cosθ2*sinθ2-1*θ2*sinθ2_square)))+cosθ1*sinθ1*
         (η1*η2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square)+η1^2*(θ2*sinθ2_square+cosθ2*sinθ2+θ2*coseθ2_square)))*Vr0+(sinθ1_square*
         (θ1*(η1^2*coseθ2_square-1*η1*η2*coseθ2_square)+η1*η2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square))+cosθ1_square*
         (θ1*(η1^2*coseθ2_square-1*η1*η2*coseθ2_square)+η1^2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square))+cosθ1*sinθ1*(η1^2*coseθ2_square-1*η1*η2*coseθ2_square))*Vt0+
         cosθ1*sinθ1*(η1^2*(-1*θ2*cosθ2*Vr2-1*cosθ2*Vt2+θ2*sinθ2*Vt2)+η1*η2*(θ2*cosθ2*Vr2+cosθ2*Vt2-1*θ2*sinθ2*Vt2))+sinθ1_square*(θ1*
         (η1*η2*(-1*θ2*sinθ2*Vt2+cosθ2*Vt2+θ2*cosθ2*Vr2)+η1^2*(θ2*sinθ2*Vt2-1*cosθ2*Vt2-1*θ2*cosθ2*Vr2))+θ1^2*
         (η1^2*(cosθ2*Vr2-1*sinθ2*Vt2)+η2^2*(cosθ2*Vr2-1*sinθ2*Vt2)+η1*η2*(2*sinθ2*Vt2-2*cosθ2*Vr2))+η2^2*(sinθ2*Vt2-1*cosθ2*Vr2)+η1*η2*
         ((θ2*sinθ2+cosθ2)*Vr2+θ2*cosθ2*Vt2))+cosθ1_square*(θ1*
         (η1*η2*(-1*θ2*sinθ2*Vt2+cosθ2*Vt2+θ2*cosθ2*Vr2)+η1^2*(θ2*sinθ2*Vt2-1*cosθ2*Vt2-1*θ2*cosθ2*Vr2))+θ1^2*
         (η1^2*(cosθ2*Vr2-1*sinθ2*Vt2)+η2^2*(cosθ2*Vr2-1*sinθ2*Vt2)+η1*η2*(2*sinθ2*Vt2-2*cosθ2*Vr2))+η1^2*(θ2*sinθ2*Vr2+θ2*cosθ2*Vt2+sinθ2*Vt2)))/(
         sinθ1_square*(2*η1*η2*coseθ2_square+η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2));
    C2 = @fastpow  -(1*((sinθ1_square*(η1*η2*(sinθ2_square+coseθ2_square)-1*η1^2*coseθ2_square)+η1^2*cosθ1_square*sinθ2_square)*Vr0+(sinθ1_square*
         (η1^2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square)+θ1*(η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1_square*
         (η1^2*(-1*θ2*sinθ2_square-1*cosθ2*sinθ2-1*θ2*coseθ2_square)+θ1*(η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1*sinθ1*
         (η1*η2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(sinθ2_square+coseθ2_square)))*Vt0+cosθ1*sinθ1*(η1*η2*(sinθ2*Vr2+cosθ2*Vt2)+η1^2*(-1*sinθ2*Vr2-1*cosθ2*Vt2))+
         sinθ1_square*
         (θ1*(η1^2*(-1*cosθ2*Vt2-1*sinθ2*Vr2)+η1*η2*(cosθ2*Vt2+sinθ2*Vr2))+η1*η2*(sinθ2*Vt2-1*cosθ2*Vr2)+η1^2*((θ2*sinθ2+cosθ2)*Vr2+θ2*cosθ2*Vt2))
         +cosθ1_square*(θ1*(η1^2*(-1*cosθ2*Vt2-1*sinθ2*Vr2)+η1*η2*(cosθ2*Vt2+sinθ2*Vr2))+η1^2*(θ2*sinθ2*Vr2+θ2*cosθ2*Vt2+sinθ2*Vt2))))/(sinθ1_square*(2*η1*η2*
         coseθ2_square+η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2));
    D2 = @fastpow -(1*((sinθ1_square*(η1^2*(θ2*sinθ2_square-1*cosθ2*sinθ2+θ2*coseθ2_square)+θ1*(η1*η2*(coseθ2_square+sinθ2_square)+η1^2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1_square*
         (η1^2*(θ2*sinθ2_square-1*cosθ2*sinθ2+θ2*coseθ2_square)+θ1*(η1*η2*(coseθ2_square+sinθ2_square)+η1^2*(-1*coseθ2_square-1*sinθ2_square)))+cosθ1*sinθ1*
         (η1*η2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(sinθ2_square+coseθ2_square)))*Vr0+(sinθ1_square*(η1^2*coseθ2_square+η1*η2*(-1*sinθ2_square-1*coseθ2_square))-1*η1^2*cosθ1_square*sinθ2_square)*Vt0+
         cosθ1*sinθ1*(η1^2*(sinθ2*Vt2-1*cosθ2*Vr2)+η1*η2*(cosθ2*Vr2-1*sinθ2*Vt2))+cosθ1_square*
         (θ1*(η1^2*(cosθ2*Vr2-1*sinθ2*Vt2)+η1*η2*(sinθ2*Vt2-1*cosθ2*Vr2))+η1^2*((sinθ2-1*θ2*cosθ2)*Vr2+θ2*sinθ2*Vt2))+sinθ1_square*(θ1*
         (η1^2*(cosθ2*Vr2-1*sinθ2*Vt2)+η1*η2*(sinθ2*Vt2-1*cosθ2*Vr2))+η1^2*(-1*θ2*cosθ2*Vr2-1*cosθ2*Vt2+θ2*sinθ2*Vt2)+η1*η2*(sinθ2*Vr2+cosθ2*Vt2))
         ))/(sinθ1_square*(2*η1*η2*coseθ2_square+η2^2*(-1*sinθ2_square-1*coseθ2_square)+η1^2*(θ2^2*sinθ2_square+(θ2^2-1)*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1_square*(η1^2*((θ2^2-1)*sinθ2_square+θ2^2*coseθ2_square)+θ1^2*
         (η2^2*(coseθ2_square+sinθ2_square)+η1^2*(coseθ2_square+sinθ2_square)+η1*η2*(-2*coseθ2_square-2*sinθ2_square))+θ1*
         (η1*η2*(2*θ2*coseθ2_square+2*θ2*sinθ2_square)+η1^2*(-2*θ2*coseθ2_square-2*θ2*sinθ2_square)))+cosθ1*sinθ1*(2*η1^2*cosθ2*sinθ2-2*η1*η2*cosθ2*sinθ2));
    return (A1, B1, C1, D1, A2, B2, C2, D2)
end

function ComputeLocation( X, c, params )
    @unpack θ1, θ2, η1, η2, Vr0, Vt0, Vr2, Vt2 = params
    r = sqrt(X[1]^2 + X[2]^2)
    θ = atan(X[2]/X[1])
    # Corner facing up or
    if sign(θ1) == -1
        if (θ ≥  0) θ -= π end
    else
        if (θ ≤ 0) θ += π end
    end
    # Define which corner is which
    ind1, ind2 = false, false
    if sign(θ1) == -1
        # if (θ≥θ1)         ind1   = true end
        # if (θ<θ1 && θ>θ2) ind2   = true end
        ind1 = θ ≥ θ1
        ind2 = θ < θ1 && θ > θ2
    else
        # if (θ≤θ1)         ind1   = true end
        # if (θ>θ1 && θ<θ2) ind2   = true end
        ind1 = θ ≤ θ1
        ind2 = θ > θ1 && θ < θ2
    end
    η, A, B, C, D = NaN, NaN, NaN, NaN, NaN
    if ind1
        η = η1
        A = c[1]
        B = c[2]
        C = c[3]
        D = c[4]
    elseif ind2
        η = η2
        A = c[5]
        B = c[6]
        C = c[7]
        D = c[8]
    end
    return r, θ, η, A, B, C, D, ind1, ind2
end

function Stokes2D_Moulas2021_p( X, c, params ) 
    @unpack θ1, θ2, η1, η2, Vr0, Vt0, Vr2, Vt2 = params
    # Compute location
    r, θ, η, A, B, C, D, ind1, ind2 = ComputeLocation( X, c, params )
    # Pressure in both corners
    p   =  -2*η*(-cos(θ)*C - sin(θ)*D)/r
    return p
end

function Stokes2D_Moulas2021_V( X, c, params ) 
    @unpack θ1, θ2, η1, η2, Vr0, Vt0, Vr2, Vt2 = params
    # Compute location
    r, θ, η, A, B, C, D, ind1, ind2 = ComputeLocation( X, c, params )
    # Velocity in both corners
    sinθ, cosθ = sincos(θ)
    Vr  = -A*sinθ + B*cosθ + C * cosθ - C * θ * sinθ + D *sinθ +D * θ * cosθ
    Vt  = -(A + C*θ)*cosθ - (D*θ + B)*sinθ
    𝑉   = @SVector([Vr; Vt])
    R   = @SMatrix([cosθ -sinθ; sinθ cosθ])
    return V = R*𝑉
end

function Stokes2D_Moulas2021_σrt( X, c, params ) 
    @unpack θ1, θ2, η1, η2, Vr0, Vt0, Vr2, Vt2 = params
    r, θ, η, A, B, C, D, ind1, ind2 = ComputeLocation( X, c, params )

    sinθ, cosθ = sincos(θ)
    σrt = η*((cosθ*A+sinθ*B+θ*cosθ*C+θ*sinθ*D)/r+
        (-cosθ*A-sinθ*B-θ*cosθ*C-2*sinθ*C+2*cosθ*D-θ*sinθ*D)/r);
    return σrt
end

@doc raw"""
    sol = Stokes2D_Moulas2021(x; params)  

Evaluates the analytical solution of [Moulas et al. (2021)](https://academic.oup.com/gji/article/227/1/576/6309899):

    x      : is the coordinate vector or tuple 
    params : optional parameter array, default = params = (θ1  =  30π/180, θ2  = 180π/180, η1  = 1e0, η2  = 1e3, Vr0 = -1., Vt0 =  0., Vr2 =  0., Vt2 =  0.,) 
and returns:

    sol    : tuple containing the solution fields p (pressure), V (velocity vector), L (velocity gratdient tensor), ε̇ (deviatoric strain rate tensor) and τ (deviatoric stress tensor)

# Examples
```julia-repl
julia> Stokes2D_Moulas2021( [1, 1] )
(p = 12.53061977287755, V = [-0.010899335607440476, 0.02773892842336736], L = [0.0028483768768197255 -0.0028483768768197255; 0.0028483768768197246 -0.0028483768768197246], ε̇ = [0.0028483768768197255 -4.336808689942018e-19; -4.336808689942018e-19 -0.0028483768768197246], τ = [5.696753753639451 -8.673617379884035e-16; -8.673617379884035e-16 -5.696753753639449], η = 1000.0)
```
```julia-repl
julia> Stokes2D_Moulas2021( (1, 1) )
(p = 12.53061977287755, V = (x = -0.010899335607440476, y = 0.02773892842336736), L = (xx = 0.0028483768768197255, xy = -0.0028483768768197255, yx = 0.0028483768768197246, yy = -0.0028483768768197246), ε̇ = (xx = 0.0028483768768197255, xy = -4.336808689942018e-19, yx = -4.336808689942018e-19, yy = -0.0028483768768197246), τ = (xx = 5.696753753639451, xy = -8.673617379884035e-16, yx = -8.673617379884035e-16, yy = -5.696753753639449), η = 1000.0)
```
"""
function Stokes2D_Moulas2021(X;
    params = (
        θ1  =  30π/180, 
        θ2  = 180π/180, 
        η1  = 1e0, 
        η2  = 1e3, 
        Vr0 = -1e0, 
        Vt0 =  0e0, 
        Vr2 =  0e0, 
        Vt2 =  0e0,) )
    # Constants
    c = ComputeConstants( params )
    # Check position
    r, θ, η, ind1, ind2 = ComputeLocation( X, c, params)
    p = Stokes2D_Moulas2021_p(X, c, params)
    v = Stokes2D_Moulas2021_V(X, c, params)
    # L = Stokes2D_Schmid2003_L(x, params) # just to check if automatic derivatives are correct
    f_cl = x -> Stokes2D_Moulas2021_V(x, c, params)
    L = ForwardDiff.jacobian(f_cl, X)
    # Postprocess deviatoric stress and strain rate
    ε̇ = 1/2*(L + L')
    τ = 2*η*ε̇
    # # CHECK: normal deviatoric stress components should vanish in polar coordinates
    # # Rotation matrix
    # R = @SMatrix([cos(θ) -sin(θ); sin(θ) cos(θ)])
    # τ_rot = R'*τ*R
    # # CHECK: shear stress in polar coordinates should equate the analytical solution
    # τ =  Stokes2D_Moulas2021_σrt( X, c, params ) 
    # τrt[i,j] = τ_rot[1,2]
    # @show τ_rot[1,1], abs(τ_rot[1,2] - τ), ε̇[1,1]+ε̇[2,2]
    return (p=p, V=v, L=L, ε̇=ε̇, τ=τ, η=η)
end

function Stokes2D_Moulas2021(coords::Union{Tuple, NamedTuple};
    params = (
        θ1  =  30π/180, 
        θ2  = 180π/180, 
        η1  = 1e0, 
        η2  = 1e3, 
        Vr0 = -1., 
        Vt0 =  0., 
        Vr2 =  0., 
        Vt2 =  0.,)  )
    X = SVector(values(coords)...)
    sol = Stokes2D_Moulas2021(X; params)
    return (p=sol.p, 
    V=(x=sol.V[1], y=sol.V[2]),
    L=(xx=sol.L[1,1], xy=sol.L[1,2], yx=sol.L[2,1], yy=sol.L[2,2]), 
    ε̇=(xx=sol.ε̇[1,1], xy=sol.ε̇[1,2], yx=sol.ε̇[2,1], yy=sol.ε̇[2,2]), 
    τ=(xx=sol.τ[1,1], xy=sol.τ[1,2], yx=sol.τ[2,1], yy=sol.τ[2,2]),
    η=sol.η) 
end